import os
from pathlib import Path

import numpy as np
from numpy.fft import fft2, ifft2
from scipy.io import loadmat
from scipy.ndimage import map_coordinates
import nibabel as nib
import matplotlib.pyplot as plt

# ---------- CONFIG ----------
IN_VOL_DIR    = Path("MiniVess/post")      # 3D .nii.gz
OUT_IDEAL_DIR = Path("MiniVess/out_ideal_2")
OUT_NOISY_DIR = Path("MiniVess/out_2")

PSF_MAT_PATH  = Path("Example_PSF/64_33/psf_all_1.mat")
PSF_KEY       = None                       # or "psf_thisAngle" if known

# Warp controls (static, same for all samples)
APPLY_WARP_TO_IDEAL = True    # apply warp to ideal image
APPLY_WARP_TO_NOISY = True    # apply warp also to noisy image
WARP_K              = 0.2    # radial distortion coefficient (sign/size controls strength)

# Noise / artifact parameters
ADD_NOISE      = True
JITTER_STD_PX  = 0.3    # per-line lateral jitter in pixels (0 to disable)
GAIN_STD       = 0.03   # per-frame multiplicative gain std (0 to disable)
POISSON_SCALE  = 500.0  # global intensity scale for Poisson
READ_NOISE_STD = 2.0    # Gaussian read noise in photon units
# ----------------------------


def load_psf_from_mat_single(mat_path: Path, key=None) -> np.ndarray:
    mat = loadmat(mat_path)
    if key is not None:
        psf = mat[key]
    else:
        psf = None
        for k, v in mat.items():
            if k.startswith("__"):
                continue
            arr = np.asarray(v)
            if arr.ndim == 3:
                psf = arr
                print(f"[PSF] Using variable '{k}' with shape {arr.shape}")
                break
        if psf is None:
            raise ValueError("No 3D array found in MAT file; specify PSF_KEY.")

    psf = np.asarray(psf, dtype=np.float32)
    # loaded as (X, Y, Z) -> (Z, Y, X)
    psf = np.transpose(psf, (2, 1, 0))
    psf_max = psf.max()
    if psf_max > 0:
        psf /= psf_max
    print(f"[PSF] Final PSF shape (Z, Y, X): {psf.shape}, max={psf.max():.4g}")
    return psf


def match_z_to_psf(volume_zyx: np.ndarray, psf_zyx: np.ndarray) -> np.ndarray:
    Zv, H, W = volume_zyx.shape
    Zp, _, _ = psf_zyx.shape

    if Zv == Zp:
        return volume_zyx

    if Zv > Zp:
        diff = Zv - Zp
        if diff % 2 != 0:
            Zv -= 1
            volume_zyx = volume_zyx[:Zv, :, :]
            diff = Zv - Zp
        clip = diff // 2
        print(f"  [Z] Clipping volume Z from {Zv + diff} to {Zp}, removing {clip} slices each side")
        return volume_zyx[clip:clip + Zp, :, :]
    else:
        diff = Zp - Zv
        if diff % 2 != 0:
            diff -= 1
        pad = diff // 2
        print(f"  [Z] Padding volume Z from {Zv} to {Zp}, adding {pad} slices each side")
        pad_before = np.zeros((pad, H, W), dtype=volume_zyx.dtype)
        pad_after = np.zeros((Zp - Zv - pad, H, W), dtype=volume_zyx.dtype)
        return np.concatenate([pad_before, volume_zyx, pad_after], axis=0)


def forward_proj_rl(psf_zyx: np.ndarray, volume_zyx: np.ndarray,
                    pad_xy: int = 16) -> np.ndarray:
    """
    FFT-based forward projector with reflective boundary handling in XY.

    Strategy:
      - Reflect-pad both volume slice and PSF slice in XY by pad_xy.
      - Do the same FFT-based convolution on the padded grids.
      - Crop back to the original (Ha, Wa) region.

    psf_zyx, volume_zyx: (Z, Y, X)
    pad_xy: how many pixels to reflect-pad on each side in X,Y (tradeoff: quality vs speed).
    """
    psf_zyx = np.asarray(psf_zyx, dtype=np.float32)
    volume_zyx = np.asarray(volume_zyx, dtype=np.float32)

    Za, Ha, Wa = volume_zyx.shape
    Zb, Hb, Wb = psf_zyx.shape
    if Za != Zb:
        raise ValueError(f"Z mismatch after matching: vol Z={Za}, psf Z={Zb}")

    # Reflect-pad slices in XY
    Ha_pad = Ha + 2 * pad_xy
    Wa_pad = Wa + 2 * pad_xy

    projection = np.zeros((Ha, Wa), dtype=np.float32)

    for z in range(Za):
        vol_slice = np.pad(
            volume_zyx[z, :, :],
            ((pad_xy, pad_xy), (pad_xy, pad_xy)),
            mode="reflect",
        )
        psf_slice = np.pad(
            psf_zyx[z, :, :],
            ((pad_xy, pad_xy), (pad_xy, pad_xy)),
            mode="reflect",
        )

        # sizes for FFT-based linear convolution on padded grid
        rb = Ha_pad + Hb
        cb = Wa_pad + Wb
        a1 = np.zeros((rb, cb), dtype=np.float32)
        b1 = np.zeros((rb, cb), dtype=np.float32)

        a1[:Ha_pad, :Wa_pad] = vol_slice
        b1[:Hb + 2 * pad_xy, :Wb + 2 * pad_xy] = psf_slice

        con1 = ifft2(fft2(a1) * fft2(b1))
        con1_real = con1.real.astype(np.float32)

        # crop center (Ha_pad, Wa_pad), then crop away the reflected border
        row_start = (rb - Ha_pad) // 2
        col_start = (cb - Wa_pad) // 2
        row_end = row_start + Ha_pad
        col_end = col_start + Wa_pad

        conv_pad = con1_real[row_start:row_end, col_start:col_end]
        projection += conv_pad[pad_xy:pad_xy + Ha, pad_xy:pad_xy + Wa]

    return projection


def apply_static_warp(img2d: np.ndarray, k: float = WARP_K) -> np.ndarray:
    if k == 0.0:
        return img2d

    img = np.asarray(img2d, dtype=np.float32)
    H, W = img.shape

    y, x = np.meshgrid(
        np.linspace(-1, 1, H, dtype=np.float32),
        np.linspace(-1, 1, W, dtype=np.float32),
        indexing="ij",
    )
    r = np.sqrt(x**2 + y**2)

    # only warp inside r<=0.8, fade to 0 between 0.8 and 1.0
    r0, r1 = 0.8, 1.0
    w = np.clip((r1 - r) / (r1 - r0), 0.0, 1.0)  # 1 in center, 0 at edges
    factor = 1.0 + k * (r**2) * w

    x_dist = x * factor
    y_dist = y * factor

    x_pix = ((x_dist + 1) * 0.5) * (W - 1)
    y_pix = ((y_dist + 1) * 0.5) * (H - 1)

    coords = np.vstack([y_pix.ravel(), x_pix.ravel()])
    warped = map_coordinates(img, coords, order=1, mode="reflect").reshape(H, W)
    return warped


def add_noise_and_artifacts(
    img2d: np.ndarray,
    jitter_std_px: float = JITTER_STD_PX,
    gain_std: float = GAIN_STD,
    poisson_scale: float = POISSON_SCALE,
    read_noise_std: float = READ_NOISE_STD,
    enabled: bool = True,
) -> np.ndarray:
    if not enabled:
        return np.asarray(img2d, dtype=np.float32)

    ideal = np.asarray(img2d, dtype=np.float32)
    H, W = ideal.shape

    # per-frame gain
    g = np.random.normal(loc=1.0, scale=gain_std) if gain_std > 0 else 1.0
    img = g * ideal

    # per-line jitter
    if jitter_std_px > 0:
        jittered = np.zeros_like(img)
        for y in range(H):
            dx = np.random.normal(loc=0.0, scale=jitter_std_px)
            x_coords = np.arange(W, dtype=np.float32) + dx
            x0 = np.floor(x_coords).astype(int)
            x1 = x0 + 1
            w1 = x_coords - x0
            w0 = 1.0 - w1

            x0 = np.clip(x0, 0, W - 1)
            x1 = np.clip(x1, 0, W - 1)
            line = img[y]
            jittered[y] = w0 * line[x0] + w1 * line[x1]
        img = jittered

    # Poisson + Gaussian noise
    scaled = img * poisson_scale
    scaled[scaled < 0] = 0
    noisy_counts = np.random.poisson(scaled).astype(np.float32)

    if read_noise_std > 0:
        noisy_counts += np.random.normal(
            loc=0.0, scale=read_noise_std, size=noisy_counts.shape
        ).astype(np.float32)

    noisy = noisy_counts / poisson_scale
    return noisy


def main():
    OUT_IDEAL_DIR.mkdir(parents=True, exist_ok=True)
    OUT_NOISY_DIR.mkdir(parents=True, exist_ok=True)

    psf_zyx = load_psf_from_mat_single(PSF_MAT_PATH, key=PSF_KEY)

    nii_exts = (".nii", ".nii.gz")
    
    H, W = 256, 256
    step = 16  # size of grid cells
    yy, xx = np.indices((H, W))
    grid = ((xx % step == 0) | (yy % step == 0)).astype(np.float32)
    
    k = 0.2
    warped = apply_static_warp(grid, k=k)
    
    plt.figure(figsize=(6, 3))
    plt.subplot(1, 2, 1); plt.imshow(grid, cmap="gray"); plt.title("Original"); plt.axis("off")
    plt.subplot(1, 2, 2); plt.imshow(warped, cmap="gray"); plt.title(f"Warped k={k}"); plt.axis("off")
    plt.tight_layout()
    plt.savefig("warp_grid_k_{:.3f}.png".format(k), dpi=150)
    plt.close()
    
    for fname in os.listdir(IN_VOL_DIR):
        if not fname.lower().endswith(nii_exts):
            continue

        in_path = IN_VOL_DIR / fname
        print(f"[VOLUME] {in_path}")
        img = nib.load(str(in_path))
        vol = img.get_fdata().astype(np.float32)

        # Assume (H, W, Z) -> (Z, H, W)
        if vol.ndim == 3 and vol.shape[2] != img.shape[0]:
            volume_zyx = np.transpose(vol, (2, 0, 1))
        elif vol.ndim == 3:
            volume_zyx = vol
        else:
            raise ValueError(f"Unexpected volume shape {vol.shape} for {in_path}")

        volume_zyx = match_z_to_psf(volume_zyx, psf_zyx)

        # Ideal forward
        ideal2d = forward_proj_rl(psf_zyx, volume_zyx)

        # Apply static warp to ideal if desired
        if APPLY_WARP_TO_IDEAL:
            ideal2d = apply_static_warp(ideal2d, k=WARP_K)

        # Noisy version (optionally warped again or just noise on warped ideal)
        noisy2d = add_noise_and_artifacts(
            ideal2d if APPLY_WARP_TO_NOISY else forward_proj_rl(psf_zyx, volume_zyx),
            enabled=ADD_NOISE,
        )

        stem = fname
        if stem.endswith(".nii.gz"):
            stem = stem[:-7]
        elif stem.endswith(".nii"):
            stem = stem[:-4]

        out_ideal_path = OUT_IDEAL_DIR / f"{stem}_ideal2d.nii.gz"
        out_noisy_path = OUT_NOISY_DIR / f"{stem}_noisy2d.nii.gz"

        nib.save(nib.Nifti1Image(ideal2d.astype(np.float32), affine=np.eye(4)),
                 str(out_ideal_path))
        nib.save(nib.Nifti1Image(noisy2d.astype(np.float32), affine=np.eye(4)),
                 str(out_noisy_path))

        print(f"  -> ideal: {out_ideal_path.name}, noisy: {out_noisy_path.name}")


if __name__ == "__main__":
    main()