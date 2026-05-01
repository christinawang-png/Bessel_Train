import os
from pathlib import Path

import numpy as np
from numpy.fft import fft2, ifft2
from scipy.io import loadmat
import nibabel as nib

# ---------- CONFIG ----------
IN_VOL_DIR    = Path("MiniVess/post")  # 3D .nii.gz
OUT_IDEAL_DIR = Path("MiniVess/out_ideal")
OUT_NOISY_DIR = Path("MiniVess/out")

PSF_MAT_PATH  = Path("Example_PSF/64_33/psf_all_1.mat")  # or .m if it's actually a MAT-file
PSF_KEY       = None  # or "psf_thisAngle" etc., if you know the variable name

# Noise / artifact parameters
JITTER_STD_PX   = 0.3    # per-line lateral jitter in pixels (0 to disable)
GAIN_STD        = 0.03   # per-frame multiplicative gain std (0 to disable)
POISSON_SCALE   = 500.0  # global intensity scale for Poisson
READ_NOISE_STD  = 2.0    # Gaussian read noise in photon units
# ----------------------------


def load_psf_from_mat_single(mat_path: Path, key=None) -> np.ndarray:
    """
    Load a single 3D PSF from a .mat (or .m) file.
    Returns psf as float32 with shape (Z, Y, X), normalized by max value.
    """
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

    # If needed, transpose here so that psf[z, y, x] is a z-slice:
    # e.g. if loaded as (Y, X, Z), do psf = np.transpose(psf, (2, 0, 1))
    # Adjust this once by printing psf.shape and checking.
    # Example (commented out):
    psf = np.transpose(psf, (2, 1, 0))

    # Normalize by maximum (like your MATLAB 'norm' option)
    psf_max = psf.max()
    if psf_max > 0:
        psf /= psf_max

    print(f"[PSF] Final PSF shape (Z, Y, X): {psf.shape}, max={psf.max():.4g}")
    return psf


def match_z_to_psf(volume_zyx: np.ndarray, psf_zyx: np.ndarray) -> np.ndarray:
    """
    Clip or pad the volume along Z so that its depth matches the PSF depth,
    following the spirit of your MATLAB code.
    volume_zyx, psf_zyx both shape (Z, Y, X).
    """
    Zv, H, W = volume_zyx.shape
    Zp, _, _ = psf_zyx.shape

    if Zv == Zp:
        return volume_zyx

    if Zv > Zp:
        # Clip symmetrically
        diff = Zv - Zp
        if diff % 2 != 0:
            # If odd, drop one slice to make it even (similar to your odd/even handling)
            Zv -= 1
            volume_zyx = volume_zyx[:Zv, :, :]
            diff = Zv - Zp
        clip = diff // 2
        print(f"  [Z] Clipping volume Z from {Zv + diff} to {Zp}, removing {clip} slices each side")
        return volume_zyx[clip:clip + Zp, :, :]

    else:
        # Zv < Zp: pad with zeros symmetrically
        diff = Zp - Zv
        if diff % 2 != 0:
            # make diff even
            diff -= 1
        pad = diff // 2
        print(f"  [Z] Padding volume Z from {Zv} to {Zp}, adding {pad} slices each side")
        pad_before = np.zeros((pad, H, W), dtype=volume_zyx.dtype)
        pad_after = np.zeros((Zp - Zv - pad, H, W), dtype=volume_zyx.dtype)
        return np.concatenate([pad_before, volume_zyx, pad_after], axis=0)


def forward_proj_rl(psf_zyx: np.ndarray, volume_zyx: np.ndarray) -> np.ndarray:
    """
    Python equivalent of forwardProj_RL_GPU for one PSF and one volume.
    psf_zyx, volume_zyx: shape (Z, Y, X).
    Returns 2D projection (Y, X).
    """
    psf_zyx = np.asarray(psf_zyx, dtype=np.float32)
    volume_zyx = np.asarray(volume_zyx, dtype=np.float32)

    Za, Ha, Wa = volume_zyx.shape
    Zb, Hb, Wb = psf_zyx.shape

    if Za != Zb:
        raise ValueError(f"Z mismatch after matching: vol Z={Za}, psf Z={Zb}")

    r = Ha + Hb
    c = Wa + Wb
    p1 = (r - Ha) // 2
    p2 = (c - Wa) // 2

    a1 = np.zeros((r, c), dtype=np.float32)
    b1 = np.zeros((r, c), dtype=np.float32)
    projection = np.zeros((Ha, Wa), dtype=np.float32)

    for z in range(Za):
        # top-left placement (like MATLAB a1(1:ra,1:ca) = Xguess(:,:,z))
        a1[:Ha, :Wa] = volume_zyx[z, :, :]
        b1[:Hb, :Wb] = psf_zyx[z, :, :]

        con1 = ifft2(fft2(a1) * fft2(b1))
        con1_real = con1.real.astype(np.float32)

        # Compute crop indices so that crop size is exactly (Ha, Wa)
        row_start = (r - Ha) // 2
        col_start = (c - Wa) // 2
        row_end = row_start + Ha
        col_end = col_start + Wa
        
        projection += con1_real[row_start:row_end, col_start:col_end]

    return projection


def add_noise_and_artifacts(
    img2d: np.ndarray,
    jitter_std_px: float = JITTER_STD_PX,
    gain_std: float = GAIN_STD,
    poisson_scale: float = POISSON_SCALE,
    read_noise_std: float = READ_NOISE_STD,
) -> np.ndarray:
    """
    Take an ideal 2D image and make a noisy, device-like version.
    Returns noisy_img (same shape as img2d).
    """
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

    # Load a single 3D PSF and normalize by max
    psf_zyx = load_psf_from_mat_single(PSF_MAT_PATH, key=PSF_KEY)

    nii_exts = (".nii", ".nii.gz")
    for fname in os.listdir(IN_VOL_DIR):
        if not fname.lower().endswith(nii_exts):
            continue

        in_path = IN_VOL_DIR / fname
        print(f"[VOLUME] {in_path}")
        img = nib.load(str(in_path))
        vol = img.get_fdata()  # often (Y, X, Z) or (X, Y, Z), depends on how it was saved

        # Decide how to get (Z, Y, X). If vol.shape looks like (H, W, Z):
        if vol.ndim == 3 and vol.shape[2] != img.shape[0]:
            # likely (H, W, Z) -> (Z, H, W)
            volume_zyx = np.transpose(vol, (2, 0, 1))
        elif vol.ndim == 3:
            # If it's already (Z, H, W), then:
            volume_zyx = vol
        else:
            raise ValueError(f"Unexpected volume shape {vol.shape} for {in_path}")

        # Match Z dimension to PSF
        volume_zyx = match_z_to_psf(volume_zyx, psf_zyx)

        # Forward projection: ideal 2D
        ideal2d = forward_proj_rl(psf_zyx, volume_zyx)

        # Noisy version
        noisy2d = add_noise_and_artifacts(ideal2d)

        # Save as 2D NIfTI (Y, X) with identity affine
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