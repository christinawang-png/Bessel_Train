import os
from pathlib import Path
import nibabel as nib
import numpy as np

# --------- CONFIG ---------
ORIG_ROOT = Path("MiniVess/raw/")     # original data (read-only)
PROC_ROOT = Path("MiniVess/post/")    # where processed volumes go
TARGET_Z = 33                                # desired z-depth per chunk
MIN_Z = TARGET_Z                                    # drop volumes thinner than this
# --------------------------


def process_volume(in_path: Path, orig_root: Path, proc_root: Path):
    img = nib.load(str(in_path))
    data = img.get_fdata()        # shape (X, Y, Z) or (X, Y, Z, T)
    affine = img.affine
    header = img.header

    if data.ndim == 4:
        # assume (X, Y, Z, T), take first timepoint for simplicity
        data = data[..., 0]

    X, Y, Z = data.shape

    if Z < MIN_Z:
        print(f"[SKIP] {in_path} has z={Z} < MIN_Z={MIN_Z}")
        return
    
    rel_path = in_path.relative_to(orig_root)
    out_dir = (proc_root / rel_path.parent).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = in_path.stem
    if stem.endswith(".nii"):
        stem = stem[:-4]  # handle .nii.gz stem quir

    if Z == TARGET_Z:
        # Already the right depth, just save as-is
        out_name = f"{stem}_z{TARGET_Z}.nii.gz"
        out_path = out_dir / out_name
        nib.save(nib.Nifti1Image(data, affine, header), str(out_path))
        print(f"[SAVE] {out_path} with original z={Z}")
        
    elif MIN_Z < Z < 2 * TARGET_Z:
        # Between MIN_Z and 2*MIN_Z: take a single centered chunk
        z_start = (Z - TARGET_Z) // 2
        z_end = z_start + TARGET_Z
        chunk = data[:, :, z_start:z_end]
        assert chunk.shape[2] == TARGET_Z
    
        out_name = f"{stem}_z{TARGET_Z}_chunk0.nii.gz"
        out_path = out_dir / out_name
        nib.save(nib.Nifti1Image(chunk, affine, header), str(out_path))
        print(f"[SAVE] {out_path} from z={Z} (center slices {z_start}:{z_end})")
    
    else:
        # Z >= 2*TARGET_Z: make multiple non-overlapping chunks of depth TARGET_Z
        n_chunks = Z // TARGET_Z  # number of full chunks you can fit
        for chunk_idx in range(n_chunks):
            z_start = chunk_idx * TARGET_Z
            z_end = z_start + TARGET_Z
            if z_end > Z:
                break  # safety, though it shouldn't happen
    
            chunk = data[:, :, z_start:z_end]
            assert chunk.shape[2] == TARGET_Z
    
            out_name = f"{stem}_z{TARGET_Z}_chunk{chunk_idx}.nii.gz"
            out_path = out_dir / out_name
            nib.save(nib.Nifti1Image(chunk, affine, header), str(out_path))
            print(f"[SAVE] {out_path} from z={Z} (slices {z_start}:{z_end})")


def main():
    nii_exts = (".nii", ".nii.gz")

    for root, dirs, files in os.walk(ORIG_ROOT):
        root_path = Path(root)
        for fname in files:
            if fname.lower().endswith(nii_exts):
                in_path = root_path / fname
                process_volume(in_path, ORIG_ROOT, PROC_ROOT)


if __name__ == "__main__":
    main()