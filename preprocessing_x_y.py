import os
from pathlib import Path
import nibabel as nib
import numpy as np

IN_DIR = Path("MiniVess/post_z/")  # input volumes (512x512xZ)
OUT_DIR = Path("MiniVess/post/")  # where cropped volumes go
OUT_DIR.mkdir(parents=True, exist_ok=True)

nii_exts = (".nii", ".nii.gz")

def crop_4_quadrants(data):
    """
    data: numpy array (X, Y, Z) with X=Y=512
    returns list of 4 arrays (256, 256, Z)
    """
    X, Y, Z = data.shape
    assert X == 512 and Y == 512, f"Expected 512x512, got {X}x{Y}"

    crops = []
    # top-left
    crops.append(data[0:256, 0:256, :])
    # top-right
    crops.append(data[0:256, 256:512, :])
    # bottom-left
    crops.append(data[256:512, 0:256, :])
    # bottom-right
    crops.append(data[256:512, 256:512, :])
    return crops

for fname in os.listdir(IN_DIR):
    if not fname.lower().endswith(nii_exts):
        continue

    in_path = IN_DIR / fname
    img = nib.load(str(in_path))
    data = img.get_fdata()
    affine = img.affine
    header = img.header

    if data.ndim == 4:
        data = data[..., 0]  # (X, Y, Z, T) -> (X, Y, Z)

    if data.shape[0] != 512 or data.shape[1] != 512:
        print(f"[SKIP] {in_path} has shape {data.shape}, not 512x512")
        continue

    crops = crop_4_quadrants(data)

    stem = fname
    if stem.endswith(".nii.gz"):
        stem = stem[:-7]
    elif stem.endswith(".nii"):
        stem = stem[:-4]

    for i, crop in enumerate(crops):
        out_name = f"{stem}_xy256_q{i}.nii.gz"
        out_path = OUT_DIR / out_name
        nib.save(nib.Nifti1Image(crop, affine, header), str(out_path))
        print(f"[SAVE] {out_path} with shape {crop.shape}")