import os
from pathlib import Path
import nibabel as nib
from collections import Counter

DATA_DIR = Path("MiniVess/post/")   # change this

nii_exts = (".nii", ".nii.gz")
xy_counts = Counter()
shapes = Counter()  # full shapes for quick overview

for root, dirs, files in os.walk(DATA_DIR):
    root_path = Path(root)
    for fname in files:
        if not fname.lower().endswith(nii_exts):
            continue

        fpath = root_path / fname
        try:
            img = nib.load(str(fpath))
            data = img.get_fdata()
        except Exception as e:
            print(f"[ERROR] {fpath}: {e}")
            continue

        if data.ndim == 4:
            # assume (X, Y, Z, T); just look at spatial
            data = data[..., 0]

        if data.ndim != 3:
            print(f"[WARN] {fpath} has unexpected shape {data.shape}")
            continue

        X, Y, Z = data.shape
        xy_counts[(X, Y)] += 1
        shapes[data.shape] += 1

        print(f"{fpath}")
        print(f"  shape: {data.shape} (X={X}, Y={Y}, Z={Z})")

print("\nSummary of (X, Y) sizes across dataset:")
for (X, Y), count in sorted(xy_counts.items()):
    print(f"  X={X:4d}, Y={Y:4d} : {count} volume(s)")

print("\nSummary of full shapes (X, Y, Z):")
for shape, count in sorted(shapes.items()):
    print(f"  {shape} : {count} volume(s)")