import os
import nibabel as nib
from collections import Counter

DATA_DIR = "MiniVess/raw/"  # change this to your dataset root

nii_exts = (".nii", ".nii.gz")
z_counts = Counter()

for root, dirs, files in os.walk(DATA_DIR):
    for fname in files:
        if fname.lower().endswith(nii_exts):
            fpath = os.path.join(root, fname)
            try:
                img = nib.load(fpath)
                data = img.get_fdata()
                shape = data.shape

                # For 3D images: (X, Y, Z)
                # For 4D images: (X, Y, Z, T) – we take Z = shape[2]
                if len(shape) < 3:
                    print(f"[WARN] {fpath} has shape {shape}, skipping")
                    continue

                z = shape[2]
                z_counts[z] += 1

                print(f"{fpath}")
                print(f"  shape: {shape}, z-slices: {z}")

            except Exception as e:
                print(f"[ERROR] Failed to load {fpath}: {e}")

print("\nSummary of z-slice counts across dataset:")
for z, count in sorted(z_counts.items()):
    print(f"  z = {z:3d} : {count} volume(s)")