## Blank mask fix
import os
import cv2
import numpy as np
from pathlib import Path
from PIL import Image

MASK_FOLDER = r"C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_segmentation\masks"  # Update this to your mask folder path

for mask_file in Path(MASK_FOLDER).glob("*.tif"):
    mask = cv2.imread(str(mask_file), cv2.IMREAD_GRAYSCALE)
    
    if mask is None:
        print(f"Skipping unreadable: {mask_file}")
        continue

    # Check if all pixels are white (255)
    if np.all(mask == 255):
        print(f"Fixing white mask: {mask_file.name}")
        black = np.zeros_like(mask)
        cv2.imwrite(str(mask_file), black)


mask_dir = Path(MASK_FOLDER)
for tif in mask_dir.glob("mask_*.tif"):
    img = Image.open(tif).convert("L")   # force grayscale, 8‑bit
    png_path = tif.with_suffix(".png")
    img.save(png_path, compress_level=0) # lossless
    tif.unlink()                         # delete original
print("✅ Masks converted to PNG")

