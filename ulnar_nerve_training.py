#!/usr/bin/env python
"""
Lightweight U‑Net trainer for ulnar‑nerve segmentation.

Images: <DATA_ROOT>/images/1.png, 2.jpg, ...
Masks : <DATA_ROOT>/masks/mask_1.png, mask_2.png, ...
Blank (all‑black) masks are ignored automatically.
"""

import os, glob, re, argparse
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models
from tensorflow.keras.applications import MobileNetV2
from tensorflow.keras.callbacks import EarlyStopping
from PIL import Image

# --------------------------- CLI arguments --------------------------- #
p = argparse.ArgumentParser()
p.add_argument("--data_root", required=True)
p.add_argument("--epochs", type=int, default=20)
p.add_argument("--batch",  type=int, default=16)
p.add_argument("--img_size", type=int, default=128)
args = p.parse_args()

ROOT      = os.path.expanduser(args.data_root)
IMG_SIZE  = args.img_size
BATCH     = args.batch
EPOCHS    = args.epochs

# ------------------ grab file paths & pair --------------------------- #
img_dir  = os.path.join(ROOT, "images")
mask_dir = os.path.join(ROOT, "masks")

# Accept images with .png, .jpg, .jpeg, .bmp, .tif
image_paths = []
for ext in ("*.png", "*.jpg", "*.jpeg", "*.bmp", "*.tif"):
    image_paths.extend(glob.glob(os.path.join(img_dir, ext)))

# Accept masks with .png, *.tif
mask_paths = []
for ext in ("mask_*.png", "mask_*.tif"):
    mask_paths.extend(glob.glob(os.path.join(mask_dir, ext)))

def num(path: str):
    m = re.search(r'(\d+)(?=\.(?:png|jpg|jpeg|bmp)$)', os.path.basename(path), re.I)
    return int(m.group(1)) if m else None

mask_by_num = {num(m): m for m in mask_paths}

paired_images, paired_masks = [], []
skipped_blank = 0

print("Pairing & checking masks …")
for img in image_paths:
    n = num(img)
    mask_path = mask_by_num.get(n)
    if not mask_path:
        continue

    # ---- quick blank‑mask test (all pixels 0) ----
    try:
        arr = np.array(Image.open(mask_path).convert("L"))
        if arr.max() == 0:          # all‑black
            skipped_blank += 1
            continue
    except Exception as e:
        print(f"[skip] {mask_path}: {e}")
        continue

    paired_images.append(img)
    paired_masks.append(mask_path)

print(f"✅  kept {len(paired_images)} pairs   |   🌓 skipped {skipped_blank} blank masks")

if not paired_images:
    raise RuntimeError("No usable image–mask pairs found!")

# ------------------- tf.data pipeline ------------------------------- #
AUTOTUNE = tf.data.AUTOTUNE

def load_png(path, channels):
    data = tf.io.read_file(path)
    img  = tf.io.decode_png(data, channels=channels)
    img  = tf.image.resize(img, [IMG_SIZE, IMG_SIZE])
    img  = tf.cast(img, tf.float32) / 255.0
    return img

def load_pair(img_path, mask_path):
    return load_png(img_path, 3), load_png(mask_path, 1)

ds = tf.data.Dataset.from_tensor_slices((paired_images, paired_masks))\
        .map(load_pair, num_parallel_calls=AUTOTUNE)\
        .cache().shuffle(512).batch(BATCH).prefetch(AUTOTUNE)

train_ds = ds.take(int(0.8*len(paired_images)))
val_ds   = ds.skip(int(0.8*len(paired_images)))

# ------------------- U‑Net with MobileNetV2 encoder ------------------ #
def up(filters):  # small helper
    return models.Sequential([
        layers.UpSampling2D(), layers.Conv2D(filters, 3, padding="same", activation="relu"), layers.BatchNormalization()
    ])

def build_unet():
    base = MobileNetV2(input_shape=(IMG_SIZE, IMG_SIZE, 3), include_top=False, weights="imagenet")
    base.trainable = False
    skips = [base.get_layer(n).output for n in ("block_1_expand_relu","block_3_expand_relu","block_6_expand_relu","block_13_expand_relu")]
    x = base.output
    for f, s in zip([512,256,128,64], reversed(skips)):
        x = up(f)(x); x = layers.Concatenate()([x, s])
    x = layers.UpSampling2D()(x)
    out = layers.Conv2D(1,1,activation="sigmoid")(x)
    return models.Model(base.input, out)

model = build_unet()
model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])
model.summary()

# ------------------------ train ------------------------ #
early = EarlyStopping(monitor="val_loss", patience=5, restore_best_weights=True)
model.fit(train_ds, validation_data=val_ds, epochs=EPOCHS, callbacks=[early])

# ------------------------ save ------------------------- #
# choose an output folder 
out_dir = ROOT  # or any folder you like

# 1) Native Keras format  (.keras)
keras_path = os.path.join(out_dir, "ulnar_unet.keras")
model.save(keras_path)
print(f"✅  Saved native Keras file → {keras_path}")

# 2) Legacy HDF5 format  (.h5)
h5_path = os.path.join(out_dir, "ulnar_unet.h5")
model.save(h5_path)
print(f"✅  Saved HDF5 file        → {h5_path}")

# 3) TensorFlow SavedModel directory (for TF Serving / TFLite)
savedmodel_dir = os.path.join(out_dir, "ulnar_unet_savedmodel")
model.export(savedmodel_dir)          # Keras 3 API
print(f"✅  Exported SavedModel    → {savedmodel_dir}")

print(f"🏁  done!  model saved")
