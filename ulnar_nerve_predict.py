import os
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from PIL import Image

# --------- CONFIG ---------
MODEL_PATH       = r"C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_segmentation\ulnar_unet.keras"  # Change to your saved model file (.keras or .h5)
IMAGES_FOLDER    = r"C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_segmentation\images"  # Change to test ultrasound images
GT_MASK_FOLDER   = r" "        # leave "" if none (ground truth masks)
OUTPUT_PRED_DIR  = r"C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_segmentation\predicitons" # Change to save predicted masks
IMG_SIZE         = 128  # Must match training image size
THRESH           = 0.5     # sigmoid threshold

# --------- Load model ---------
model = tf.keras.models.load_model(MODEL_PATH)
    # OR: model = tf.keras.models.load_model("PATH/TO/ulnar_unet.h5")

# ------------------- HELPERS ------------------- #
def load_image(path):
    img = Image.open(path).convert("RGB").resize((IMG_SIZE, IMG_SIZE))
    arr = np.array(img) / 255.0
    return arr

def load_mask_png(path):
    if not os.path.exists(path):    # ground‑truth may be missing
        return None
    m = Image.open(path).convert("L").resize((IMG_SIZE, IMG_SIZE))
    return (np.array(m) > 127).astype(np.uint8)  # binarize 0/1

def dice_coeff(pred, gt):
    inter = np.sum(pred * gt)
    return 2*inter / (np.sum(pred) + np.sum(gt) + 1e-7)

def iou(pred, gt):
    inter = np.sum(pred * gt)
    union = np.sum(pred) + np.sum(gt) - inter
    return inter / (union + 1e-7)

# ------------------- PREDICT LOOP -------------- #
image_files = sorted([f for f in os.listdir(IMAGES_FOLDER)
                      if f.lower().endswith((".png", ".jpg", ".jpeg", ".bmp"))])

results = []  # store metrics
overlays = [] # for interactive viewer
titles   = []

for fname in image_files:
    img_path = os.path.join(IMAGES_FOLDER, fname)
    img_arr  = load_image(img_path)
    inp      = np.expand_dims(img_arr, 0)

    pred     = model.predict(inp, verbose=0)[0, ..., 0]
    pred_bin = (pred > THRESH).astype(np.uint8)

    # ----- save predicted mask -----
    pred_img = (pred_bin * 255).astype(np.uint8)
    out_path = os.path.join(OUTPUT_PRED_DIR, f"pred_{os.path.splitext(fname)[0]}.png")
    Image.fromarray(pred_img).save(out_path)

    # ----- metrics if GT present -----
    gt_path  = os.path.join(GT_MASK_FOLDER, f"mask_{os.path.splitext(fname)[0]}.png") if GT_MASK_FOLDER else ""
    gt_bin   = load_mask_png(gt_path) if GT_MASK_FOLDER else None
    dice = iou_ = None
    if gt_bin is not None:
        dice = dice_coeff(pred_bin, gt_bin)
        iou_ = iou(pred_bin, gt_bin)
        results.append((fname, dice, iou_))
    else:
        results.append((fname, None, None))

    # ----- prepare overlay for viewer -----
    overlay = img_arr.copy()
    overlay[pred_bin == 1] = [1.0, 0.0, 0.0]  # red overlay
    overlays.append((img_arr, pred_bin, overlay))
    titles.append(fname)

print(f"✅  Saved {len(image_files)} predicted masks to {OUTPUT_PRED_DIR}")

# ----- Print metrics summary -----
valid = [r for r in results if r[1] is not None]
if valid:
    mean_dice = np.mean([r[1] for r in valid])
    mean_iou  = np.mean([r[2] for r in valid])
    print(f"Dice mean: {mean_dice:.4f} | IoU mean: {mean_iou:.4f}   "
          f"over {len(valid)} images with GT masks")
else:
    print("No ground‑truth masks found – metrics skipped.")

# ------------------- INTERACTIVE VIEWER -------- #
print("\nInteractive viewer: n=next, p=prev, q=quit")

idx = 0
fig, axes = plt.subplots(1, 3, figsize=(12,4))
plt.tight_layout()

def display(i):
    axes[0].imshow(overlays[i][0]); axes[0].set_title("Image");      axes[0].axis("off")
    axes[1].imshow(overlays[i][1], cmap="gray"); axes[1].set_title("Pred Mask"); axes[1].axis("off")
    axes[2].imshow(overlays[i][2]); axes[2].set_title("Overlay");    axes[2].axis("off")
    fig.suptitle(titles[i])
    fig.canvas.draw_idle()

display(idx)

def on_key(event):
    global idx
    if event.key == "n":
        idx = (idx + 1) % len(overlays)
        display(idx)
    elif event.key == "p":
        idx = (idx - 1) % len(overlays)
        display(idx)
    elif event.key == "q":
        plt.close(fig)
        sys.exit(0)

fig.canvas.mpl_connect("key_press_event", on_key)
plt.show()