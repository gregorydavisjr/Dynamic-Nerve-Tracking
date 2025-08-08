> READ ME
# HOW TO USE TENSORFLOW PREDICTIONS
3 step process for using the US images to train a U-net model and make predictions
1. Process the masks to make ensure they are compatible
2. Train the U-net model using tensorflow in python
3. Run a prediction script using the .keras model in python
## Processing Masks
> mask_fix.py 
>
Before using the masks they must be "fixed" using the mask_fix.py script
- the blank masks must be inverted to work properly *(black instead of white; UN is nowhere instead of everywhere)*
- the masks must be converted to a high quality .png as well
>
Once opened the directory for the masks must be configured before processing them
- MASK_FOLDER = set the path that the masks are stored in 
>
run 'mask_fix.py' to fix the empty masks
-  in the terminal enter python mask_fix.py
>
*ensure you are in the proper directory*


## Running the training program
> ulnar_nerve_training.py
> 
  ensure you are in the proper directory using "cd"
>
  then in the terminal enter:
> 
python ulnar_nerve_training.py --data_root "C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_segmentation" --epochs 20 --batch 16 --img_size 128

- epoch = how many times it cycles through the images
- batch = how many images at a time
- img size = must match size of images and masks

*if not already in the proper directory you may need to enter the full path to the code, replacing "ulnar_nerve_training.py"*

python "C:\Users\ddavi\OneDrive\Desktop\School\MEMS FIRE\Robot Code\Unet\ulnar_nerve_training.py"

## Using the model to predict
> ulnar_nerve_predict.py
>
once the program is open, configure the save directory to the right path
- MODEL_PATH = set the path that the model was saved to
- IMAGES_FOLDER = set the path that the US images are in
- GT_MASK_FOLDER = set the path that the masks are in
- OUTPUT_PRED_DIR = set a path for the prediction images to save to
- IMG_SIZE = must match the size set during training
- THRESH = sigmoid threshold
>

### Going through predictions
>
Once the model has made its predictions the program will open a GUI
- On your keyboard press p to go to the previous figure and n to go to the next
- Press q to quit the program
>
Congrats!







