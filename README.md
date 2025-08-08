# Niryo Pose Logger

This application logs the pose and joint data of a Niryo Ned2 robot over time while it's in learning mode, and saves the data to an Excel file. It also generates a 3D plot of the recorded trajectory.

## 🧠 Key Features

- Logs robot TCP pose and joint angles at 20Hz
- Saves to Excel `.xlsx` file with timestamps
- Displays a 3D plot of the robot’s TCP trajectory
- GUI interface using `tkinter` for easy control
- Designed for Niryo Ned2 in learning mode

## 🛠️ Requirements

Install the following Python packages:
*Using windows powershell or equivalent*

```bash
pip install pandas matplotlib pyniryo2
```

Ensure the Niryo Ned2 robot is connected and reachable via the configured IP.
(whether using ethernet or hotspot)

## ⚙️ Configuration

In the script:

```python
ROBOT_IP = "x.x.x.x" #For ethernet: "169.254.200.200", for hotspot: "10.10.10.10"
SAVE_DIR = "C:\\Users\\ddavi\\OneDrive\\Desktop\\School\\MEMS FIRE\\Robot Code\\Data" #Change to the proper save directory for the .xlsx file
```

- Replace `ROBOT_IP` with your robot’s IP address. #ethernet or hotspot
- Update `SAVE_DIR` to your desired output directory.

## 🚀 How to Use
1. Connect to the robot:  
   Ensure the robot is on
   ETHERNET
   - plug in ethernet cord
   - ensure 'ROBOT_IP' is set to "169.254.200.200"
   HOTSPOT
   - connect to the robot's hotspot "Niryo Hotspot 23-e36-4c9" 
      - Password: niryorobot
   - ensure 'ROBOT_IP' is set to "10.10.10.10"

2. Run the script:  
   *ensure the correct directory is set using /cd followed by the directory* 
   ```bash
   python "Pose Logger.py" #Script can also be run by clicking the run button at the top right of vscode
   ```

3. Use the GUI:
   - Press **Start Recording** and move the robot by hand.
   - Press **Stop Recording** to save the data and generate a 3D plot.

4. Data will be saved in the specified `SAVE_DIR` as an Excel file.

## 📊 Output

- Excel file with:
  - Real Time
  - Timestamp
  - TCP Pose (X, Y, Z, Roll, Pitch, Yaw)
  - Joint Angles

- 3D matplotlib plot showing the TCP path.

## 🔒 Cleanup

The robot automatically exits learning mode and disconnects when you close the GUI.

## 🧑‍💻 Author

Script created by Gregory Davis for research in the MEMS FIRE program at the University of Pittsburgh.

---

### ⚠️ Disclaimer

This script requires a Niryo Ned2 robot. Always ensure the robot is in learning mode before moving it by hand.
