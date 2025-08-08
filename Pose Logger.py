#log_pose w/ times
import os
import time
from datetime import datetime
import tkinter as tk
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyniryo2 import NiryoRobot  # type: ignore

# === CONFIGURATION ===
FREQUENCY = 20  # Hz
PERIOD = 1 / FREQUENCY
ROBOT_IP = "169.254.200.200"
SAVE_DIR = os.path.join("C:\\Users\\ddavi\\OneDrive\\Desktop\\School\\MEMS FIRE\\Robot Code\\Data") #CHANGE THIS TO SAVE DIRECTORY
os.makedirs(SAVE_DIR, exist_ok=True)

class PoseLoggerApp:
    def __init__(self, master):
        self.master = master
        master.title("Niryo Pose Logger")
        master.geometry("400x200")

        self.start_button = tk.Button(master, text="Start Recording", command=self.start_recording, font=("Arial", 16))
        self.start_button.pack(pady=10)

        self.stop_button = tk.Button(master, text="Stop Recording", command=self.stop_recording, font=("Arial", 16), state=tk.DISABLED)
        self.stop_button.pack(pady=10)

        self.status_label = tk.Label(master, text="Press 'Start Recording' to begin.", font=("Arial", 12))
        self.status_label.pack(pady=10)

        self.recording = False
        self.data = []
        self.start_time = None

        # Optional: Handle window close to clean up robot connection
        self.master.protocol("WM_DELETE_WINDOW", self.on_close)

    def start_recording(self):
        self.start_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)
        self.status_label.config(text="Recording poses & joints... Move robot. \nPress 'Stop Recording' to finish.")
        self.data = []
        self.start_time = time.time()
        self.recording = True
        self.log_pose_joints()

    def log_pose_joints(self):
        if not self.recording:
            return
        pose = robot.arm.get_pose()
        joints = robot.arm.get_joints()
        #pose2 = robot.arm.forward_kinematics(joints)
        if pose and joints is not None:
            current_time = time.time() - self.start_time
            self.data.append({
                "Real Time": datetime.now().strftime('%H:%M:%S.%f')[:-3],
                "Timestamp (s)": round(current_time, 4),
                "X(m)": round(pose.x, 5),
                "Y(m)": round(pose.y, 5),
                "Z(m)": round(pose.z, 5),
                "Roll(radians)": round(pose.roll, 5),
                "Pitch(radians)": round(pose.pitch, 5),
                "Yaw(radians)": round(pose.yaw, 5),
                "Base (radians)": round(joints[0], 5),
                "Shoulder (radians)": round(joints[1], 5),
                "Elbow pitch (radians)": round(joints[2], 5),
                "Elbow roll (radians)": round(joints[3], 5),
                "Wrist pitch (radians)": round(joints[4], 5),
                "Wrist roll (radians)": round(joints[5], 5),
                #"Pitch' (radians)": round(pose2.pitch, 5),
                #"Starting Pitch (radians)": round(starting_pose.pitch, 5)
            })
            print(f"Logged pose at {current_time:.3f}s: {pose}")
            print(f"Logged joints at {current_time:.3f}s: {joints}")
            ##print(f"Logged pitch at {current_time:.3f}s: {pose2.pitch}")
        self.master.after(int(PERIOD * 1000), self.log_pose_joints)

    def stop_recording(self):
        self.recording = False
        self.stop_button.config(state=tk.DISABLED)
        self.status_label.config(text="Saving data...")
        self.master.update()
        
       # === SAVE DATA TO EXCEL ===
        try:
            save_file = os.path.join(SAVE_DIR, f"pose_joint_data_{datetime.now().strftime('%Y.%m.%d_%H.%M.%S')}.xlsx")
            df = pd.DataFrame(self.data)
            print(f"Total poses recorded: {len(self.data)}")
            print(f"Total joints recorded: {len(self.data)}")
            if df.empty:
                print("❌ No pose or joint data collected. Nothing will be saved.")
                self.status_label.config(text="⚠️ No data collected. Try again.")
                return
            df.to_excel(save_file, index=False)
            self.status_label.config(text=f"✅ Data saved to:\n{save_file}\nPress 'Start Recording' for another session.")
        except Exception as e:
            self.status_label.config(text=f"❌ Error: {e}")
        self.start_button.config(state=tk.NORMAL)
        
        # === 3D PLOT TRAJECTORY ===
        try:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')

            xs = df['X(m)']
            ys = df['Y(m)']
            zs = df['Z(m)']

            ax.plot(xs, ys, zs, label='TCP Trajectory', marker='o')
            ax.set_title(f"3D TCP Trajectory: {datetime.now().strftime('%Y.%m.%d %H.%M.%S')}")
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_zlabel("Z (m)")
            ax.legend()
            plt.tight_layout()

            # Save plot as PNG
           # plot_path = os.path.join(SAVE_DIR, f"pose_joint_data_{datetime.now().strftime('%Y.%m.%d_%H.%M.%S')}.png")
           # plt.savefig(plot_path, dpi=300)
            plt.show()

           # print(f"✅ 3D trajectory plot saved to: {plot_path}")
        except Exception as e:
            print(f"❌ Failed to generate 3D plot: {e}")

    def on_close(self):
        # Clean up robot connection when closing the GUI
        try:
            robot.arm.set_learning_mode(False)
            robot.end()
        except Exception:
            pass
        self.master.destroy()

if __name__ == "__main__":
    try:
        # === CONNECT TO ROBOT ===
        print("Connecting to robot...")
        robot = NiryoRobot(ROBOT_IP)
        robot.arm.calibrate_auto()
        robot.arm.set_learning_mode(True)
        time.sleep(0.5) # pause to let mode settle
        print("⚠️ Make sure to move robot by hand — learning mode is ON")
        print("Connected to Ned^2! Opening GUI.")
        
        #starting_pose = robot.arm.get_pose()
        #print(f"Starting pitch: {starting_pose.pitch} (radians)")
        

        root = tk.Tk()
        app = PoseLoggerApp(root)
        root.mainloop()
        print("Robot connection closed. Exiting application.")
    except Exception as e:
        print("\n❌ An error occurred:")
        print(e)
        import traceback
        traceback.print_exc()
    input("\nPress Enter to exit...")
# === END OF LOG_POSE SCRIPT ===