import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # ensure 3D support
import glob
import os
from matplotlib.animation import FuncAnimation

# Get sorted solution files based on the timestep number in the filename.
files = sorted(glob.glob("data/solution_*.dat"),
               key=lambda f: int(os.path.basename(f).split('_')[1].split('.')[0]))

if not files:
    print("No solution files found in the 'data/' directory.")
    exit(1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("z", fontsize=12)  # second spatial coordinate, renamed as z
ax.set_zlabel("Heat", fontsize=12)
title = ax.set_title("")

# Global variable for the current surface plot.
surf = None

def init():
    global surf
    # Read the first file to initialize the plot.
    with open(files[0]) as f:
        header_line = f.readline().strip()  # e.g., "# Time: 0.0000"
        time_val = float(header_line.split()[-1])
        data = np.genfromtxt(f)
    
    # Filter out any potential NaN rows.
    data = data[~np.isnan(data).any(axis=1)]
    num_points = data.shape[0]
    n = int(np.sqrt(num_points))
    if n * n != num_points:
        print("Data does not form a square grid.")
        exit(1)
    
    # Reshape the columns into n x n grids.
    X = data[:, 0].reshape(n, n)
    Z = data[:, 1].reshape(n, n)  # second column used as z
    U = data[:, 2].reshape(n, n)  # heat value
    
    # Create the surface plot.
    surf = ax.plot_surface(X, Z, U, cmap='viridis')
    title.set_text(f"Time t = {time_val:.4f}")
    
    # Fix axis limits (using the range from the first frame).
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Z), np.max(Z))
    ax.set_zlim(np.min(U), np.max(U))
    return surf, title

def update(frame):
    global surf
    # Remove the old surface before plotting the new one.
    if surf is not None:
        surf.remove()
        
    with open(files[frame]) as f:
        header_line = f.readline().strip()
        time_val = float(header_line.split()[-1])
        data = np.genfromtxt(f)
    
    data = data[~np.isnan(data).any(axis=1)]
    num_points = data.shape[0]
    n = int(np.sqrt(num_points))
    
    X = data[:, 0].reshape(n, n)
    Z = data[:, 1].reshape(n, n)
    U = data[:, 2].reshape(n, n)
    
    surf = ax.plot_surface(X, Z, U, cmap='viridis')
    title.set_text(f"Time t = {time_val:.4f}")
    return surf, title

# Create the animation with a 1000ms (1 second) interval between frames,
# and add repeat=False so the animation stops after the final frame.
ani = FuncAnimation(fig, update, frames=len(files), init_func=init,
                    interval=10, blit=False, repeat=False)

plt.tight_layout()
plt.show()
