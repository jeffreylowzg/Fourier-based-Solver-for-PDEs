import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob, os
from matplotlib.animation import FuncAnimation

# collect and sort your data files
files = sorted(glob.glob("data/solution_*.dat"),
               key=lambda f: int(os.path.basename(f).split('_')[1].split('.')[0]))
if not files:
    raise RuntimeError("No solution files found")

# 1) Pre–scan to find global heat extrema:
global_min = np.inf
global_max = -np.inf
for fn in files:
    data = np.genfromtxt(fn, skip_header=1)
    data = data[~np.isnan(data).any(axis=1)]
    U = data[:,2]
    global_min = min(global_min, U.min())
    global_max = max(global_max, U.max())

fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")    # call it “y” rather than “z”
ax.set_zlabel("Heat")
title = ax.set_title("")

surf = None

def init():
    global surf
    # load first frame
    with open(files[0]) as f:
        t = float(f.readline().split()[-1])
        data = np.genfromtxt(f)
    data = data[~np.isnan(data).any(axis=1)]
    n = int(np.sqrt(len(data)))
    X = data[:,0].reshape(n,n)
    Y = data[:,1].reshape(n,n)
    U = data[:,2].reshape(n,n)

    surf = ax.plot_surface(X, Y, U, cmap='viridis')
    title.set_text(f"Time t = {t:.4f}")

    # fix all three axes once and for all
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_zlim(global_min, global_max)

    return surf, title

def update(i):
    global surf
    if surf:
        surf.remove()
    with open(files[i]) as f:
        t = float(f.readline().split()[-1])
        data = np.genfromtxt(f)
    data = data[~np.isnan(data).any(axis=1)]
    n = int(np.sqrt(len(data)))
    X = data[:,0].reshape(n,n)
    Y = data[:,1].reshape(n,n)
    U = data[:,2].reshape(n,n)

    surf = ax.plot_surface(X, Y, U, cmap='viridis')
    title.set_text(f"Time t = {t:.4f}")
    # no need to reset limits here
    return surf, title

# Create the animation
ani = FuncAnimation(fig, update, init_func=init,
                    frames=len(files), interval=50, blit=False, repeat=False)

# Save the animation to a video file
ani.save('animation.mp4', writer='ffmpeg', fps=20)
