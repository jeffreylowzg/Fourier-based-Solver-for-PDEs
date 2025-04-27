import numpy as np
import matplotlib.pyplot as plt
import glob, os
from matplotlib.animation import FuncAnimation

# User settings
threshold_percent = 1.0  # Error threshold in %

# Load files
pseudo_files = sorted(glob.glob("data/pseudo_solution_*.dat"),
                      key=lambda f: int(os.path.basename(f).split('_')[2].split('.')[0]))
full_files = sorted(glob.glob("data/full_solution_*.dat"),
                    key=lambda f: int(os.path.basename(f).split('_')[2].split('.')[0]))

if not pseudo_files or not full_files:
    raise RuntimeError("Missing pseudo or full solution files.")

# Precompute global color limits for consistent z-axis
global_min = np.inf
global_max = -np.inf
for pseudo_fn, full_fn in zip(pseudo_files, full_files):
    for fn in [pseudo_fn, full_fn]:
        data = np.genfromtxt(fn, skip_header=1)
        data = data[~np.isnan(data).any(axis=1)]
        U = data[:,2]
        global_min = min(global_min, U.min())
        global_max = max(global_max, U.max())

fig = plt.figure(figsize=(14,6))

ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

for ax in (ax1, ax2):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(global_min, global_max)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')

ax1.set_title("Pseudo-Spectral RK4")
ax2.set_title("Full-Spectral RK4")

surf1 = None
surf2 = None
error_scatter1 = None
error_scatter2 = None

def init():
    global surf1, surf2, error_scatter1, error_scatter2
    for ax in (ax1, ax2):
        ax.clear()
    return []

def update(i):
    global surf1, surf2, error_scatter1, error_scatter2

    ax1.clear()
    ax2.clear()

    # Load pseudo
    with open(pseudo_files[i]) as f:
        t_pseudo = float(f.readline().split()[-1])
        data_pseudo = np.genfromtxt(f)
    data_pseudo = data_pseudo[~np.isnan(data_pseudo).any(axis=1)]
    n = int(np.sqrt(len(data_pseudo)))
    Xp = data_pseudo[:,0].reshape(n,n)
    Yp = data_pseudo[:,1].reshape(n,n)
    Up = data_pseudo[:,2].reshape(n,n)

    # Load full
    with open(full_files[i]) as f:
        t_full = float(f.readline().split()[-1])
        data_full = np.genfromtxt(f)
    data_full = data_full[~np.isnan(data_full).any(axis=1)]
    Xf = data_full[:,0].reshape(n,n)
    Yf = data_full[:,1].reshape(n,n)
    Uf = data_full[:,2].reshape(n,n)

    # Compute relative error
    eps = 1e-12
    rel_error = np.abs(Uf - Up) / (np.abs(Up) + eps)
    mask = (rel_error * 100.0) > threshold_percent

    # Plot surfaces
    surf1 = ax1.plot_surface(Xp, Yp, Up, cmap='viridis')
    surf2 = ax2.plot_surface(Xf, Yf, Uf, cmap='viridis')

    # Highlight errors with red dots
    error_indices = np.argwhere(mask)
    if error_indices.size > 0:
        xi = error_indices[:,0]
        yi = error_indices[:,1]
        ax1.scatter(Xp[xi,yi], Yp[xi,yi], Up[xi,yi]+0.02, color='red', s=10, label='> threshold')
        ax2.scatter(Xf[xi,yi], Yf[xi,yi], Uf[xi,yi]+0.02, color='red', s=10)

    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    ax1.set_zlim(global_min, global_max)
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_zlim(global_min, global_max)

    ax1.set_title("Pseudo-Spectral RK4")
    ax2.set_title("Full-Spectral RK4")
    fig.suptitle(f"Time t = {t_pseudo:.4f}  |  Highlight: > {threshold_percent:.1f}% error", fontsize=16)

    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u')

    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u')

    return surf1, surf2

ani = FuncAnimation(fig, update, init_func=init,
                    frames=len(pseudo_files), interval=50, blit=False, repeat=False)

plt.tight_layout()
plt.show()
