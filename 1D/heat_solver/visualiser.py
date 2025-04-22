import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob

# Get sorted solution files
files = sorted(glob.glob("data/solution_*.dat"), key=lambda f: int(f.split('_')[1].split('.')[0]))

# Prepare data lists
time_steps = []
solution_matrix = []

# Read data from each solution file
for file in files:
    with open(file) as f:
        time_line = f.readline()
        time = float(time_line.strip().split()[-1])
        data = np.loadtxt(file, skiprows=1)
        x = data[:, 0]
        u = data[:, 1]
        time_steps.append(time)
        solution_matrix.append(u)

# Convert lists to arrays for plotting
X, T = np.meshgrid(x, time_steps)
U = np.array(solution_matrix)

# Plotting the 3D surface plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
surf = ax.plot_surface(X, T, U, cmap='viridis', edgecolor='none')

# Labels and title
ax.set_xlabel('Spatial coordinate $x$', fontsize=12)
ax.set_ylabel('Time $t$', fontsize=12)
ax.set_zlabel('Solution $u(x,t)$', fontsize=12)
ax.set_title('Heat Equation Solution over Space and Time', fontsize=14)

# Add color bar for clarity
fig.colorbar(surf, shrink=0.5, aspect=5, label='$u(x,t)$')

plt.tight_layout()
plt.show()