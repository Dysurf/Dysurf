import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator

# Set global font properties to Arial, size 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

# Assume sqebin is your data; here we use random data for simulation
# Replace with your actual sqebin data
sqebin = np.loadtxt("SQEBIN_300K.dat")
sqebin = sqebin.reshape([101, 20, 50, 50], order='F')  # Reshape to (energy, ql, qk, qh)
print("sqebin shape:", sqebin.shape)  # Expected shape: (101, 20, 50, 50)

# Sum along axis=1 (dimension of size 20)
sqesum = np.sum(sqebin[:,0:5,:,:], axis=1)
print("sqesum shape:", sqesum.shape)  # Expected shape: (101, 50, 50)

# Define coordinate arrays
x = np.linspace(-0.5, 0.5, 50)
y = np.linspace(-3, -1, 50)
z = np.linspace(0, 50, 101)

# Downsample: take every 1 point (step size 1)
step = 1
x_sub = x[::step]  # 50 -> 25
y_sub = y[::step]  # 50 -> 25
z_sub = z[::step]  # 101 -> 51

# Create coordinate grid for voxel vertices (vertex count = voxel count + 1)
x_vertices = np.linspace(-0.5, 0.5, len(x_sub) + 1)
y_vertices = np.linspace(-3, -1, len(y_sub) + 1)
z_vertices = np.linspace(0, 50, len(z_sub) + 1)

# Create meshgrid for vertex coordinates
X, Y, Z = np.meshgrid(x_vertices, y_vertices, z_vertices, indexing='ij')
print("X shape:", X.shape)  # Expected shape: (26, 26, 52)

# Transpose and downsample sqesum to match (x, y, z) order
voxel_data = np.log10(np.transpose(sqesum, (2, 1, 0)) + 1e-9)  # Shape: (50, 50, 101)
voxel_data_sub = voxel_data[::step, ::step, ::step]  # Downsample, shape: (25, 25, 51)
print("voxel_data_sub shape:", voxel_data_sub.shape)  # Expected shape: (25, 25, 51)

# Normalize intensity for color mapping
norm = plt.Normalize(vmin=-5, vmax=3)
colors = plt.cm.jet(norm(voxel_data_sub))

# Create filled array for voxels
filled = np.ones(voxel_data_sub.shape, dtype=bool)  # Shape: (25, 25, 51), display all voxels
# Optional: Set a threshold to reduce voxel count, e.g., filled = np.abs(voxel_data_sub) > threshold

# Create figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot voxel plot
ax.voxels(X, Y, Z, filled, facecolors=colors, edgecolors=None)

# Set axis labels with Arial font, size 10
ax.set_xlabel('H H 0', fontsize=14, fontfamily='Arial')
ax.set_ylabel('H -H 0', fontsize=14, fontfamily='Arial')
ax.set_zlabel('Energy (meV)', fontsize=14, fontfamily='Arial')

# Set tick labels font to Arial, size 4
ax.tick_params(axis='x', labelsize=14, which='both', labelfontfamily='Arial')
ax.tick_params(axis='y', labelsize=14, which='both', labelfontfamily='Arial')
ax.tick_params(axis='z', labelsize=14, which='both', labelfontfamily='Arial')

# Set perspective effect
ax.dist = 14  # Simulate perspective projection

# Set box aspect ratio
ax.set_box_aspect([1, 1, 1])

# Add colorbar with Arial font for tick labels
sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
cbar.ax.tick_params(labelsize=14, labelfontfamily='Arial')

# Set axis ranges
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(-3, -1)
ax.set_zlim(0, 50)

# Set view angle
ax.view_init(elev=30, azim=48, roll=0)

# Set axis tick intervals
ax.xaxis.set_major_locator(MultipleLocator(0.25))  # x-axis tick interval: 0.25
ax.yaxis.set_major_locator(MultipleLocator(0.5))   # y-axis tick interval: 0.5
ax.zaxis.set_major_locator(MultipleLocator(10))    # z-axis tick interval: 12

plt.savefig('4dgan.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('4dgan.eps', facecolor='white')
plt.show()