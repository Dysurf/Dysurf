import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib import cm

# Set global font properties to Arial, size 14
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

# Assume sqebin is your data; replace with actual sqebin data
sqebin = np.loadtxt("SQEBIN_300K.dat")
sqebin = sqebin.reshape([101, 20, 50, 50], order='F')  # Reshape to (energy, ql, qk, qh)
print("sqebin shape:", sqebin.shape)  # Expected shape: (101, 20, 50, 50)

# Extract slice at qh=25 (corresponding to qh=-2 plane)
sqebin_slice = sqebin[:, :, 25, :]
print("sqebin_slice shape:", sqebin_slice.shape)  # Expected shape: (101, 20, 50)

# Sum along axis=1 (ql, indices 0:4)
sqebin_slice_sum = np.sum(sqebin_slice[:, 0:4, :], axis=1)
print("sqebin_slice_sum shape:", sqebin_slice_sum.shape)  # Expected shape: (101, 50)

# Define coordinate arrays
x = np.linspace(-0.5, 0.5, 50)  # qk
y = np.ones((50)) * -2  # Fixed at qh=-2
z = np.linspace(0, 50, 101)  # energy

# Create meshgrid for surface plot, matching voxel_data shape
Z, X = np.meshgrid(z, x, indexing='ij')  # Shape: (101, 50)
Y = np.ones_like(Z) * -2  # Shape: (101, 50), y fixed at -2

# Apply log transformation to intensity data to match previous color mapping
voxel_data = np.log10(sqebin_slice_sum + 1e-9)  # Shape: (101, 50)

# Create figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot surface with voxel_data for color mapping
norm = plt.Normalize(vmin=-5, vmax=3)  # Consistent with previous color range
surf = ax.plot_surface(X, Y, Z, facecolors=cm.jet(norm(voxel_data)), shade=False)

# Set axis labels with Arial font, size 14
ax.set_xlabel('H H 0', fontsize=14, fontfamily='Arial')
ax.set_ylabel('H -H 0', fontsize=14, fontfamily='Arial')
ax.set_zlabel('Energy (meV)', fontsize=14, fontfamily='Arial')

# Set tick labels font to Arial, size 14
ax.tick_params(axis='x', labelsize=14, which='both', labelfontfamily='Arial')
ax.tick_params(axis='y', labelsize=14, which='both', labelfontfamily='Arial')
ax.tick_params(axis='z', labelsize=14, which='both', labelfontfamily='Arial')

# Set axis ranges
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(-3, -1)  # Narrow y-axis range to focus on qh=-2 plane
ax.set_zlim(0, 50)

# Set axis tick intervals
ax.xaxis.set_major_locator(MultipleLocator(0.25))  # x-axis tick interval: 0.25
ax.yaxis.set_major_locator(MultipleLocator(0.5))   # y-axis tick interval: 0.5 (though fixed at -2)
ax.zaxis.set_major_locator(MultipleLocator(10))     # z-axis tick interval: 5

# Set perspective effect
ax.dist = 10  # Simulate perspective projection

# Set box aspect ratio
ax.set_box_aspect([1, 1, 1])  # Adjusted ratio for better plane visualization

# Set view angle
ax.view_init(elev=12, azim=50, roll=0)

# Add colorbar with Arial font for tick labels
sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
# cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
# cbar.ax.tick_params(labelsize=14, labelfontfamily='Arial')

plt.savefig('slice1.png',dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('slice1.eps', facecolor='white')
plt.show()