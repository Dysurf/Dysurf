import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib import cm

# Set global font properties to Arial, size 14
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

# Load data
sqebin = np.loadtxt("SQEBIN_300K.dat")
sqebin = sqebin.reshape([101, 20, 50, 50], order='F')  # Reshape to (energy, ql, qk, qh)
print("sqebin shape:", sqebin.shape)  # Expected shape: (101, 20, 50, 50)

# Extract qh slice from indices 12:36 (qh from -2.5 to -1.5, assuming qh: -5 to 5)
sqebin_slice = sqebin[:, :, :, :]  # Shape: (101, 20, 50, 24)
print("sqebin_slice shape:", sqebin_slice.shape)

# Sum along ql (axis=1, indices 0:4)
sqebin_slice_sum = np.sum(np.sum(sqebin_slice[:, 0:4, :, :], axis=1)[:,:,0:2], axis=2)  # Shape: (101, 50, 24)
print("sqebin_slice_sum shape:", sqebin_slice_sum.shape)

# Define coordinates
# Define start and end points
start_point = np.array([0.5, -2.5])
end_point = np.array([-0.5, -1.5])

# Generate 50 equally spaced points
points = np.linspace(start_point, end_point, 50)  # ([X,Y])
x = points[:,0]
y = points[:,1]
z = np.linspace(0, 50, 101)  # energy

# Create meshgrid for surface plot, matching voxel_data shape
Z, X = np.meshgrid(z, x, indexing='ij')  # Shape: (101, 50)
Y = np.ones_like(Z) * y  # Shape: (101, 50), y fixed for each point

print("X shape:", X.shape)
print("Y shape:", Y.shape)
print("Z shape:", Z.shape)

# Log-transform intensity data
voxel_data = np.log10(sqebin_slice_sum + 1e-9)  # Shape: (101, 50, 24)
print("voxel_data shape:", voxel_data.shape)

# Create figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot surface for each qk slice
norm = plt.Normalize(vmin=-5, vmax=3)

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
ax.set_ylim(-3, -1)
ax.set_zlim(0, 50)

# Set axis tick intervals
ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.zaxis.set_major_locator(MultipleLocator(10))

# Set aspect ratio
ax.set_box_aspect([1, 1, 1])  # Adjust for better visualization

# Set view angle
ax.view_init(elev=12, azim=50, roll=0)

# Add colorbar with Arial font for tick labels
sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
cbar.ax.tick_params(labelsize=14, labelfontfamily='Arial')

plt.savefig('slice2.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('slice2.eps', facecolor='white')
plt.show()