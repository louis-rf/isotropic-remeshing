import mesh_utils

v = mesh_utils.load("data/bunny.ply")
vv = mesh_utils.smooth_point_set(v)
vvv = mesh_utils.reconstruct_Co3Ne(vv)
vvvv = mesh_utils.remesh_smooth(vvv)

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extracting the x, y, z coordinates
x = vvvv["vertices"][:, 0]
y = vvvv["vertices"][:, 1]
z = vvvv["vertices"][:, 2]

faces = vvvv["facets"]

# Plotting the surface
ax.plot_trisurf(x, y, z, triangles=faces, cmap='viridis', edgecolor='none')

# Adding labels and showing the plot
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()