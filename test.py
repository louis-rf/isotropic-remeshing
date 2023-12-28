import pygeogram

v = pygeogram.load("data/bunny.ply")
vv = pygeogram.smooth_point_set(
    v,     # mesh_data
    2,     # nb_iterations
    30     # nb_neighbors
)
vvv = pygeogram.reconstruct_Co3Ne(
    vv,
    5.0,   # radius
    0,     # nb_iterations
    30,    # nb_neighbors
    1      # mesh_repair_colocate # MESH_REPAIR_DEFAULT
)
vvvv = pygeogram.remesh_smooth(
    vvv,   # mesh_data
    1000,  # nb_points
    1.0,   # tri_shape_adapt
    0.0,   # tri_size_adapt
    3,     # normal_iter
    5,     # Lloyd_iter
    30,    # Newton_iter
    7,     # Newton_m
    10000  # LFS_samples
)

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