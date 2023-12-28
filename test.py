import matplotlib.pyplot as plt

import pygeogram

mesh = pygeogram.load("data/bunny.ply")
mesh = pygeogram.smooth_point_set(
    mesh,  # mesh_data
    2,     # nb_iterations
    30     # nb_neighbors
)
mesh = pygeogram.reconstruct_Co3Ne(
    mesh,  # mesh_data
    5.0,   # radius
    0,     # nb_iterations
    30,    # nb_neighbors
    1      # mesh_repair_colocate # MESH_REPAIR_DEFAULT
)
mesh = pygeogram.remesh_smooth(
    mesh,  # mesh_data
    1000,  # nb_points
    1.0,   # tri_shape_adapt
    0.0,   # tri_size_adapt
    3,     # normal_iter
    5,     # Lloyd_iter
    30,    # Newton_iter
    7,     # Newton_m
    10000  # LFS_samples
)
mesh = pygeogram.remesh_smooth(
    mesh,  # mesh_data
    100,  # nb_points
    1.0,   # tri_shape_adapt
    0.0,   # tri_size_adapt
    3,     # normal_iter
    5,     # Lloyd_iter
    30,    # Newton_iter
    7,     # Newton_m
    10000  # LFS_samples
)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extracting the x, y, z coordinates
vertices = mesh["vertices"]
faces = mesh["facets"]

x = vertices[:, 0]
y = vertices[:, 1]
z = vertices[:, 2]

# Plotting the surface
ax.plot_trisurf(x, y, z, triangles=faces, cmap='viridis', edgecolor='none')

# Adding labels and showing the plot
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()