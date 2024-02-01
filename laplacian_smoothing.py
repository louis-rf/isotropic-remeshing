import imageio
import matplotlib.pyplot as plt
import numpy as np

import pygeogram


def save_new_bunny(name):
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
    np.savez_compressed(name, **mesh)


def plot(mesh):
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


name = "data/bunny_remeshed_1000.npz"
# save_new_bunny(name)
mesh = dict(np.load(name).items())



from scipy.spatial import KDTree


def compute_vertex_normals(vertices, faces):
    # Initialize normals array
    normals = np.zeros_like(vertices)
    
    # Compute normals for each face
    for face in faces:
        v1, v2, v3 = vertices[face]
        normal = np.cross(v2 - v1, v3 - v1)
        normal /= np.linalg.norm(normal)
        for vertex in face:
            normals[vertex] += normal
    
    # Normalize the vertex normals
    for i in range(len(normals)):
        normals[i] /= np.linalg.norm(normals[i])
    
    return normals

def bilateral_smoothing(vertices, faces, sigma_s=1.0, sigma_r=0.1, iterations=10):
    # Compute vertex normals
    normals = compute_vertex_normals(vertices, faces)
    
    # Build a KD-Tree for efficient nearest neighbor queries
    kd_tree = KDTree(vertices)
    
    for _ in range(iterations):
        for i, vertex in enumerate(vertices):
            # Find the neighbors within the sigma_s spatial range
            indices = kd_tree.query_ball_point(vertex, sigma_s)
            W = 0.0
            new_position = np.zeros_like(vertex)
            
            for j in indices:
                if i == j:
                    continue
                neighbor_vertex = vertices[j]
                neighbor_normal = normals[j]
                
                # Calculate the spatial weight
                spatial_weight = np.exp(-np.linalg.norm(vertex - neighbor_vertex)**2 / (2 * sigma_s**2))
                
                # Calculate the range weight
                range_weight = np.exp(-np.arccos(np.clip(np.dot(normals[i], neighbor_normal), -1.0, 1.0))**2 / (2 * sigma_r**2))
                
                # Calculate the bilateral weight
                bilateral_weight = spatial_weight * range_weight
                
                # Update the new position
                new_position += bilateral_weight * neighbor_vertex
                W += bilateral_weight
            
            # Normalize the new position and update the vertex
            if W > 0:
                vertices[i] = new_position / W
    
    return vertices

# def compute_laplacian(vertices, faces):
#     """
#     Compute the Laplacian of each vertex in the mesh.
    
#     :param vertices: numpy array of shape [n_vertices, 3] representing the 3D coordinates of each vertex.
#     :param faces: numpy array of shape [n_faces, 3] representing the indices of the vertices for each triangular face.
#     :return: numpy array representing the Laplacian of each vertex.
#     """
#     n_vertices = vertices.shape[0]
#     laplacian = np.zeros_like(vertices)

#     # Create an adjacency list for the vertices
#     adjacency_list = [[] for _ in range(n_vertices)]
#     for face in faces:
#         for i in range(3):
#             adjacency_list[face[i]].extend([face[(i+1)%3], face[(i+2)%3]])

#     # Calculate the Laplacian for each vertex
#     for i in range(n_vertices):
#         adjacent_vertices = np.unique(adjacency_list[i])
#         # print(adjacent_vertices)
#         laplacian[i] = np.mean(vertices[adjacent_vertices], axis=0) - vertices[i]

#     return laplacian


# def laplacian_smoothing(vertices, faces, lambda_factor=0.1, iterations=10):
#     """
#     Apply Laplacian smoothing to the mesh.

#     :param vertices: numpy array of shape [n_vertices, 3] representing the 3D coordinates of each vertex.
#     :param faces: numpy array of shape [n_faces, 3] representing the indices of the vertices for each triangular face.
#     :param lambda_factor: Smoothing factor.
#     :param iterations: Number of iterations for smoothing.
#     :return: numpy array of the smoothed vertices.
#     """
#     for _ in range(iterations):
#         laplacian = compute_laplacian(vertices, faces)
#         vertices += lambda_factor * laplacian
#         # vertices *= 0.999
#     return vertices


def compute_laplacian(vertices, faces):
    """
    Compute the Laplacian of each vertex in the mesh.
    
    :param vertices: numpy array of shape [n_vertices, 3] representing the 3D coordinates of each vertex.
    :param faces: numpy array of shape [n_faces, 3] representing the indices of the vertices for each triangular face.
    :return: numpy array representing the Laplacian of each vertex.
    """
    n_vertices = vertices.shape[0]
    laplacian = np.zeros_like(vertices)

    # Create an adjacency list for the vertices
    adjacency_list = [[] for _ in range(n_vertices)]
    for face in faces:
        for i in range(3):
            adjacency_list[face[i]].extend([face[(i+1)%3], face[(i+2)%3]])

    # Calculate the Laplacian for each vertex
    for i in range(n_vertices):
        adjacent_vertices = np.unique(adjacency_list[i])
        laplacian[i] = np.mean(vertices[adjacent_vertices], axis=0) - vertices[i]

    return laplacian

def taubin_smoothing(vertices, faces, lambda_factor=0.1, mu_factor=-0.11, iterations=10):
    """
    Apply Taubin smoothing to the mesh.
    
    :param vertices: numpy array of shape [n_vertices, 3] representing the 3D coordinates of each vertex.
    :param faces: numpy array of shape [n_faces, 3] representing the indices of the vertices for each triangular face.
    :param lambda_factor: Smoothing factor for the Laplacian smoothing step.
    :param mu_factor: Smoothing factor for the shrinkage step, typically slightly larger than lambda and negative.
    :param iterations: Number of iterations for smoothing.
    :return: numpy array of the smoothed vertices.
    """
    for _ in range(iterations):
        # Smoothing step
        laplacian = compute_laplacian(vertices, faces)
        vertices += lambda_factor * laplacian
        
        # Shrinkage step
        laplacian = compute_laplacian(vertices, faces)
        vertices += mu_factor * laplacian

    return vertices


def inflate_to_sphere(vertices, iterations=10, inflation_factor=0.1):
    """
    Smooth and inflate a mesh towards a spherical shape.

    :param vertices: numpy array of shape [n_vertices, 3] representing the 3D coordinates of each vertex.
    :param iterations: Number of iterations for the inflation process.
    :param inflation_factor: Factor that controls how much vertices are moved towards the spherical shape in each iteration.
    :return: numpy array of the vertices inflated towards a spherical shape.
    """
    # Estimate the center and target radius of the mesh
    center = np.mean(vertices, axis=0)
    target_radius = np.mean(np.linalg.norm(vertices - center, axis=1))

    for _ in range(iterations):
        # Apply spherical projection and inflation
        for i, vertex in enumerate(vertices):
            # Calculate the direction from the center to the current vertex
            direction = vertex - center
            current_distance = np.linalg.norm(direction)
            direction /= current_distance  # Normalize the direction

            # Calculate the target position on the sphere surface
            target_position = center + direction * target_radius

            # Move vertex towards the target position
            vertices[i] += inflation_factor * (target_position - vertex)

    return vertices

def smoothing_and_remeshing(
    mesh, nb_points, lambda_factor=0.1, smoothing_iterations=10, total_iterations=10,
):
    vertices, faces = mesh["vertices"], mesh["facets"]
    width = vertices.max() - vertices.min()
    history = []
    for i in range(total_iterations):
        # vertices = laplacian_smoothing(
        #     vertices,
        #     faces,
        #     lambda_factor=lambda_factor,
        #     iterations=smoothing_iterations,
        # )
        vertices = inflate_to_sphere(vertices, iterations=3)
        # vertices = taubin_smoothing(vertices, faces, .5, -.55)
        # vertices = bilateral_smoothing(vertices, faces, sigma_s=width * 0.05, sigma_r=0.05)
        mesh = {"vertices": vertices, "edges": mesh["edges"], "facets": faces}
        mesh = pygeogram.remesh_smooth(
            mesh,  # mesh_data
            nb_points,  # nb_points
            1.0,   # tri_shape_adapt
            0.0,   # tri_size_adapt
            3,     # normal_iter
            5,     # Lloyd_iter
            30,    # Newton_iter
            7,     # Newton_m
            10000  # LFS_samples
        )
        vertices, faces = mesh["vertices"], mesh["facets"]
        history.append({"vertices": np.array(vertices), "facets": np.array(faces)})
        print(f"Done: {i + 1} / {total_iterations}")
    return {"vertices": vertices, "facets": faces}, history


vertices = mesh["vertices"]
faces = mesh["facets"]
# plot({"vertices": vertices, "facets": faces})

# vertices = laplacian_smoothing(vertices, faces, lambda_factor=0.5, iterations=100)
nb_points = 1000
mesh = pygeogram.remesh_smooth(
            mesh,  # mesh_data
            nb_points,  # nb_points
            1.0,   # tri_shape_adapt
            0.0,   # tri_size_adapt
            3,     # normal_iter
            5,     # Lloyd_iter
            30,    # Newton_iter
            7,     # Newton_m
            10000  # LFS_samples
        )
mesh, history = smoothing_and_remeshing(
    mesh, nb_points, lambda_factor=0.1, smoothing_iterations=10, total_iterations=10,
)

# # plot({"vertices": vertices, "facets": faces})
# for mesh in history:
#     plot(mesh)


def plot_meshes(mesh_list, output_filename):
    # Set up the GIF writer
    with imageio.get_writer(output_filename, mode='I', duration=0.03) as writer:
        for i, mesh in enumerate(mesh_list, start=1):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Extracting the x, y, z coordinates
            vertices = mesh["vertices"]
            faces = mesh["facets"]
            x = vertices[:, 0]
            y = vertices[:, 1]
            z = vertices[:, 2]

            # Plotting the surface
            ax.plot_trisurf(y, x, z, triangles=faces, cmap='viridis', edgecolor='none')

            # Adding labels  
            ax.set_axis_off()
            ax.set_title(f"Mesh: {i}")
            # Rotate and save each frame
            for angle in range(0, 360, 10):
                ax.view_init(30, angle)
                plt.draw()
                frame = plt_to_image(fig)
                writer.append_data(frame)
            plt.close(fig)
            print(f"Rendered mesh {i}")

def plt_to_image(fig):
    """ Convert a Matplotlib figure to an RGB image. """
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return image


plot_meshes(history, 'output_movie.gif')