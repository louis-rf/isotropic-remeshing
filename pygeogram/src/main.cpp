#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_remesh.h"
#include "geogram/mesh/mesh_repair.h"
#include "geogram/mesh/mesh_preprocessing.h"
#include "geogram/mesh/mesh_geometry.h"
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/points/co3ne.h>
#include <iostream>
#include <string> // Include for std::string
#include <vector>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

py::array_t<double> get_vertices(GEO::MeshVertices& vertices) {
    size_t num_vertices = vertices.nb();
    GEO::index_t dim = vertices.dimension();

    // Ensure the vertices are in double precision mode
    if (!vertices.double_precision()) {
        throw std::runtime_error("Vertices are not in double precision mode.");
    }

    auto result = py::array_t<double>({num_vertices, static_cast<size_t>(dim)});
    py::buffer_info buf_info = result.request();
    double *ptr = static_cast<double *>(buf_info.ptr);

    for (size_t i = 0; i < num_vertices; i++) {
        const double* vertex_ptr = vertices.point_ptr(i);
        for (GEO::index_t c = 0; c < dim; ++c) {
            ptr[i * dim + c] = vertex_ptr[c];
        }
    }

    return result;
}

py::array_t<unsigned int> get_edges(GEO::MeshEdges& edges) {
    size_t num_edges = edges.nb();

    // Creating an uninitialized NumPy array with the correct shape
    auto result = py::array_t<GEO::index_t>({num_edges, static_cast<size_t>(2)});
    py::buffer_info buf_info = result.request();
    GEO::index_t *ptr = static_cast<GEO::index_t *>(buf_info.ptr);

    for (size_t i = 0; i < num_edges; i++) {
        // Each edge has two vertices, accessed by the vertex method
        ptr[i * 2 + 0] = edges.vertex(i, 0); // Vertex index of the first endpoint
        ptr[i * 2 + 1] = edges.vertex(i, 1); // Vertex index of the second endpoint
    }
    return result;
}

py::array_t<unsigned int> get_facets(const GEO::MeshFacets& facets) {
    size_t num_facets = facets.nb();

    // Determine the maximum number of vertices per facet to allocate sufficient space
    size_t max_vertices_per_facet = 0;
    for (size_t i = 0; i < num_facets; ++i) {
        max_vertices_per_facet = std::max(max_vertices_per_facet, static_cast<size_t>(facets.nb_vertices(i)));
    }

    // Creating an uninitialized NumPy array with the correct shape
    auto result = py::array_t<GEO::index_t>({num_facets, max_vertices_per_facet});
    py::buffer_info buf_info = result.request();
    GEO::index_t* ptr = static_cast<GEO::index_t*>(buf_info.ptr);

    for (size_t i = 0; i < num_facets; ++i) {
        size_t num_vertices = facets.nb_vertices(i);
        for (size_t j = 0; j < num_vertices; ++j) {
            ptr[i * max_vertices_per_facet + j] = facets.vertex(i, j);
        }
        // Fill the remaining space with a sentinel value if the facet has fewer vertices
        for (size_t j = num_vertices; j < max_vertices_per_facet; ++j) {
            ptr[i * max_vertices_per_facet + j] = GEO::NO_VERTEX; // Define NO_VERTEX as per your context
        }
    }
    return result;
}

py::dict mesh2dict(GEO::Mesh& mesh) {
    py::array_t<double> vertices_array = get_vertices(mesh.vertices);
    py::array_t<GEO::index_t> edges_array = get_edges(mesh.edges);
    py::array_t<unsigned int> facets = get_facets(mesh.facets);

    // Create a dictionary to return the arrays
    py::dict result;
    result["vertices"] = vertices_array;
    result["edges"] = edges_array;
    result["facets"] = facets;

    return result;
}

void begin(GEO::Mesh& mesh) {
    mesh.vertices.set_double_precision();            
}

void _smooth_point_set(
    GEO::Mesh& mesh, GEO::index_t nb_iterations=2, GEO::index_t nb_neighbors=30
) {
    begin(mesh);
    if(nb_iterations != 0) {
        GEO::Co3Ne_smooth(mesh, nb_neighbors, nb_iterations);
    }
}

void _reconstruct_Co3Ne(
    GEO::Mesh& mesh,
    double radius=5.0,
    GEO::index_t nb_iterations=0,
    GEO::index_t nb_neighbors=30,
    GEO::MeshRepairMode mesh_repair_colocate=GEO::MESH_REPAIR_DEFAULT
) {
    begin(mesh);
    double R = GEO::bbox_diagonal(mesh);
    GEO::mesh_repair(mesh, mesh_repair_colocate, 1e-6*R);
    radius *= 0.01 * R;
    if(nb_iterations != 0) {
        GEO::Co3Ne_smooth(mesh, nb_neighbors, nb_iterations);
    }
    GEO::Co3Ne_reconstruct(mesh, radius);
    GEO::orient_normals(mesh);
}

void _remesh_smooth(
    GEO::Mesh& mesh,
    GEO::index_t nb_points = 30000,
    double tri_shape_adapt = 1.0,
    double tri_size_adapt = 0.0,
    GEO::index_t normal_iter = 3,
    GEO::index_t Lloyd_iter = 5,
    GEO::index_t Newton_iter = 30,
    GEO::index_t Newton_m = 7,
    GEO::index_t LFS_samples = 10000
) {
    if (mesh.facets.nb() == 0) {
        GEO::Logger::err("Remesh") << "mesh has no facet" << std::endl;
        return;
    }

    if (!mesh.facets.are_simplices()) {
        GEO::Logger::err("Remesh") << "mesh need to be simplicial, use repair" << std::endl;
        return;
    }

    begin(mesh);
    GEO::Mesh remesh;

    if (tri_shape_adapt != 0.0) {
        tri_shape_adapt *= 0.02;
        compute_normals(mesh);
        if (normal_iter != 0) {
            GEO::Logger::out("Nsmooth") << "Smoothing normals, " << normal_iter << " iteration(s)" << std::endl;
            simple_Laplacian_smooth(mesh, normal_iter, true);
        }
        set_anisotropy(mesh, tri_shape_adapt);
    } else {
        mesh.vertices.set_dimension(3);
    }

    if (tri_size_adapt != 0.0) {
        compute_sizing_field(mesh, tri_size_adapt, LFS_samples);
    } else {
        GEO::AttributesManager& attributes = mesh.vertices.attributes();
        if (attributes.is_defined("weight")) {
            attributes.delete_attribute_store("weight");
        }
    }

    GEO::remesh_smooth(mesh, remesh, nb_points, 0, Lloyd_iter, Newton_iter, Newton_m);

    GEO::MeshElementsFlags what = GEO::MeshElementsFlags(GEO::MESH_VERTICES | GEO::MESH_EDGES | GEO::MESH_FACETS | GEO::MESH_CELLS);
    mesh.clear();
    mesh.copy(remesh, true, what);

    GEO::orient_normals(mesh);
}

py::dict load(const std::string& file_path) {
    GEO::Mesh mesh;
    mesh.clear(false, false);
    GEO::MeshIOFlags flags;
    GEO::mesh_load(file_path, mesh, flags);
    std::cout << "Mesh loaded successfully!" << std::endl;

    return mesh2dict(mesh);
}

GEO::Mesh* py2mesh(const py::dict& mesh_data) {
    // GEO::Mesh mesh;
    auto* mesh = new GEO::Mesh();
    mesh->vertices.clear();
    mesh->edges.clear();
    mesh->facets.clear();  // Clear existing facets

    // Assuming mesh_data contains numpy arrays for vertices, edges, and facets
    py::array_t<double> vertices_array = mesh_data["vertices"].cast<py::array_t<double>>();
    py::array_t<GEO::index_t> edges_array = mesh_data["edges"].cast<py::array_t<GEO::index_t>>();
    py::array_t<GEO::index_t> facets_array = mesh_data["facets"].cast<py::array_t<GEO::index_t>>();

    // Populate vertices
    auto v_buf = vertices_array.unchecked<2>();
    for (ssize_t i = 0; i < v_buf.shape(0); ++i) {
        mesh->vertices.create_vertex(&v_buf(i, 0));
    }

    // Populate edges
    auto e_buf = edges_array.unchecked<2>();
    for (ssize_t i = 0; i < e_buf.shape(0); ++i) {
        mesh->edges.create_edge(e_buf(i, 0), e_buf(i, 1));
    }

    // Populate facets
    auto f_buf = facets_array.unchecked<2>();
    for (ssize_t i = 0; i < f_buf.shape(0); ++i) {
        ssize_t num_vertices = f_buf.shape(1);
        std::vector<GEO::index_t> facet_vertices(num_vertices);
        for (ssize_t j = 0; j < num_vertices; ++j) {
            facet_vertices[j] = f_buf(i, j);
        }

        // Create the facet based on the number of vertices
        if (num_vertices == 3) {
            mesh->facets.create_triangle(facet_vertices[0], facet_vertices[1], facet_vertices[2]);
        } else if (num_vertices == 4) {
            mesh->facets.create_quad(facet_vertices[0], facet_vertices[1], facet_vertices[2], facet_vertices[3]);
        } else {
            mesh->facets.create_polygon(facet_vertices.size(), facet_vertices.data());
        }
    }
    return mesh;
}

py::dict test(const py::dict& mesh_data) {
    auto* mesh = py2mesh(mesh_data);
    return mesh2dict(*mesh);
}

py::dict smooth_point_set(
    const py::dict& mesh_data,
    int nb_iterations = 2,
    int nb_neighbors = 30
) {
    auto* mesh = py2mesh(mesh_data);

    // Smooth point set with additional arguments
    _smooth_point_set(
        *mesh, 
        static_cast<GEO::index_t>(nb_iterations),
        static_cast<GEO::index_t>(nb_neighbors)
    );

    return mesh2dict(*mesh);
}

py::dict reconstruct_Co3Ne(
    const py::dict& mesh_data,
    float radius = 5.0,
    int nb_iterations = 0,
    int nb_neighbors = 30,
    int mesh_repair_colocate = static_cast<int>(GEO::MESH_REPAIR_DEFAULT)
) {
    auto* mesh = py2mesh(mesh_data);

    // Reconstruct using additional arguments
    _reconstruct_Co3Ne(
        *mesh, 
        static_cast<double>(radius),
        static_cast<GEO::index_t>(nb_iterations),
        static_cast<GEO::index_t>(nb_neighbors),
        static_cast<GEO::MeshRepairMode>(mesh_repair_colocate)
    );

    return mesh2dict(*mesh);
}

py::dict remesh_smooth(
    const py::dict& mesh_data,
    int nb_points = 30000,
    float tri_shape_adapt = 1.0,
    float tri_size_adapt = 0.0,
    int normal_iter = 3,
    int Lloyd_iter = 5,
    int Newton_iter = 30,
    int Newton_m = 7,
    int LFS_samples = 10000
) {
    auto* mesh = py2mesh(mesh_data);

    // Remesh with additional arguments
    _remesh_smooth(
        *mesh, 
        static_cast<GEO::index_t>(nb_points),
        static_cast<double>(tri_shape_adapt),
        static_cast<double>(tri_size_adapt),
        static_cast<GEO::index_t>(normal_iter),
        static_cast<GEO::index_t>(Lloyd_iter),
        static_cast<GEO::index_t>(Newton_iter),
        static_cast<GEO::index_t>(Newton_m),
        static_cast<GEO::index_t>(LFS_samples)
    );

    return mesh2dict(*mesh);
}

PYBIND11_MODULE(pygeogram, m) {
    // Initialize the GEOgram library
    GEO::initialize();
    // Setting up a default logger
    GEO::Logger::instance()->set_quiet(false);

    // Import necessary argument groups (needed for smooth_point_set and some others)
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("co3ne");

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: pygeogram

        .. autosummary::
           :toctree: _generate

           load
           smooth_point_set
           reconstruct_poisson
           remesh_smooth
    )pbdoc";

    m.def("load", &load, R"pbdoc(
        Load a mesh file.
    )pbdoc");

    m.def("smooth_point_set", &smooth_point_set, R"pbdoc(
        Smooth the points.
    )pbdoc");

    m.def("reconstruct_Co3Ne", &reconstruct_Co3Ne, R"pbdoc(
        Add the faces.
    )pbdoc");

    m.def("remesh_smooth", &remesh_smooth, R"pbdoc(
        Recompute the faces.
    )pbdoc");

    m.def("test", &test, R"pbdoc(
        Convert a mesh in python to a mesh in C++.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
