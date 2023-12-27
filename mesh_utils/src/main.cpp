#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_remesh.h"
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
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

// GEO::Mesh mesh_loader(const std::string& file_path) {
//     // Initialize the GEOgram library
//     GEO::initialize();

//     // Setting up a default logger
//     GEO::Logger::instance()->set_quiet(false);

//     GEO::Mesh mesh;
//     mesh.clear(false, false);
//     GEO::MeshIOFlags flags;
//     GEO::mesh_load(file_path, mesh, flags);
//     std::cout << "Mesh loaded successfully!" << std::endl;

//     return mesh;
// }


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


py::dict load(const std::string& file_path) {
    // GEO::Mesh _mesh = mesh_loader(file_path);
    // Initialize the GEOgram library
    GEO::initialize();

    // Setting up a default logger
    GEO::Logger::instance()->set_quiet(false);

    GEO::Mesh mesh;
    mesh.clear(false, false);
    GEO::MeshIOFlags flags;
    GEO::mesh_load(file_path, mesh, flags);
    std::cout << "Mesh loaded successfully!" << std::endl;

    return mesh2dict(mesh);
}

py::dict py2mesh(const py::dict& mesh_data) {
    GEO::Mesh mesh;
    mesh.vertices.clear();
    mesh.edges.clear();
    mesh.facets.clear();  // Clear existing facets

    // Assuming mesh_data contains numpy arrays for vertices, edges, and facets
    py::array_t<double> vertices_array = mesh_data["vertices"].cast<py::array_t<double>>();
    py::array_t<GEO::index_t> edges_array = mesh_data["edges"].cast<py::array_t<GEO::index_t>>();
    py::array_t<GEO::index_t> facets_array = mesh_data["facets"].cast<py::array_t<GEO::index_t>>();

    // Populate vertices
    auto v_buf = vertices_array.unchecked<2>();
    for (ssize_t i = 0; i < v_buf.shape(0); ++i) {
        mesh.vertices.create_vertex(&v_buf(i, 0));
    }

    // Populate edges
    auto e_buf = edges_array.unchecked<2>();
    for (ssize_t i = 0; i < e_buf.shape(0); ++i) {
        mesh.edges.create_edge(e_buf(i, 0), e_buf(i, 1));
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
            mesh.facets.create_triangle(facet_vertices[0], facet_vertices[1], facet_vertices[2]);
        } else if (num_vertices == 4) {
            mesh.facets.create_quad(facet_vertices[0], facet_vertices[1], facet_vertices[2], facet_vertices[3]);
        } else {
            mesh.facets.create_polygon(facet_vertices.size(), facet_vertices.data());
        }
    }
    // Verify data by retrieving vertices, edges, and facets
    py::array_t<double> new_vertices = get_vertices(mesh.vertices);
    py::array_t<GEO::index_t> new_edges = get_edges(mesh.edges);
    py::array_t<GEO::index_t> new_facets = get_facets(mesh.facets);  // Assuming get_facets is implemented similarly

    // Optional: Print or otherwise validate new_vertices, new_edges, and new_facets

    // Create a dictionary to return the arrays
    py::dict result;
    result["vertices"] = new_vertices;
    result["edges"] = new_edges;
    result["facets"] = new_facets;

    return result;
}


PYBIND11_MODULE(mesh_utils, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: mesh_utils

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

    m.def("py2mesh", &py2mesh, R"pbdoc(
        Convert a mesh in python to a mesh in C++.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
