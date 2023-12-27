#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_remesh.h"
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <iostream>
#include <string> // Include for std::string


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

py::dict load(const std::string& file_path) {
    // Initialize the GEOgram library
    GEO::initialize();

    // Setting up a default logger
    GEO::Logger::instance()->set_quiet(false);

    GEO::Mesh mesh;
    mesh.clear(false, false);
    GEO::MeshIOFlags flags;
    GEO::mesh_load(file_path, mesh, flags);
    std::cout << "Mesh loaded successfully!" << std::endl;

    // Get vertices and edges
    py::array_t<double> vertices_array = get_vertices(mesh.vertices);
    py::array_t<GEO::index_t> edges_array = get_edges(mesh.edges);

    // Create a dictionary to return the arrays
    py::dict result;
    result["vertices"] = vertices_array;
    result["edges"] = edges_array;

    return result;
}

py::dict py2mesh(const py::dict& mesh_data) {
    GEO::Mesh mesh;
    mesh.vertices.clear();
    mesh.edges.clear();

    // Assuming mesh_data contains numpy arrays for vertices and edges
    py::array_t<double> vertices_array = mesh_data["vertices"].cast<py::array_t<double>>();
    py::array_t<GEO::index_t> edges_array = mesh_data["edges"].cast<py::array_t<GEO::index_t>>();

    // Populate vertices
    auto v_buf = vertices_array.unchecked<2>();
    for (ssize_t i = 0; i < v_buf.shape(0); ++i) {
        // Here, you need a method to add each vertex to mesh.vertices
        // Assuming create_vertex(const double* coords) can be used
        mesh.vertices.create_vertex(&v_buf(i, 0));
    }

    // Populate edges
    auto e_buf = edges_array.unchecked<2>();
    for (ssize_t i = 0; i < e_buf.shape(0); ++i) {
        // Here, you need a method to add each edge to mesh.edges
        // Assuming create_edge(index_t v1, index_t v2) can be used
        mesh.edges.create_edge(e_buf(i, 0), e_buf(i, 1));
    }

    // Verify data by retrieving vertices and edges
    py::array_t<double> new_vertices = get_vertices(mesh.vertices);
    py::array_t<GEO::index_t> new_edges = get_edges(mesh.edges);

    // Optional: Print or otherwise validate new_vertices and new_edges

    // Create a dictionary to return the arrays
    py::dict result;
    result["vertices"] = new_vertices;
    result["edges"] = new_edges;

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
