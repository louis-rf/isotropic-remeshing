#include <pybind11/pybind11.h>
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_remesh.h"
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <iostream>
#include <string> // Include for std::string


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


// namespace py = pybind11;

// PYBIND11_MODULE(cmake_example, m) {
//     py::class_<GEO::Mesh>(m, "Mesh")
//         .def(py::init<GEO::index_t, bool>())
//         .def_readwrite("vertices", &GEO::Mesh::vertices)
//         // .def_readwrite("edges", &GEO::Mesh::edges)
//         // .def_readwrite("facets", &GEO::Mesh::facets)
//         // Add bindings for other mesh attributes as needed
//         ;
// }


int add(int i, int j) {
    return i + j;
}


// py::array_t<double> get_vertices_as_numpy(GEO::Mesh& mesh) {
//     auto vertices = mesh.vertices; // Assuming this gives you access to vertex data
//     size_t num_vertices = mesh.nb_v;/* Number of vertices */;
//     size_t dimension = mesh.dim;/* Dimension of each vertex, e.g., 3 for 3D vertices */;

//     // Create a NumPy array with the same data
//     return py::array_t<double>(
//         {num_vertices, dimension}, // shape of the array
//         {dimension * sizeof(double), sizeof(double)}, // strides
//         vertices.data(), // the data pointer
//         m // Python object that owns the data
//     );
// }


void load(const std::string& file_path) {
    // Initialize the GEOgram library
    GEO::initialize();

    // Setting up a default logger (if required)
    GEO::Logger::instance()->set_quiet(false);

    GEO::Mesh mesh;
    mesh.clear(false,false);
    GEO::MeshIOFlags flags;
    GEO::mesh_load(file_path, mesh, flags);
    std::cout << "Mesh loaded successfully!" << std::endl;

    // Your code to process the mesh...

    // // Finalize the GEOgram library (if there is a corresponding finalize function)
    // GEO::terminate();
    // return get_vertices_as_numpy(mesh)
}

namespace py = pybind11;

PYBIND11_MODULE(mesh_utils, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: mesh_utils

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("load", &load, R"pbdoc(
        Load a mesh file.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
