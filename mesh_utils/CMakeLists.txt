cmake_minimum_required(VERSION 3.4...3.18)
project(mesh_utils)

add_subdirectory(pybind11)
pybind11_add_module(mesh_utils src/main.cpp)

# Use the GEOGRAM paths provided from setup.py
include_directories(${GEOGRAM_INCLUDE_DIR})

# Link GEOgram library
target_link_libraries(mesh_utils PRIVATE ${GEOGRAM_LIBRARY})

# Pass version info from setup.py
target_compile_definitions(mesh_utils
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
