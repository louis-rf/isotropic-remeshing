cmake_minimum_required(VERSION 3.4...3.18)
project(pygeogram)

add_subdirectory(pybind11)
pybind11_add_module(pygeogram src/main.cpp)

# Use the GEOGRAM paths provided from setup.py
include_directories(${GEOGRAM_INCLUDE_DIR})

# Link GEOgram library
target_link_libraries(pygeogram PRIVATE ${GEOGRAM_LIBRARY})

# Pass version info from setup.py
target_compile_definitions(pygeogram
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
