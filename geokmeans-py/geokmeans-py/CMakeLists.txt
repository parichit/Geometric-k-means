cmake_minimum_required(VERSION 3.15...3.30)
project(geokmeans LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# scikit-build-core provides the right Python; we only need the Development
# component for the extension module.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(CPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/geokmeans/_cpp")
if(NOT EXISTS "${CPP_DIR}/geokmeans.hpp")
    message(FATAL_ERROR
        "Vendored C++ sources not found in ${CPP_DIR}.\n"
        "Run `python sync_sources.py` before building.")
endif()

pybind11_add_module(_core src/geokmeans/_core.cpp)

target_include_directories(_core PRIVATE
    "${CPP_DIR}"
    "${CPP_DIR}/eigen")

# Optimise the compiled clustering routines.
if(MSVC)
    target_compile_options(_core PRIVATE /O2 /EHsc)
else()
    target_compile_options(_core PRIVATE -O3 -funroll-loops)
endif()

install(TARGETS _core DESTINATION geokmeans)
