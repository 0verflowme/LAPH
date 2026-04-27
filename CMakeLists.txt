cmake_minimum_required(VERSION 3.16)

project(laph
    VERSION 0.3.0
    DESCRIPTION "Lifted Affine-Phase Hypergraph exact quantum simulation kernel"
    LANGUAGES CXX
)

option(LAPH_BUILD_TESTS "Build native LAPH tests" ON)
option(LAPH_BUILD_EXAMPLES "Build native LAPH examples" ON)

add_library(laph_core
    src/solver.cpp
    src/phase.cpp
    src/optimizer.cpp
    src/gf2_dense.cpp
    src/density_kernel.cpp
    src/clifford_tableau.cpp
    src/laph.cpp
)
target_compile_features(laph_core PUBLIC cxx_std_17)
set_target_properties(laph_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(laph_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(LAPH_BUILD_EXAMPLES)
    add_executable(laph_demo examples/demo.cpp)
    target_link_libraries(laph_demo PRIVATE laph_core)
endif()

if(LAPH_BUILD_TESTS)
    enable_testing()
    add_executable(laph_tests tests/unit_tests.cpp)
    target_link_libraries(laph_tests PRIVATE laph_core)
    add_test(NAME laph_tests COMMAND laph_tests)
endif()

find_package(pybind11 CONFIG QUIET)

if(pybind11_FOUND)
    message(STATUS "Building Python bindings with pybind11")
    pybind11_add_module(_laph bindings/python_module.cpp)
    target_link_libraries(_laph PRIVATE laph_core)
    install(TARGETS _laph LIBRARY DESTINATION laph)
else()
    message(STATUS "pybind11 not found. Skipping Python bindings.")
endif()
