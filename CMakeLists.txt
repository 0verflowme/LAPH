cmake_minimum_required(VERSION 3.16)

project(laph
    VERSION 0.3.0
    DESCRIPTION "Lifted Affine-Phase Hypergraph exact quantum simulation kernel"
    LANGUAGES CXX
)

add_library(laph INTERFACE)
target_compile_features(laph INTERFACE cxx_std_17)
target_include_directories(laph INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

add_executable(laph_demo examples/demo.cpp)
target_link_libraries(laph_demo PRIVATE laph)

enable_testing()
add_executable(laph_tests tests/unit_tests.cpp)
target_link_libraries(laph_tests PRIVATE laph)
add_test(NAME laph_tests COMMAND laph_tests)
