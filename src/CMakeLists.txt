cmake_minimum_required(VERSION 3.25)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
add_executable(finite-element-solver main.cpp
        ${INCLUDE_DIR}/io_parser.h io_parser.cpp
        ${INCLUDE_DIR}/solver.h solver.cpp
        ${INCLUDE_DIR}/beam_element.h beam_element.cpp
        ${INCLUDE_DIR}/planar_element.h planar_element.cpp
        ${INCLUDE_DIR}/material.h material.cpp
        ${INCLUDE_DIR}/node.h node.cpp
        ${INCLUDE_DIR}/exceptions.h exceptions.cpp)

target_include_directories(finite-element-solver PUBLIC ${INCLUDE_DIR})
target_precompile_headers(finite-element-solver PRIVATE ${INCLUDE_DIR}/pch.h)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(finite-element-solver PRIVATE Eigen3::Eigen)