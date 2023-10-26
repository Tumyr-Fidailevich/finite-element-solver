cmake_minimum_required(VERSION 3.25)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
add_executable(finite-element-solver main.cpp ${INCLUDE_DIR}/lin-alg.h lin-alg.cpp
        ${INCLUDE_DIR}/io_parser.h io_parser.cpp ${INCLUDE_DIR}/solver.h solver.cpp
        ${INCLUDE_DIR}/physics.h physics.cpp)
target_include_directories(finite-element-solver PUBLIC ${INCLUDE_DIR})
target_precompile_headers(finite-element-solver PRIVATE ${INCLUDE_DIR}/pch.h)