cmake_minimum_required(VERSION 3.10)
project(CHM_MKE VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Все исходники — жёстко прописаны
set(SOURCES
    src/mke_solver.cpp
    src/mesh_generator.cpp
    src/plot_utils.cpp
    src/globals.cpp
    main.cpp
)

add_executable(fem ${SOURCES})

target_include_directories(fem PRIVATE ${CMAKE_SOURCE_DIR}/include)