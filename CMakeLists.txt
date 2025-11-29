cmake_minimum_required(VERSION 3.10)
project(CHM_MKE VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Источники
file(GLOB SOURCES "src/*.cpp" "main.cpp")

# Исполняемый файл
add_executable(fem ${SOURCES})

# Заголовки
target_include_directories(fem PRIVATE ${CMAKE_SOURCE_DIR}/include)