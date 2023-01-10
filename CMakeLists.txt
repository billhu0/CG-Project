cmake_minimum_required(VERSION 3.16)
project(CS171-final-project LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)

add_subdirectory(libs)
add_subdirectory(src)

add_executable(${PROJECT_NAME}-main main.cpp)

target_link_libraries(${PROJECT_NAME}-main PRIVATE renderer)