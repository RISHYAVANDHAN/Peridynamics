cmake_minimum_required(VERSION 3.29)
project(1D_CPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 CONFIG REQUIRED)

# Define the executable with only implementation (.cpp) files
add_executable(1D_CPP
        main_1d.cpp
        mesh_1D.cpp
)

# Link the Eigen3 library
target_link_libraries(1D_CPP PRIVATE Eigen3::Eigen)
