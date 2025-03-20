cmake_minimum_required(VERSION 3.29)
project(CPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 CONFIG REQUIRED)

# Define the executable before linking libraries
add_executable(CPP main.cpp
        Points.h
        mesh.h
        Neighbour.h
        solver.h
        Debug.h
        Assemble.h
        VTKExport.h
        update.h
        mesh.cpp
        UnitTest.h
        UnitTest.cpp
)

target_link_libraries(CPP PRIVATE Eigen3::Eigen)
