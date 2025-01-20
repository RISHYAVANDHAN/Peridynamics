cmake_minimum_required(VERSION 3.29)
project(CPP)

set(CMAKE_CXX_STANDARD 26)

add_executable(CPP main.cpp
        Points.h
        mesh.h
        Neighbour.h
        solver.h
        compute.h
)
