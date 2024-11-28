#include <iostream>
#include "Points.h"
#include "mesh.h"


int main()
{
    std::cout << "Starting Peridynamics Simulation" << std::endl;

    Points point;
    int PD = 1; // Problem Definition
    int Partition; // we already have number of points = number of partitions
    double domain_size;
    double horizon_size;
    double grid_space;
    int degree = 1; // 1 - linear, 2 - quadratic i´m still not sure if this is needed as it's not used in peridynamics, but in FEM its used, linear , quadratic
    int number_of_patches = static_cast<int>(horizon_size / grid_space); // to the left and the right
    double d = 0.5; // magnitude of deformation
    double no_of_steps;
    //calculate step_size later = 1/no_of_steps
    double deformed_steps = d/no_of_steps;
    double C1{}; //material constant - resistance against change in length - the spring constant
    double NN{}; // material power law
    char flag = 'F'; // "F" or "D" for force or displacement
    int C_flag = 0; // for the CC function 0 - Constant, 1 - linear, 2 - Quadratic, 3 - Cubic, 4 - Quadratic, 5 - Exponential
    // the C_flag is used for analysing the stiffness between A-B, A-C (A, B, C are the points in the peridynamic framework)

    generate_mesh(PD,point number_of_points, int Partition, int degree, double domain_size);


    std::cout << std::endl;
    return 0;
}
