#include <iostream>
#include <cmath>
#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"


int main()
{
    std::cout << "Starting Peridynamics Simulation" << std::endl;


    int PD = 1; // Problem Definition
    int Partition = 10; // we already have number of points = number of partitions
    double domain_size = 10.0;
    double delta = 10.0; //horizon size
    double Delta = 10.0; // grid space
    int degree = 1; // 1 - linear, 2 - quadratic i´m still not sure if this is needed as it's not used in peridynamics, but in FEM its used, linear , quadratic
    int number_of_patches = std::round(static_cast<double>(delta) / Delta); // to the left
    int number_of_right_patches = 1; // this is user input 1/2/3
    int number_of_neighbours = 3;
    double d = 0.5; // magnitude of deformation
    double no_of_steps = 1000.0;
    //calculate step_size later = 1/no_of_steps
    double deformed_steps = d/no_of_steps;
    double C1{}; //material constant - resistance against change in length - the spring constant
    double NN{}; // material power law
    char flag = 'F'; // "F" or "D" for force or displacement
    int C_flag = 0; // for the CC function 0 - Constant, 1 - linear, 2 - Quadratic, 3 - Cubic, 4 - Quadratic, 5 - Exponential
    // the C_flag is used for analysing the stiffness between A-B, A-C (A, B, C are the points in the peridynamic framework)
    int number_of_points = 10;

    Points point;
    std::vector<Points> point_list = generate_mesh(PD, number_of_points, Partition, degree, domain_size, number_of_patches, Delta, number_of_right_patches);
    generate_neighbour_list(PD, number_of_points, point_list, number_of_patches,number_of_right_patches, delta);

    std::cout << std::endl;
    return 0;
}

