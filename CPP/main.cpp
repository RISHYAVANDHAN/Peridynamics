//PERIDYNAMICS SIMULATION

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "solver.h"
#include "Debug.h"
#include "VTKExport.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <fstream>


int main(int argc, char *argv[])
{
    std::cout << "Starting Peridynamics Simulation" << std::endl;

    // Problem Definition Variables
    int PD = 2; // Problem Dimension (2D)
    int Partition; // Number of partitions for mesh (not used in the provided code)
    double domain_size = 10.0;
    double delta = 3.0; // Horizon size
    double Delta = 1.0; // Grid spacing
    int degree = 1; // Degree of approximation (linear or quadratic)
    int number_of_patches = std::floor(delta / Delta); // Number of left patches (based on horizon size)
    int number_of_right_patches = 1; // Right patches (user-defined)
    int number_of_neighbours = 3; // Number of neighbors for peridynamics
    double d = 0.5; // Magnitude of deformation
    double no_of_steps = 1000.0; // Total number of simulation steps
    double deformed_steps = d / no_of_steps; // Step size for deformation
    double C1 = 0.0; // Material constant (spring constant)
    double NN = 0.0; // Material power law
    const std::string& DEF_flag = "EXP"; // "EXT" - Extension; "EXT" - Expansion; "SHR" - Shear
    int C_flag = 0; // Material behavior flag (0 - constant, 1 - linear, etc.)

    const int number_of_points = std::floor(domain_size / Delta); // Number of mesh points

    // Debug information about parameters
    std::cout << "PD: " << PD << std::endl;
    std::cout << "Number of left patches: " << number_of_patches << std::endl;
    std::cout << "Number of right patches: " << number_of_right_patches << std::endl;
    std::cout << "Horizon size: " << delta << std::endl;

    // Generate mesh based on problem definition
    std::vector<Points> point_list = generate_mesh(PD, d, domain_size, number_of_patches, Delta, number_of_right_patches, DEF_flag); // Assuming 'EXT' for deformation flag

    // Write mesh to VTK file
    write_vtk(point_list, "C:/Users/srini/Downloads/FAU/Semwise Course/Programming Project/Peridynamics/coloured_mesh.vtk");

    // Generate neighbor list and calculate R values (for peridynamics)
    generate_neighbour_list(PD, point_list, delta);
    calculate_r(PD, point_list, NN, C1, delta);

    // Debugging mesh
    debug_it(PD, point_list);

    // Additional functionality: Simulate the peridynamics steps
    for (int step = 0; step < no_of_steps; step++) {
        std::cout << "Simulation Step: " << step + 1 << "/" << no_of_steps << std::endl;

        // Here, you can simulate the peridynamics step based on your model
        // Update points, calculate forces/displacements, etc.
    }

    std::cout << "Simulation Completed!" << std::endl;
    return 0;
}


