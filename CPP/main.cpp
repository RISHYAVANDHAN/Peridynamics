//PERIDYNAMICS SIMULATION

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "solver.h"
#include "Debug.h"
#include "VTKExport.h"
#include "update.h"
#include "Assemble.h"


#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <fstream>
#include <chrono>

int main(int argc, char *argv[])
{
    std::cout << "Starting Peridynamics Simulation" << std::endl;

    // Problem Definition Variables
    int PD = 2; // Problem Dimension (2D)
    //int Partition; // Number of partitions for mesh (not used in the provided code)
    double domain_size = 10.0;
    double delta = 3.0; // Horizon size
    double Delta = 1.0; // Grid spacing
    //int degree = 1; // Degree of approximation (linear or quadratic)
    int number_of_patches = std::floor(delta / Delta); // Number of left patches (based on horizon size)
    int number_of_right_patches = 1; // Right patches (user-defined)
    //int number_of_neighbours = 3; // Number of neighbors for peridynamics
    double d = 0.5; // Magnitude of deformation
    double no_of_steps = 1000.0; // Total number of simulation steps
    double loadStep = 1.0 / no_of_steps;
    double deformed_steps = d / no_of_steps; // Step size for deformation
    double C1 = 0.0; // Material constant (spring constant)
    double NN = 0.0; // Material power law
    const std::string& DEF_flag = "EXP"; // "EXT" - Extension; "EXT" - Expansion; "SHR" - Shear
    //int C_flag = 0; // Material behavior flag (0 - constant, 1 - linear, etc.)
    const int number_of_points = std::floor(domain_size / Delta); // Number of mesh points
    double tol = 1e-11; // Tolerance
    int LF = 0; // LoadFactor
    int counter = 0;
    int DOFs = 0;
    double normnull = 0.0;
    Eigen::VectorXd dx;
    int min_try = 0;
    int max_try = 20;

    // Debug information about parameters
    std::cout << "PD: " << PD << std::endl;
    std::cout << "Number of left patches: " << number_of_patches << std::endl;
    std::cout << "Number of right patches: " << number_of_right_patches << std::endl;
    std::cout << "Horizon size: " << delta << std::endl;

    // Generate mesh based on problem definition
    std::vector<Points> point_list = generate_mesh(PD, d, domain_size, number_of_points, number_of_patches, Delta, number_of_right_patches, DEF_flag, DOFs); // Assuming 'EXT' for deformation flag
    generate_neighbour_list(PD, point_list, delta);
    calculate_r(PD, point_list, NN, C1, delta);

    while (LF <= (1.0 + 1e-8))
    {
        std::cout << "LF: " << LF << std::endl;
        point_list = update_prescribed(point_list, LF, delta, PD);
        int error_counter = 1;
        bool isNotAccurate = true;

        while (isNotAccurate)
        {
            Eigen::VectorXd R = assembleResidual(point_list, DOFs);

            if (error_counter == 1)
            {
                normnull = R.norm();
                std::cout << "Initial residual norm: " << normnull << " at iteration " << counter << std::endl;
            }

            // Print the residual norm at each iteration
            std::cout << "Residual norm: " << R.norm() << " at iteration " << counter << std::endl;

            // Check stopping conditions
            if (error_counter > 1)
            {
                double relative_norm = R.norm() / normnull;
                double absolute_norm = R.norm();
                // Stopping conditions:
                // 1. Relative residual norm is below tolerance
                // 2. Absolute residual norm is below tolerance
                // 3. Maximum number of iterations is reached
                if ((relative_norm < tol || absolute_norm < tol || error_counter > max_try) && error_counter > min_try)
                {
                    isNotAccurate = false; // Exit the loop
                    break;
                }
            }

            Eigen::MatrixXd K = assembleStiffness(point_list, DOFs);
            Eigen::MatrixXd A = Eigen::MatrixXd(K);
            dx = -K.fullPivLu().solve(R);
            Eigen::MatrixXd KKtmp = K - K.transpose();

            point_list = update_displaced(point_list, dx, delta, PD);

            error_counter = error_counter + 1;
        }

        if (error_counter > max_try) {
            std::cout << "Convergence is not obtained!" << std::endl;
            break;
        }
        LF = LF + loadStep;
        counter = counter + 1;
    }
    // Write mesh to VTK file
    write_vtk(point_list, "C:/Users/srini/Downloads/FAU/Semwise Course/Programming Project/Peridynamics/coloured_mesh.vtk");
    // Debugging mesh
    debug_it(PD, point_list);

    return 0;
}

// End of file

