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
#include <numeric>
#include <memory>


int main(int argc, char *argv[]) {
    std::cout << "Starting Peridynamics Simulation" << std::endl;

    int PD = 1;
    double domain_size = 10.0;
    double delta = 3.0;
    double Delta = 1.0;
    int number_of_patches = std::floor(delta / Delta);
    int number_of_right_patches = 1;
    double d = 0.5;
    double no_of_steps = 1000.0;
    double loadStep = 1.0 / no_of_steps;
    double deformed_steps = d / no_of_steps;
    double C1 = 0.0;
    double NN = 0.0;
    const std::string DEF_flag = "EXP";
    const int number_of_points = std::floor(domain_size / Delta);
    double tol = 1e-11;
    int LF = 0;
    int counter = 0;
    int DOFs = 0;
    double normnull = 0.0;
    Eigen::VectorXd dx;
    int min_try = 0;
    int max_try = 20;

    std::cout << "PD: " << PD << std::endl;
    std::cout << "Number of left patches: " << number_of_patches << std::endl;
    std::cout << "Number of right patches: " << number_of_right_patches << std::endl;
    std::cout << "Horizon size: " << delta << std::endl;

    std::vector<Points> point_list = generate_mesh(PD, d, domain_size, number_of_points, number_of_patches, Delta, number_of_right_patches, DEF_flag, DOFs);
    generate_neighbour_list(PD, point_list, delta);

    // Set C1 to a non-zero value to ensure system is well-conditioned
    C1 = 1.0;  // Use an appropriate material constant for your simulation

    calculate_r(PD, point_list, NN, C1, delta);

    while (LF <= (1.0 + 1e-8)) {
        std::cout << "LF: " << LF << std::endl;
        point_list = update_prescribed(point_list, LF, PD, delta);
        int error_counter = 1;
        bool isNotAccurate = true;

        while (isNotAccurate) {
            Eigen::VectorXd R = assembleResidual(point_list, DOFs);

            if (error_counter == 1) {
                normnull = R.norm();
                std::cout << "Initial residual norm: " << normnull << " at iteration " << counter << std::endl;
            }

            std::cout << "Residual norm: " << R.norm() << " at iteration " << counter << std::endl;

            if (error_counter > 1) {
                double relative_norm = R.norm() / normnull;
                double absolute_norm = R.norm();
                if ((relative_norm < tol || absolute_norm < tol || error_counter > max_try) && error_counter > min_try) {
                    isNotAccurate = false;
                    break;
                }
            }

            // Skip solving if residual is already close to zero
            if (R.norm() < tol) {
                std::cout << "Residual already at tolerance, skipping solve step." << std::endl;
                isNotAccurate = false;
                break;
            }

            Eigen::SparseMatrix<double> K = assembleStiffness(point_list, DOFs, PD);

            // Check if K is empty or all zeros
            if (K.nonZeros() == 0) {
                std::cout << "Warning: Stiffness matrix is empty!" << std::endl;
                isNotAccurate = false;
                break;
            }

            // Use a more robust solver with regularization
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(K);

            if (solver.info() != Eigen::Success) {
                std::cout << "Warning: Factorization failed! Adding regularization..." << std::endl;
                // Add small regularization to diagonal
                for (int i = 0; i < K.rows(); ++i) {
                    K.coeffRef(i, i) += 1e-10;
                }
                solver.compute(K);

                if (solver.info() != Eigen::Success) {
                    std::cout << "Factorization still failed after regularization!" << std::endl;
                    isNotAccurate = false;
                    break;
                }
            }

            dx = -solver.solve(R);

            if (solver.info() != Eigen::Success) {
                std::cout << "Solve failed!" << std::endl;
                isNotAccurate = false;
                break;
            }

            point_list = update_displaced(point_list, dx, PD, delta);

            error_counter = error_counter + 1;
        }

        if (error_counter > max_try) {
            std::cout << "Convergence is not obtained!" << std::endl;
            break;
        }
        LF = LF + loadStep;
        counter = counter + 1;
    }

    write_vtk(point_list, "coloured_mesh.vtk");
    debug_it(PD, point_list);

    return 0;
}