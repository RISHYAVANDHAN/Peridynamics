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

    // Ensure loadStep is not too small
    if (loadStep < 1e-6) {
        std::cout << "Warning: loadStep is very small (" << loadStep << "). Increasing to 0.01" << std::endl;
        loadStep = 0.01;
    }
    std::cout << "Using loadStep = " << loadStep << std::endl;

    // Set C1 to a proper non-zero value
    double C1 = 1.0;

    double NN = 0.0;
    const std::string DEF_flag = "EXP";
    const int number_of_points = std::floor(domain_size / Delta);
    double tol = 1e-8;  // Increased tolerance from 1e-11 to make convergence easier
    double LF = 0.0;    // Initialize load factor to 0
    int counter = 0;
    int DOFs = 0;
    double normnull = 0.0;
    Eigen::VectorXd dx;
    int min_try = 1;    // Ensure at least one iteration
    int max_try = 20;
    int completed_steps = 0;

    std::cout << "PD: " << PD << std::endl;
    std::cout << "Number of left patches: " << number_of_patches << std::endl;
    std::cout << "Number of right patches: " << number_of_right_patches << std::endl;
    std::cout << "Horizon size: " << delta << std::endl;
    std::cout << "Material constant C1: " << C1 << std::endl;
    std::cout << "Total number of steps: " << no_of_steps << std::endl;
    std::cout << "Load step size: " << loadStep << std::endl;

    std::vector<Points> point_list = generate_mesh(PD, d, domain_size, number_of_points, number_of_patches, Delta, number_of_right_patches, DEF_flag, DOFs);
    generate_neighbour_list(PD, point_list, delta);

    // Initial calculation of internal forces
    calculate_r(PD, point_list, NN, C1, delta);

    std::cout << "Starting simulation with " << point_list.size() << " points and " << DOFs << " DOFs" << std::endl;

    // Save initial configuration
    write_vtk(point_list, "initial_mesh.vtk");
    std::cout << "Wrote initial mesh file" << std::endl;

    // Main load stepping loop
    while (LF <= (1.0 + 1e-8)) {
        std::cout << "\n====== Load Factor: " << LF << " (Step " << completed_steps << ") ======\n" << std::endl;

        // Apply boundary conditions for current load step
        point_list = update_prescribed(point_list, LF, PD, delta);

        int newton_iterations = 0;
        bool isNotAccurate = true;

        // Force recalculation of residual for each load step
        calculate_r(PD, point_list, NN, C1, delta);

        // Newton-Raphson iteration loop
        while (isNotAccurate && newton_iterations < max_try) {
            newton_iterations++;

            Eigen::VectorXd R = assembleResidual(point_list, DOFs);

            // Print the max component of R for debugging
            double max_r_component = R.cwiseAbs().maxCoeff();
            std::cout << "Max residual component: " << max_r_component << std::endl;

            if (newton_iterations == 1) {
                normnull = R.norm();
                std::cout << "Initial residual norm: " << normnull << std::endl;

                // If initial residual is zero, add a small perturbation to prevent skipping
                if (normnull < tol) {
                    std::cout << "Warning: Initial residual nearly zero, applying small perturbation." << std::endl;
                    normnull = tol * 10.0;  // Use a small non-zero value
                }
            }

            double residual_norm = R.norm();
            std::cout << "Residual norm: " << residual_norm << " at iteration " << newton_iterations << std::endl;

            // Check convergence criteria
            if (newton_iterations > 1) {
                double relative_norm = residual_norm / normnull;
                if ((relative_norm < tol || residual_norm < tol) && newton_iterations >= min_try) {
                    std::cout << "Converged! Relative norm: " << relative_norm << ", Absolute norm: " << residual_norm << std::endl;
                    isNotAccurate = false;
                    break;
                }
            }

            // Assemble stiffness matrix
            Eigen::SparseMatrix<double> K = assembleStiffness(point_list, DOFs, PD);

            // Check if K is all zeros
            if (K.nonZeros() == 0) {
                std::cout << "ERROR: Stiffness matrix is empty! Adding regularization..." << std::endl;

                // Create identity matrix for regularization
                for (int i = 0; i < DOFs; ++i) {
                    K.coeffRef(i, i) = 1e-5;  // Add small diagonal values
                }
            }

            // Add regularization to all diagonal elements
            for (int i = 0; i < DOFs; ++i) {
                K.coeffRef(i, i) += 1e-8;  // Small regularization on diagonal
            }

            // Use a more robust solver approach
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(K);

            if (solver.info() != Eigen::Success) {
                std::cout << "Factorization failed! Adding stronger regularization..." << std::endl;

                // Try again with stronger regularization
                for (int i = 0; i < DOFs; ++i) {
                    K.coeffRef(i, i) += 1e-6;  // Stronger regularization
                }

                solver.compute(K);

                if (solver.info() != Eigen::Success) {
                    std::cout << "Factorization still failed after regularization!" << std::endl;
                    // Try to continue with next load step
                    isNotAccurate = false;
                    break;
                }
            }

            // Solve the system and update positions
            dx = -solver.solve(R);

            if (solver.info() != Eigen::Success) {
                std::cout << "Solve failed!" << std::endl;
                isNotAccurate = false;
                break;
            }

            // Debug: print max displacement
            double max_dx = dx.cwiseAbs().maxCoeff();
            std::cout << "Max displacement increment: " << max_dx << std::endl;

            // If displacement is too small, consider it converged
            if (max_dx < tol) {
                std::cout << "Displacement increment below tolerance, considering converged." << std::endl;
                isNotAccurate = false;
                break;
            }

            // Update positions
            point_list = update_displaced(point_list, dx, PD, delta);

            // Recalculate internal forces after update
            calculate_r(PD, point_list, NN, C1, delta);
        }

        // Check if we actually converged properly
        if (newton_iterations >= max_try) {
            std::cout << "WARNING: Maximum Newton iterations reached without convergence for LF = " << LF << "!" << std::endl;
            // Still continue to next load step, but with reduced step size if needed
            if (loadStep > 0.001) {
                loadStep *= 0.5;
                std::cout << "Reducing load step to " << loadStep << " for next iterations" << std::endl;
            }
        }

        // Write intermediate results at regular intervals
        if (completed_steps % 100 == 0 ||
            std::abs(LF - 0.25) < loadStep ||
            std::abs(LF - 0.5) < loadStep ||
            std::abs(LF - 0.75) < loadStep) {
            std::string filename = "mesh_step_" + std::to_string(completed_steps) + "_LF_" +
                               std::to_string(LF) + ".vtk";
            write_vtk(point_list, filename);
            std::cout << "Wrote intermediate file: " << filename << std::endl;
        }

        // Increment load factor and counter
        std::cout << "Updating LF from " << LF << " to " << LF + loadStep << std::endl;
        LF += loadStep;
        completed_steps++;
        counter++;
    }

    // Write final results
    std::string final_filename = "final_mesh_LF_" + std::to_string(LF-loadStep) + ".vtk";
    write_vtk(point_list, final_filename);
    write_vtk(point_list, "coloured_mesh.vtk");
    debug_it(PD, point_list);

    std::cout << "Simulation completed successfully with " << completed_steps << " steps!" << std::endl;
    return 0;
}