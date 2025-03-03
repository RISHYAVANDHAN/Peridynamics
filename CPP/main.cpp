//PERIDYNAMICS SIMULATION - Point-based Continuum Kinematics approach

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
#include <Eigen/SparseLU>
#include <string>
#include <fstream>
#include <chrono>
#include <numeric>
#include <memory>
#include <iomanip>

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "=== Starting Point-Based Peridynamics Simulation ===" << std::endl;

    // Simulation parameters
    int PD = 1;                  // Problem dimension
    double domain_size = 10.0;   // Size of physical domain
    double delta = 3.0;          // Horizon (interaction radius)
    double Delta = 1.0;          // Point spacing
    double d = 0.5;              // Total deformation magnitude

    // Material parameters
    double C1 = 1.0;             // Material constant

    // Boundary and patch parameters
    int number_of_patches = std::floor(delta / Delta);
    int number_of_right_patches = 1;
    const std::string DEF_flag = "EXP";  // Deformation mode (EXP, EXT, SHR)

    // Solver parameters
    const int max_steps = 100;           // Maximum number of load steps
    const int max_newton_iters = 20;     // Maximum Newton iterations per step
    const double tol = 1e-8;             // Convergence tolerance
    int DOFs = 0;                        // Will be calculated during mesh generation

    std::cout << "=== Configuration ===" << std::endl;
    std::cout << "Dimension: " << PD << "D" << std::endl;
    std::cout << "Domain size: " << domain_size << std::endl;
    std::cout << "Horizon (delta): " << delta << std::endl;
    std::cout << "Point spacing (Delta): " << Delta << std::endl;
    std::cout << "Deformation type: " << DEF_flag << std::endl;
    std::cout << "Material constant C1: " << C1 << std::endl;
    std::cout << "Load steps: " << max_steps << std::endl;
    std::cout << "Tolerance: " << tol << std::endl;

    // Generate mesh
    std::cout << "\n=== Generating mesh... ===" << std::endl;
    const int number_of_points = std::floor(domain_size / Delta);
    std::vector<Points> point_list = generate_mesh(PD, d, domain_size, number_of_points,
                                                  number_of_patches, Delta,
                                                  number_of_right_patches, DEF_flag, DOFs);

    // Generate neighbor lists
    std::cout << "Generating neighbor lists..." << std::endl;
    generate_neighbour_list(PD, point_list, delta);
    std::cout << "Mesh generated with " << point_list.size() << " points and " << DOFs << " DOFs" << std::endl;

    // Save initial configuration
    write_vtk(point_list, "initial_mesh.vtk");
    std::cout << "Initial mesh saved to 'initial_mesh.vtk'" << std::endl;

    // Calculate initial internal forces
    double NN = 0.0;  // Initial value for the nonlinear parameter
    calculate_r(PD, point_list, NN, C1, delta);

    // Main simulation loop - load stepping
    double load_step = 1.0 / max_steps;
    double load_factor = 0.0;
    int step_count = 0;

    std::cout << "\n=== Starting Simulation with Load Stepping ===" << std::endl;

    while (load_factor < 1.0 && step_count < max_steps) {
        // Increase load factor
        load_factor += load_step;
        if (load_factor > 1.0) load_factor = 1.0; // Ensure we don't exceed 1.0

        step_count++;
        std::cout << "\n--- Load Step " << step_count << " ---" << std::endl;
        std::cout << "Load Factor: " << std::fixed << std::setprecision(4) << load_factor << std::endl;

        // Apply boundary conditions for this load step
        point_list = update_prescribed(point_list, load_factor, PD, delta);

        // Newton-Raphson iterations
        int newton_iter = 0;
        bool converged = false;
        double residual_norm_initial = 0.0;

        // Reset internal forces for this load step
        calculate_r(PD, point_list, NN, C1, delta);

        while (!converged && newton_iter < max_newton_iters) {
            newton_iter++;

            // Assemble residual vector
            Eigen::VectorXd R = assembleResidual(point_list, DOFs);
            double residual_norm = R.norm();

            // Store initial residual for relative convergence check
            if (newton_iter == 1) {
                residual_norm_initial = residual_norm;
                std::cout << "Initial residual norm: " << residual_norm_initial << std::endl;

                // Handle case of very small initial residual
                if (residual_norm_initial < tol) {
                    std::cout << "Initial residual already below tolerance." << std::endl;
                    converged = true;
                    break;
                }
            }

            // Check convergence
            double relative_residual = residual_norm / residual_norm_initial;
            std::cout << "Iteration " << newton_iter << ": ||R|| = " << residual_norm
                      << " (relative: " << relative_residual << ")" << std::endl;

            if (residual_norm < tol || relative_residual < tol) {
                std::cout << "Converged in " << newton_iter << " iterations!" << std::endl;
                converged = true;
                break;
            }

            // Assemble stiffness matrix
            Eigen::SparseMatrix<double> K = assembleStiffness(point_list, DOFs, PD);

            if (K.nonZeros() == 0) {
                std::cerr << "ERROR: Stiffness matrix is empty!" << std::endl;
                break;
            }

            // Solve linear system without regularization
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(K);
            solver.factorize(K);

            if (solver.info() != Eigen::Success) {
                std::cerr << "ERROR: Factorization failed!" << std::endl;

                // Try to adapt by reducing load step instead of adding regularization
                load_step /= 2.0;
                std::cout << "Reducing load step to " << load_step << " and retrying" << std::endl;
                load_factor -= load_step; // Back up one step
                step_count--;
                break;
            }

            // Compute displacement increment
            Eigen::VectorXd dx = solver.solve(-R);

            if (solver.info() != Eigen::Success) {
                std::cerr << "ERROR: Solve failed!" << std::endl;
                break;
            }

            // Update positions
            point_list = update_displaced(point_list, dx, PD, delta);

            // Recalculate internal forces
            calculate_r(PD, point_list, NN, C1, delta);

            // Output energy for monitoring
            double energy = assembleEnergy(point_list);
            std::cout << "Total energy: " << energy << std::endl;
        }

        // Check if Newton-Raphson converged
        if (!converged) {
            std::cout << "WARNING: Newton-Raphson did not converge in " << max_newton_iters << " iterations" << std::endl;

            // Adapt step size for better convergence
            if (load_step > 0.001) {
                load_step /= 2.0;
                load_factor -= load_step; // Back up one step
                step_count--;
                std::cout << "Reducing load step to " << load_step << " and retrying" << std::endl;
                continue;
            }
        }

        // Save intermediate results at specific intervals
        if (step_count % 10 == 0 || std::abs(load_factor - 0.25) < load_step/2 ||
            std::abs(load_factor - 0.5) < load_step/2 || std::abs(load_factor - 0.75) < load_step/2 ||
            std::abs(load_factor - 1.0) < load_step/2) {
            std::string filename = "mesh_step_" + std::to_string(step_count) +
                                  "_LF_" + std::to_string(load_factor) + ".vtk";
            write_vtk(point_list, filename);
            std::cout << "Saved intermediate result: " << filename << std::endl;
        }
    }

    // Save final results
    write_vtk(point_list, "final_mesh.vtk");
    debug_it(PD, point_list);

    // Calculate and print execution time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\n=== Simulation completed in " << elapsed.count() << " seconds ===" << std::endl;
    std::cout << "Completed " << step_count << " load steps" << std::endl;

    return 0;
}