#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "solver.h"
#include "Debug.h"
#include "VTKExport.h"
#include "update.h"
#include "Assemble.h"
#include "UnitTest.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "\n======================================================" << std::endl;
    std::cout << "Starting Peridynamics Simulation" << std::endl;
    std::cout << "======================================================" << std::endl;

    // Problem parameters
    int PD = 1;                // 1D simulation
    double domain_size = 5.0;  // Domain size
    double delta = 3.0;         //horizon size
    double Delta = 1.0;         // grid space
    double d = 0.25;            // Prescribed displacement
    int number_of_patches = floor(delta / Delta); // for now its 3, but not always check it once as well std::round(static_cast<double>(delta) / Delta); // to the left
    int number_of_right_patches = 1; //

    // Material properties
    double C1 = 1.0, C2 = 1.0, C3 = 1.0;
    std::string DEF_flag = "EXP";

    runUnitTests(PD);

    // Simulation settings
    int steps = 1000;
    double load_step = 0.005;
    double min_load_step = 0.001;
    double tol = 1e-8;        // Convergence tolerance
    int min_iter = 2;         // Minimum iterations for convergence
    int max_iter = 20;        // Maximum iterations per load step

    // Output settings
    std::string output_dir = "./output/";
    bool save_iteration_vtk = true;
    int vtk_iteration_frequency = 2;

    // Initialize mesh
    std::cout << "\n=== Mesh Initialization ===\n";
    int DOFs = 0;
    auto points = generate_mesh(PD, domain_size, number_of_patches, Delta, number_of_right_patches, d, DEF_flag, DOFs);

    // Store the previous converged state for each point
    for (auto& point : points) {
        point.x_prev = point.x;
    }

    // Generate neighbors
    generate_neighbour_list(PD, points, delta);
    std::cout << "Mesh contains " << points.size() << " points with " << DOFs << " DOFs\n";
    std::cout << "======================================================" << std::endl;

    // Print simulation parameters for reference
    std::cout << "Simulation Parameters:" << std::endl;
    std::cout << "PD: " << PD << " | Domain Size: " << domain_size << " | Delta: " << Delta << std::endl;
    std::cout << "Steps: " << steps << " | Load Step: " << load_step << " | Tolerance: " << tol << std::endl;
    std::cout << "C1: " << C1 << " | C2: " << C2 << " | C3: " << C3 << std::endl;
    std::cout << "======================================================" << std::endl;

    // Initialize load factor
    double LF = 0.1;
    double current_load_step = load_step;

    int counter = 0;
    bool simulation_completed = false;

    // Main loading loop (similar to MATLAB implementation)
    while (LF <= 1.0 + 1e-8 && !simulation_completed && counter < steps) {
        std::cout << "\n======================================================" << std::endl;
        std::cout << "Load factor: " << LF << " (Step " << counter << " of " << steps << ")" << std::endl;

        // Store last converged state
        auto last_converged_state = points;

        // Apply boundary conditions
        points = update_prescribed(points, LF, PD, delta);

        // === NEWTON-RAPHSON SOLVER IMPLEMENTATION === //
        int iter = 0;
        bool converged = false;
        double residual_norm_initial = 0.0;

        // Store the state to revert to if needed
        auto last_state = points;

        while (iter < max_iter && !converged) {
            iter++;
            std::cout << "Iteration " << iter << " of " << max_iter << std::endl;

            // Reset all force and stiffness terms
            for (auto& point : points) {
                point.Ra_sum.setZero();
                point.Ra_1.setZero();
                point.Ra_2.setZero();
                point.Ra_3.setZero();
                point.Kab_sum.setZero();
                point.Kab_1.setZero();
                point.Kab_2.setZero();
                point.Kab_3.setZero();
                point.psi = 0.0;
            }

            try {
                // Compute forces (residuals)
                calculate_r(PD, points, C1, C2, C3, delta);

                // Assemble global residual vector
                Eigen::VectorXd R = assembleResidual(points, DOFs, PD);
                double residual_norm = R.norm();

                // Track initial residual for convergence check
                if (iter == 1) {
                    residual_norm_initial = residual_norm;
                    std::cout << "Initial residual norm: " << residual_norm << std::endl;
                }

                double normalized_residual = residual_norm_initial > 0 ?
                                             residual_norm / residual_norm_initial :
                                             residual_norm;

                std::cout << "Residual norm: " << residual_norm
                          << " (normalized: " << normalized_residual << ")" << std::endl;

                // Check for convergence
                if ((iter >= min_iter) &&
                    ((normalized_residual < tol) || (residual_norm < tol))) {
                    std::cout << "Converged at iteration " << iter << std::endl;
                    converged = true;
                    break;
                }

                // Assemble stiffness matrix
                Eigen::SparseMatrix<double> K = assembleStiffness(points, DOFs, PD);

                // Solve linear system
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
                solver.analyzePattern(K);
                solver.factorize(K);

                if (solver.info() != Eigen::Success) {
                    std::cerr << "WARNING: Factorization failed! Trying QR solver instead." << std::endl;

                    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr_solver;
                    qr_solver.compute(K);

                    if (qr_solver.info() != Eigen::Success) {
                        throw std::runtime_error("Both LU and QR factorization failed!");
                    }

                    Eigen::VectorXd dx = qr_solver.solve(-R);

                    // Check solution quality
                    if ((K * dx + R).norm() / R.norm() > 0.1) {
                        std::cerr << "WARNING: Poor solution quality from QR solver!" << std::endl;
                    }

                    // Apply displacement increment
                    points = update_displaced(points, dx, PD, delta);
                } else {
                    // LU factorization succeeded
                    Eigen::VectorXd dx = solver.solve(-R);

                    // Apply displacement increment
                    points = update_displaced(points, dx, PD, delta);
                }

                // Regenerate neighbor lists after displacement
                generate_neighbour_list(PD, points, delta);
            }
            catch (const std::exception& e) {
                std::cerr << "ERROR during iteration " << iter << ": " << e.what() << std::endl;
                converged = false;
                break;
            }
        }

        // Check if we need to reduce load step size
        if (!converged) {
            std::cout << "Failed to converge in " << max_iter << " iterations!" << std::endl;
            current_load_step /= 2.0;

            if (current_load_step < min_load_step) {
                std::cout << "Load step too small. Aborting simulation." << std::endl;
                simulation_completed = true;
            } else {
                std::cout << "Reducing load step to " << current_load_step << " and retrying." << std::endl;
                points = last_converged_state;  // Restore last stable state
                LF -= current_load_step;        // Back up one step
            }
        }

        // Export results
        if (save_iteration_vtk && (counter % vtk_iteration_frequency == 0)) {
            std::string filename = output_dir + "output_step_" + std::to_string(counter) + ".vtk";
            write_vtk(points, filename);
            std::cout << "Wrote output to: " << filename << std::endl;
        }

        // Calculate energies (optional, similar to MATLAB)
        double total_energy = 0.0;
        for (const auto& point : points) {
            total_energy += point.psi;
        }
        std::cout << "Total strain energy: " << total_energy << std::endl;

        counter++;
        LF += current_load_step;

        // Ensure we hit exactly 1.0 for the final step
        if (LF > 1.0 - current_load_step / 2.0 && LF < 1.0) {
            LF = 1.0;
        }
    }

    // Final output
    write_vtk(points, output_dir + "final_solution.vtk");
    std::cout << "\n======================================================" << std::endl;
    std::cout << "Final solution exported to VTK." << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Simulation completed in " << duration << " seconds" << std::endl;
    std::cout << "======================================================" << std::endl;

    return 0;
}
