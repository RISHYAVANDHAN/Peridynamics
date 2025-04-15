#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "solver.h"
#include "VTKExport.h"
#include "update.h"
#include "Assemble.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <string>
#include <iomanip>
#include <filesystem>

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "\n======================================================" << std::endl;
    std::cout << "Starting 1D Peridynamics Simulation" << std::endl;
    std::cout << "======================================================" << std::endl;

    // Problem parameters
    int PD = 1;                // 1D simulation
    double domain_size = 1.0;  // Domain size
    double delta = 0.301;        // Horizon size
    double Delta = 0.1;       // Grid spacing
    double d = 0.0001;            // Prescribed displacement

    // Calculate number of patches
    int number_of_patches = 3; //std::ceil(delta / Delta) but 3 for now
    int number_of_right_patches = 1;

    // Material properties
    double C1 = 0.5;          // Material constant
    double C2 = 0.0;          // Not used in 1D
    double C3 = 0.0;          // Not used in 1D
    std::string DEF_flag = "EXT";

    // Simulation settings
    int steps = 1;
    double load_step = d/steps;
    double min_load_step = 0.01;
    double tol = 1e-4;
    int min_iter = 3;
    int max_iter = 20;

    // Output settings
    std::string output_dir = "output_vtk/";
    bool save_iteration_vtk = true;
    int vtk_iteration_frequency = 5;

    // Create output directory
    std::filesystem::create_directories(output_dir);

    // Initialize mesh
    std::cout << "\n=== Mesh Initialization ===\n";
    int DOFs = 0;
    int DOCs = 0;
    auto points = generate_mesh(PD, domain_size, number_of_patches, Delta,
                              number_of_right_patches, DOFs, DOCs, d);
    std::cout << "Mesh contains " << points.size() << " points with " << DOFs << " DOFs\n";

    // Generate neighbors
    generate_neighbour_list(PD, points, delta);
    for (const auto& i : points) {
        std::cout << "Nr: " << i.Nr << ", X: [";
        for (const auto& val : i.X) {
            std::cout << val << " ";
        }
        std::cout << "], x: [";
        for (const auto& val : i.x) {
            std::cout << val << " ";
        }
        std::cout << "], Volume: " << i.volume << std::endl;
        std::cout << "BC: " << i.BC << " Flag: " << i.Flag << std::endl;

        std::cout << "Neighbours of " << i.Nr << " are: [";
        if (PD == 1) {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\nNumber of neighbours for point " << i.Nr << ": " << i.n1;
        }
        std::cout << "]\n";
        std::cout << std::endl;
    }

    // Export initial mesh
    std::string initial_mesh_filename = output_dir + "initial_mesh.vtk";
    write_vtk(points, initial_mesh_filename);

    // Print parameters
    std::cout << "Simulation Parameters:" << std::endl;
    std::cout << "Domain Size: " << domain_size << " | Delta: " << Delta
              << " | Horizon: " << delta << std::endl;
    std::cout << "Steps: " << steps << " | Load Step: " << load_step
              << " | Tolerance: " << tol << std::endl;
    std::cout << "Material constant C1: " << C1 << std::endl;
    std::cout << "======================================================" << std::endl;

    // Initialize load factor
    double LF = 0.0;
    double current_load_step = load_step;
    int counter = 0;
    bool simulation_success = true;

    // Main loading loop
    while (LF <= 1.0 + 1e-8 && counter < steps && simulation_success) {
        std::cout << "\nLoad Step " << counter + 1 << "/" << steps
                  << " | LF: " << LF << std::endl;

        auto last_converged_state = points;
        bool converged = false;
        int iter = 0;
        double residual_norm_initial = 0.0;

        // Apply boundary conditions
        points = update_prescribed(points, LF, PD, d);

        // Newton-Raphson iterations
        while (iter < max_iter && !converged) {
            iter++;

            // Reset force and stiffness terms
            for (auto& point : points) {
                point.Ra_sum.setZero();
                point.Ra_1.setZero();
                point.Kab_sum.setZero();
                point.Kab_1.setZero();
                point.psi = 0.0;
            }

            // Compute forces and stiffness
            calculate_r(PD, points, C1, C2, C3, delta);
            Eigen::VectorXd R = assembleResidual(points, DOFs, PD);
            double residual_norm = R.norm();

            if (iter == 1) {
                residual_norm_initial = residual_norm;
            }

            double normalized_residual = residual_norm_initial > 0 ?
                                      residual_norm / residual_norm_initial :
                                      residual_norm;

            std::cout << "  Iter " << iter << ": Residual = " << residual_norm
                     << " (norm: " << normalized_residual << ")" << std::endl;

            // Check convergence
            if (iter >= min_iter && normalized_residual < tol) {
                converged = true;
                std::cout << "  Converged in " << iter << " iterations." << std::endl;
                break;
            }

            // Assemble and solve linear system
            Eigen::SparseMatrix<double> K = assembleStiffness(points, DOFs, PD);
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(K);

            if (solver.info() != Eigen::Success) {
                std::cerr << "  Matrix factorization failed!" << std::endl;
                converged = false;
                break;
            }

            Eigen::VectorXd dx = solver.solve(-R);

            // Apply displacement increment
            points = update_displaced(points, dx, PD, delta);
            generate_neighbour_list(PD, points, delta);
        }

        // Handle convergence failure
        if (!converged) {
            current_load_step -= 0.05 ;
            std::cout << "  Reducing load step to " << current_load_step << std::endl;

            if (current_load_step < min_load_step) {
                std::cout << "  Minimum load step reached. Saving final state." << std::endl;
                simulation_success = false;
            } else {
                points = last_converged_state;
                LF -= current_load_step;
                continue;
            }
        }

        // Successful step
        if (save_iteration_vtk && (counter % vtk_iteration_frequency == 0)) {
            std::string filename = output_dir + "step_" + std::to_string(counter) + ".vtk";
            write_vtk(points, filename);
        }

        // Adjust load step adaptively
        if (iter < 5) {
            current_load_step = std::min(current_load_step * 1.5, load_step);
        }

        counter++;
        LF += current_load_step;

        // Ensure we hit exactly 1.0 at the end
        if (LF > 1.0 - current_load_step/2.0 && LF < 1.0) {
            LF = 1.0;
        }
    }

    // Final output
    write_vtk(points, output_dir + "final_solution.vtk");

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "\nSimulation completed in " << duration << " seconds" << std::endl;
    std::cout << "Final load factor reached: " << LF << std::endl;
    std::cout << "======================================================" << std::endl;

    return 0;
}