// 1D PERIDYNAMICS SIMULATION - Point-based Continuum Kinematics approach

#include "Points_1D.h"
#include "mesh_1D.h"
#include "Neighbour_1D.h"
#include "solver_1D.h"
#include "Debug_1D.h"
#include "update_1D.h"
#include "Assemble_1D.h"

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

int main() {
    std::ofstream outFile("output.txt");
    std::cout.rdbuf(outFile.rdbuf());

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "=== Starting 1D Peridynamics Simulation ===" << std::endl;

    // --- Simulation Parameters ---
    int PD = 1;
    double domain_size = 10.0;
    double delta = 3.0;
    double Delta = 1.0;
    double d = 0.5;  // Deformation factor
    double C1 = 1.0; // Material constant

    int number_of_patches = std::floor(delta / Delta);
    int number_of_right_patches = 1;
    const std::string DEFflag = "EXP";

    const int max_steps = 1000;
    const int max_newton_iters = 20;
    const double tol = 1e-8;
    int DOFs = 0;

    std::cout << "Domain: " << domain_size << ", Delta: " << Delta << ", Horizon: " << delta << std::endl;

    // --- Mesh Generation ---
    int number_of_points = std::floor(domain_size / Delta);
    std::vector<Points> point_list = generate_mesh(d, domain_size, number_of_points,
                                                  number_of_patches, Delta,
                                                  number_of_right_patches, DEFflag, DOFs);

    std::cout << "Generated mesh with " << point_list.size() << " points and " << DOFs << " DOFs.\n";

    // --- Neighbor List ---
    generate_neighbour_list_1D(point_list, delta);

    // --- Apply Initial Boundary Conditions ---
    point_list = update_prescribed(point_list, 0.01);

    // --- Initial Internal Forces ---
    double NN = 2.0;
    cal_pw_residual_tangent(point_list, NN, C1, delta);

    // --- Debug Check ---
    for (size_t i = 0; i < point_list.size(); ++i) {
        if (point_list[i].Kab_1 == 0.0 && !point_list[i].neighbour_list_1N.empty()) {
            std::cout << "Warning: Point " << i << " has neighbors but zero stiffness.\n";
        }
    }

    // --- Simulation Loop ---
    double load_step = 1.0 / max_steps;
    double load_factor = 0.1;
    int step_count = 0;

    while (load_factor < 1.0 && step_count < max_steps) {
        load_factor += load_step;
        if (load_factor > 1.0) load_factor = 1.0;

        step_count++;
        std::cout << "\n--- Load Step " << step_count << ", Load Factor: " << load_factor << " ---" << std::endl;

        point_list = update_prescribed(point_list, load_factor);
        cal_pw_residual_tangent(point_list, NN, C1, delta);

        int newton_iter = 0;
        bool converged = false;
        double residual_norm_initial = 0.0;

        while (!converged && newton_iter < max_newton_iters) {
            newton_iter++;
            Eigen::VectorXd R = assembleResidual(point_list, DOFs);
            double residual_norm = R.norm();

            Eigen::SparseMatrix<double> K = assembleStiffness(point_list, DOFs);
            debug_it(PD, point_list, R, K);

            if (newton_iter == 1) {
                residual_norm_initial = residual_norm;
                if (residual_norm_initial < tol) {
                    converged = true;
                    break;
                }
            }

            double relative_residual = residual_norm / residual_norm_initial;
            if (residual_norm < tol || relative_residual < tol) {
                converged = true;
                break;
            }
            std::cout << "Residual norm: " << residual_norm << std::endl;
            std::cout << "Relative residual: " << relative_residual << std::endl;
            std::cout << "Load factor: " << load_factor << std::endl;

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(K);
            solver.factorize(K);
            for (const auto& point : point_list) {
                if (std::abs(point.Kab_sum) < 1e-12) {
                    std::cout << "WARNING: Zero stiffness for point " << point.Nr << std::endl;
                }
            }


            if (solver.info() != Eigen::Success) {
                std::cerr << "Factorization failed. Trying smaller step." << std::endl;
                load_step /= 2.0;
                load_factor -= load_step;
                step_count--;
                break;
            }

            Eigen::VectorXd dx = solver.solve(-R);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Solve failed." << std::endl;
                break;
            }

            point_list = update_displaced(point_list, dx);
            cal_pw_residual_tangent(point_list, NN, C1, delta);
        }


        if (!converged) {

            std::cerr << "Did not converge. Backtracking." << std::endl;
            load_step /= 2.0;
            load_factor -= load_step;
            step_count--;
            continue;
        }
    }

    // --- Final Output ---
    std::cout << "\n=== Simulation complete. Steps: " << step_count << " ===" << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}
