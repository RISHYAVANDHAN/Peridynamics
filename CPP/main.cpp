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
#include <fstream>

int main() {
    // Run unit tests
    runUnitTests();

    auto start_time = std::chrono::high_resolution_clock::now();

    int PD = 1;
    double domain_size = 10.0, delta = 2.0, Delta = 1.0, d = 0.5;
    double C1 = 0.5, C2 = 0.5, C3 = 0.5;
    std::string DEF_flag = "EXT";

    std::cout << "\n=== Mesh Initialization ===\n";
    int DOFs = 0;
    auto points = generate_mesh(PD, d, domain_size, domain_size / Delta, delta / Delta, Delta, 1, DEF_flag, DOFs);

    generate_neighbour_list(PD, points, delta);
    std::cout << "Mesh contains " << points.size() << " points with " << DOFs << " DOFs\n";

    points = update_prescribed(points, 0.0, PD, delta);
    write_vtk(points, "initial_mesh.vtk");

    // Output file setup
    std::ofstream outfile("output.txt");
    outfile << "\n======================================================\n";
    outfile << "Number of DOFs: " << DOFs << "\n";
    outfile << "Number of Points: " << points.size() << "\n";
    outfile << "======================================================\n";

    // Number of iterations to run
    int num_iterations = 10;

    for (int iter = 0; iter < num_iterations; iter++) {
        std::cout << "\n=== Iteration: " << iter << " ===\n";
        outfile << "\n=== Iteration: " << iter << " ===\n";

        // Apply a fixed displacement increment to all free DOFs
        for (auto& point : points) {
            for (int dim = 0; dim < PD; ++dim) {
                if (point.BC[dim] == 1) {  // Only update free nodes
                    point.x[dim] += 0.5;  // Add fixed increment of 0.5
                }
            }
        }

        // Calculate residual and energy
        calculate_r(PD, points, C1, C2, C3, delta);

        // Assemble residual and stiffness matrix to check if they're working
        Eigen::VectorXd R = assembleResidual(points, DOFs, PD);
        Eigen::SparseMatrix<double> K = assembleStiffness(points, DOFs, PD);

        // Compute global energy
        double global_energy = assembleEnergy(points);

        // Print statistics to verify everything is working
        std::cout << "Residual Norm: " << R.norm() << " | Global Energy: " << global_energy << "\n";
        outfile << "Residual Norm: " << R.norm() << " | Global Energy: " << global_energy << "\n";

        std::cout << "Stiffness Matrix: " << K.rows() << "x" << K.cols()
                  << " with " << K.nonZeros() << " non-zeros\n";
        outfile << "Stiffness Matrix: " << K.rows() << "x" << K.cols()
                << " with " << K.nonZeros() << " non-zeros\n";

        // Force re-computation of neighbor interactions
        for (auto& point : points) {
            point.Ra_sum.setZero();
            point.Kab_sum.setZero();
            point.psi = 0.0;
        }

        // Output VTK file for visualization
        if (iter % 2 == 0) {
            write_vtk(points, "output_" + std::to_string(iter) + ".vtk");
        }
    }

    write_vtk(points, "final_solution.vtk");
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "\nCompleted in " << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << " seconds\n";
    outfile << "\nCompleted in " << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << " seconds\n";

    outfile.close();
    return 0;

}