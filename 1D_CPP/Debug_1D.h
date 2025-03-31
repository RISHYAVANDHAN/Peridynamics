#ifndef DEBUG_H
#define DEBUG_H

#include <vector>
#include "Points_1D.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>

void debug_it(int PD, std::vector<Points>& point_list, const Eigen::VectorXd& R, const Eigen::SparseMatrix<double>& K) {
    std::cout << "\n=== Debug Output (1D) ===" << std::endl;

    std::cout << "Residual Vector Size (R): " << R.size() << std::endl;
    std::cout << "Stiffness Matrix Size (K): " << K.rows() << " x " << K.cols() << std::endl;
    std::cout << "Non-zero Entries in K: " << K.nonZeros() << std::endl;

    std::cout << "\nFirst 10 Entries of R:" << std::endl;
    for (int i = 0; i < std::min(10, static_cast<int>(R.size())); ++i) {
        std::cout << "R[" << i << "] = " << R(i) << std::endl;
    }

    std::cout << "\nFirst 5x5 Block of K:" << std::endl;
    for (int i = 0; i < std::min(5, static_cast<int>(K.rows())); ++i) {
        for (int j = 0; j < std::min(5, static_cast<int>(K.cols())); ++j) {
            std::cout << std::setw(10) << K.coeff(i, j) << " ";
        }
        std::cout << std::endl;
    }


    for (const auto& i : point_list) {
        std::cout << "\nPoint Nr: " << i.Nr << std::endl;
        std::cout << "  Reference Position (X): " << i.X << std::endl;
        std::cout << "  Current Position (x): " << i.x << std::endl;
        std::cout << "  Volume: " << i.volume << std::endl;

        std::cout << "  Boundary Condition (BC): " << i.BC << ", Flag: " << i.Flag << std::endl;

        std::cout << "  Stiffness Kab_1: " << i.Kab_1 << std::endl;
        std::cout << "  Stiffness Kab_sum: " << i.Kab_sum << std::endl;

        std::cout << "  Residual Ra_sum: " << i.Ra_sum << std::endl;

        std::cout << "  1-Neighbors: ";
        for (const auto& n : i.neighbour_list_1N) {
            std::cout << "{ ";
            for (const int val : n) std::cout << val << " ";
            std::cout << "} ";
        }
        std::cout << "\n  Number of 1-neighbors: " << i.n1 << std::endl;
    }

    std::cout << "\n=== End of Debug Output (1D) ===" << std::endl;
}

#endif // DEBUG_H
