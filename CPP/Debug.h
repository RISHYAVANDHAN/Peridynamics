#ifndef DEBUG_H
#define DEBUG_H

#include <vector>
#include <iostream>
#include <iomanip> // For better formatting
#include <Eigen/Dense>
#include <Eigen/Sparse>

void debug_it(int PD, std::vector<Points>& point_list, const Eigen::VectorXd& R, const Eigen::SparseMatrix<double>& K) {
    std::cout << "\n=== Debug Output ===" << std::endl;

    // Print sizes of residual vector and stiffness matrix
    std::cout << "Size of Residual Vector (R): " << R.size() << " (Expected: " << R.size() << " DOFs)" << std::endl;
    std::cout << "Size of Stiffness Matrix (K): " << K.rows() << "x" << K.cols()
              << " (Expected: " << K.rows() << "x" << K.cols() << ")" << std::endl;
    std::cout << "Number of Non-Zero Entries in K: " << K.nonZeros() << std::endl;

    // Print residual vector (first 10 entries for brevity)
    std::cout << "\nFirst 10 Entries of Residual Vector (R):" << std::endl;
    for (int i = 0; i < std::min(10, static_cast<int>(R.size())); ++i) {
        std::cout << "R[" << i << "] = " << R(i) << std::endl;
    }

    // Print stiffness matrix (first 5x5 block for brevity)
    std::cout << "\nFirst 5x5 Block of Stiffness Matrix (K):" << std::endl;
    for (Eigen::Index i = 0; i < std::min(static_cast<Eigen::Index>(5), K.rows()); ++i) {
        for (Eigen::Index j = 0; j < std::min(static_cast<Eigen::Index>(5), K.cols()); ++j) {
            std::cout << std::setw(10) << K.coeff(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Print detailed information for each point
    for (const auto& i : point_list) {
        std::cout << "\nPoint Nr: " << i.Nr << std::endl;
        std::cout << "  Reference Position (X): [" << i.X.transpose() << "]" << std::endl;
        std::cout << "  Current Position (x): [" << i.x.transpose() << "]" << std::endl;
        std::cout << "  Volume: " << i.volume << std::endl;

        // Print boundary conditions and flag
        std::cout << "  BC: " << i.BC.transpose() << ", Flag: " << i.Flag << std::endl;

        // Print stiffness matrices
        std::cout << "  Stiffness Matrix (Kab_1):\n" << i.Kab_1 << std::endl;
        std::cout << "  Stiffness Matrix (Kab_sum):\n" << i.Kab_sum << std::endl;

        // Print residual vector
        std::cout << "  Residual Vector (Ra_sum): [" << i.Ra_sum.transpose() << "]" << std::endl;

        // Print neighbor lists
        std::cout << "  Neighbors of Point " << i.Nr << ":" << std::endl;
        if (PD == 1) {
            std::cout << "    1-Neighbors: ";
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    Number of 1-neighbors: " << i.n1 << std::endl;
        } else if (PD == 2) {
            std::cout << "    1-Neighbors: ";
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    2-Neighbors: ";
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    Number of 1-neighbors: " << i.n1 << ", Number of 2-neighbors: " << i.n2 << std::endl;
        } else if (PD == 3) {
            std::cout << "    1-Neighbors: ";
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    2-Neighbors: ";
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    3-Neighbors: ";
            for (const auto& n : i.neighbour_list_3N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n    Number of 1-neighbors: " << i.n1 << ", Number of 2-neighbors: " << i.n2
                      << ", Number of 3-neighbors: " << i.n3 << std::endl;
        }
    }

    std::cout << "\n=== End of Debug Output ===" << std::endl;
}


#endif // DEBUG_H