#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"

// Function to compute the total energy of the system
inline double assembleEnergy(const std::vector<Points>& point_list) {
    double Psi = 0.0;

    for (const auto& point : point_list) {
        Psi += point.volume * point.psi; // Sum up energy contributions from each point
    }

    return Psi;
}

// Function to assemble the residual vector
inline Eigen::VectorXd assembleResidual(const std::vector<Points>& point_list, const int DOFs) {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs); // Initialize residual vector

    for (size_t i = 0; i < point_list.size(); i++) {
        // Get residual for this point
        Eigen::VectorXd R_P = point_list[i].Ra_sum;

        // Get dimension from BC vector
        int PD = point_list[i].BC.size();

        for (int ii = 0; ii < PD; ii++) {
            if (point_list[i].BC(ii) == 1) { // Check if this DOF is active
                int dof_idx = static_cast<int>(point_list[i].DOF(ii)); // Get DOF index
                if (dof_idx > 0 && dof_idx <= DOFs) { // Ensure DOF index is valid
                    R(dof_idx - 1) += R_P(ii); // Add contribution to residual vector
                }
            }
        }
    }

    return R;
}

// Function to assemble the stiffness matrix
inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& point_list, int DOFs, int PD) {
    std::vector<Eigen::Triplet<double>> triplets; // For sparse matrix storage

    // Iterate over each point in the point list
    for (size_t p = 0; p < point_list.size(); ++p) {
        const Eigen::Matrix3d& K_P = point_list[p].Kab_sum; // Total stiffness matrix for point p
        const Eigen::Vector3d& BCflg_p = point_list[p].BC; // Boundary condition flags for point p
        const Eigen::Vector3d& DOF_p = point_list[p].DOF;  // DOF for point p

        // Loop over each dimension (x, y, z) depending on PD (1D, 2D, or 3D)
        for (int pp = 0; pp < PD; ++pp) {
            if (BCflg_p(pp) == 1) { // Check if this DOF is active
                int row = static_cast<int>(DOF_p(pp)) - 1; // Convert to 0-based indexing

                if (row >= 0 && row < DOFs) { // Ensure row index is valid
                    // Self contribution (diagonal entry)
                    if (pp < K_P.rows() && pp < K_P.cols()) {
                        triplets.emplace_back(row, row, K_P(pp, pp));
                    }

                    // Check neighbors
                    for (size_t j = 0; j < point_list[p].neighbour_list.size(); ++j) {
                        int q = point_list[p].neighbour_list[j][0]; // Neighbor index

                        const Eigen::Vector3d& BCflg_q = point_list[q].BC; // Boundary condition flags for neighbor q
                        const Eigen::Vector3d& DOF_q = point_list[q].DOF;  // DOF for neighbor q

                        // Loop over each dimension (x, y, z) for neighbor q
                        for (int qq = 0; qq < PD; ++qq) {
                            if (BCflg_q(qq) == 1) { // Check if this DOF for neighbor q is active
                                int col = static_cast<int>(DOF_q(qq)) - 1; // Convert to 0-based indexing

                                if (col >= 0 && col < DOFs) { // Ensure column index is valid
                                    // Add off-diagonal contribution if available
                                    if (pp < K_P.rows() && qq < K_P.cols()) {
                                        triplets.emplace_back(row, col, K_P(pp, qq));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Construct the sparse matrix using the triplets collected
    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}

#endif // ASSEMBLE_H