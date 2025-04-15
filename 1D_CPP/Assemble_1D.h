#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points_1D.h"

// Function to compute the total energy of the system
inline double assembleEnergy(const std::vector<Points>& point_list) {
    double Psi = 0.0;
    for (const auto& point : point_list) {
        std::cout<< point.volume << std::endl;
        Psi += point.volume * point.psi; // Sum up energy contributions from each point
    }

    return Psi;
}

// Function to assemble the residual vector
inline Eigen::VectorXd assembleResidual(const std::vector<Points>& point_list, const int DOFs) {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs); // Initialize residual vector

    for (const auto& point : point_list) {
        const double& R_P = point.Ra_sum;
        int PD = 1; // Problem dimension (PD)

        for (int ii = 0; ii < PD; ii++) {
            if (point.BC == 1.0) { // Check if DOF is active
                int dof_idx = static_cast<int>(point.DOF); // Get DOF index

                if (dof_idx >= 0 && dof_idx < DOFs) { // Ensure valid DOF index
                    R(dof_idx) += R_P; // Add contribution to residual
                }
                // If the DOF index is invalid, skip without printing a warning
            }
        }
    }

    return R; // Return the assembled residual vector
}


// Function to assemble the stiffness matrix
inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& point_list, int DOFs) {
    std::vector<Eigen::Triplet<double>> triplets;

    for (size_t p = 0; p < point_list.size(); ++p) {
        const double& K_P = point_list[p].Kab_sum;
        const int BCflg_p = point_list[p].BC;
        const int DOF_p   = point_list[p].DOF;

        if (BCflg_p == 1) {
            int row = DOF_p - 1;
            if (row >= 0 && row < DOFs) {
                triplets.emplace_back(row, row, K_P);  // Diagonal contribution
            }
        }
    }

    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}


#endif // ASSEMBLE_H