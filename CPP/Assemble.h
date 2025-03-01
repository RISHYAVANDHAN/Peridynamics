//
// Created by srini on 26/02/2025.
//

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"

inline double assembleEnergy(const std::vector<Points>& point_list) {
    double Psi = 0.0;

    for (const auto& point : point_list) {
        Psi += point.volume * point.psi;
    }

    return Psi;
}

inline Eigen::VectorXd assembleResidual(const std::vector<Points>& point_list, const int DOFs) {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs);

    for (size_t i = 0; i < point_list.size(); i++) {
        Eigen::Vector3d R_P = point_list[i].Ra_sum;

        for (int ii = 0; ii < 3; ii++) {
            if (point_list[i].BC(ii) == 1) {
                int dof_idx = static_cast<int>(point_list[i].DOF(ii));
                // Check if dof_idx is zero-indexed (0 to DOFs-1)
                if (dof_idx >= 0 && dof_idx < DOFs) {
                    R(dof_idx) += R_P(ii);
                }
            }
        }
    }

    return R;
}


inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& point_list, int DOFs, int PD) {
    std::vector<Eigen::Triplet<double>> triplets;  // For sparse matrix storage

    // Iterate over each point in the point list
    for (auto p = 0; p < point_list.size(); ++p) {
        Eigen::MatrixXd K_P = point_list[p].Kab_sum;  // Total stiffness matrix for point p
        Eigen::Vector3d BCflg_p = point_list[p].BC;   // Boundary condition flags for point p
        Eigen::Vector3d DOF_p = point_list[p].DOF;    // DOF for point p

        // Loop over each dimension (x, y, z) depending on PD (1D, 2D or 3D)
        for (int pp = 0; pp < PD; ++pp) {
            if (BCflg_p(pp) == 1) {  // Check if this DOF is active
                int dof_p = static_cast<int>(DOF_p(pp));  // Global DOF index for point p

                // Make sure dof_p is valid (zero-indexed)
                if (dof_p >= 0 && dof_p < DOFs) {
                    // Copy neighbour list (2D vector) and append `p` to its own list of neighbors
                    std::vector<std::vector<int>> nbrL = point_list[p].neighbour_list;
                    nbrL.push_back({p});  // Append the current point to its neighbors

                    // Loop over each neighbor (including self)
                    for (size_t q_idx = 0; q_idx < nbrL.size(); ++q_idx) {
                        int q = nbrL[q_idx][0];  // Get the neighbor index

                        // Get the boundary condition flags and DOFs for neighbor `q`
                        Eigen::Vector3d BCflg_q = point_list[q].BC;
                        Eigen::Vector3d DOF_q = point_list[q].DOF;

                        // Make sure we're not attempting to extract a block larger than the matrix dimensions
                        if (pp < K_P.rows() && q_idx < K_P.cols()) {
                            // Extract the stiffness contribution K_PQ from the total stiffness matrix
                            double K_value = 0.0;

                            // For 1D, simply get the direct value
                            if (PD == 1) {
                                K_value = K_P(pp, q_idx);

                                if (BCflg_q(0) == 1) {  // Check if this DOF for neighbor `q` is active
                                    int dof_q = static_cast<int>(DOF_q(0));  // Global DOF index for neighbor `q`

                                    // Add to the triplet list (for sparse matrix construction)
                                    if (dof_q >= 0 && dof_q < DOFs) {
                                        triplets.emplace_back(dof_p, dof_q, K_value);
                                    }
                                }
                            }
                            // For higher dimensions, we need to extract blocks
                            else {
                                // Loop over each dimension (x, y, z) for neighbor `q`
                                for (int qq = 0; qq < PD; ++qq) {
                                    if (BCflg_q(qq) == 1) {  // Check if this DOF for neighbor `q` is active
                                        int dof_q = static_cast<int>(DOF_q(qq));  // Global DOF index for neighbor `q`

                                        // Safely get the value from K_P
                                        if (pp < K_P.rows() && q_idx * PD + qq < K_P.cols()) {
                                            K_value = K_P(pp, q_idx * PD + qq);

                                            // Add to the triplet list (for sparse matrix construction)
                                            if (dof_q >= 0 && dof_q < DOFs) {
                                                triplets.emplace_back(dof_p, dof_q, K_value);
                                            }
                                        }
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