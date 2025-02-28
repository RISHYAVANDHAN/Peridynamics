//
// Created by srini on 26/02/2025.
//

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"


// Assembly function for total energy
inline double assembleEnergy(const std::vector<Points>& point_list) {
    double Psi = 0.0;

    for (const auto& point : point_list) {
        // Using psi and volume from your Points class
        Psi += point.volume * point.psi;
    }

    return Psi;
}

// Assembly function for residual
inline Eigen::VectorXd assembleResidual(const std::vector<Points>& point_list, const int DOFs) {

    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs);
    for (size_t i = 0; i < point_list.size(); i++) {
        // Get residual for this point

        Eigen::VectorXd R_P = point_list[i].Ra_sum;

        // Get dimension from BC vector
        int PD = point_list[i].BC.size();

        for (int ii = 0; ii < PD; ii++) {
            if (point_list[i].BC(ii) == 1) {
                int dof_idx = static_cast<int>(point_list[i].DOF(ii));
                if (dof_idx > 0 && dof_idx <= DOFs) {
                    R(dof_idx-1) += R_P(ii);
                }
            }
        }
    }

    return R;
}

// Assembly function for stiffness matrix
inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& point_list, int DOFs) {
    // Triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    // Estimate size - adjust as needed
    triplets.reserve(DOFs * 20);

    // Dimension of the problem
    int PD = point_list[0].BC.size();

    for (size_t p = 0; p < point_list.size(); p++) {
        // Process 1-neighbor interactions
        for (int pp = 0; pp < PD; pp++) {
            if (point_list[p].BC(pp) == 1) {
                int dof_p = static_cast<int>(point_list[p].DOF(pp));

                // Process all neighbors including the point itself
                std::vector<int> all_neighbors;
                for (const auto& nbrData : point_list[p].neighbour_list_1N) {
                    all_neighbors.push_back(nbrData[0]);
                }
                all_neighbors.push_back(p); // Add self (as done in MATLAB)

                for (const auto& q : all_neighbors) {
                    for (int qq = 0; qq < PD; qq++) {
                        if (point_list[q].BC(qq) == 1) {
                            int dof_q = static_cast<int>(point_list[q].DOF(qq));

                            // Extract stiffness value
                            double k_value = 0.0;
                            if (PD == 2) {
                                // For 2D
                                if (pp == 0 && qq == 0) k_value = point_list[p].Kab_1(0, 0);
                                else if (pp == 0 && qq == 1) k_value = point_list[p].Kab_1(0, 1);
                                else if (pp == 1 && qq == 0) k_value = point_list[p].Kab_1(1, 0);
                                else if (pp == 1 && qq == 1) k_value = point_list[p].Kab_1(1, 1);
                            }
                            else if (PD == 3) {
                                // For 3D
                                k_value = point_list[p].Kab_1(pp, qq);
                            }

                            // Add to triplets
                            if (dof_p > 0 && dof_q > 0) {
                                triplets.emplace_back(dof_p-1, dof_q-1, k_value);
                            }
                        }
                    }
                }

                // Process 2-neighbor interactions
                for (const auto& nbrData : point_list[p].neighbour_list_2N) {
                    int q1 = nbrData[0];
                    int q2 = nbrData[1];

                    for (int qq = 0; qq < PD; qq++) {
                        if (point_list[q1].BC(qq) == 1) {
                            int dof_q = static_cast<int>(point_list[q1].DOF(qq));

                            // Extract stiffness value for 2-neighbor interaction
                            double k_value = point_list[p].Kab_2(pp, qq);

                            // Add to triplets
                            if (dof_p > 0 && dof_q > 0) {
                                triplets.emplace_back(dof_p-1, dof_q-1, k_value);
                            }
                        }

                        if (point_list[q2].BC(qq) == 1) {
                            int dof_q = static_cast<int>(point_list[q2].DOF(qq));

                            // Extract stiffness value for 2-neighbor interaction
                            double k_value = point_list[p].Kab_2(pp, qq);

                            // Add to triplets
                            if (dof_p > 0 && dof_q > 0) {
                                triplets.emplace_back(dof_p-1, dof_q-1, k_value);
                            }
                        }
                    }
                }

                // Process 3-neighbor interactions (if PD == 3)
                if (PD == 3) {
                    for (const auto& nbrData : point_list[p].neighbour_list_3N) {
                        int q1 = nbrData[0];
                        int q2 = nbrData[1];
                        int q3 = nbrData[2];

                        for (int qq = 0; qq < PD; qq++) {
                            // Process each neighbor
                            const std::vector<int> nbrs = {q1, q2, q3};
                            for (const auto& q : nbrs) {
                                if (point_list[q].BC(qq) == 1) {
                                    int dof_q = static_cast<int>(point_list[q].DOF(qq));

                                    // Extract stiffness value for 3-neighbor interaction
                                    double k_value = point_list[p].Kab_3(pp, qq);

                                    // Add to triplets
                                    if (dof_p > 0 && dof_q > 0) {
                                        triplets.emplace_back(dof_p-1, dof_q-1, k_value);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Create sparse matrix from triplets
    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(triplets.begin(), triplets.end());

    return K;
}

#endif // ASSEMBLE_H