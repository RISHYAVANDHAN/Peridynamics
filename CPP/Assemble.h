#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"
#include <iostream>

// Assemble the global energy
inline double assembleEnergy(const std::vector<Points>& points) {
    double global_energy = 0.0;

    for (const auto& p : points) {
        global_energy += p.psi;  // Accumulate energy from each point
    }

    // Debug output
    std::cout << "Total energy calculated in assembleEnergy: " << global_energy << std::endl;
    return global_energy;
}

// Assemble the global residual vector
inline Eigen::VectorXd assembleResidual(const std::vector<Points>& points, int DOFs, int PD) {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs);

    // Debug count
    int free_dofs_count = 0;

    for (size_t p_idx = 0; p_idx < points.size(); ++p_idx) {
        const Points& p = points[p_idx];

        // Only consider free DOFs (BC == 1)
        for (int dim = 0; dim < PD; ++dim) {
            if (p.BC[dim] == 1) {  // Free DOF
                int idx = p.DOF[dim] - 1;  // Convert 1-based DOF to 0-based index
                if (idx >= 0 && idx < DOFs) {
                    R(idx) += p.Ra_sum(dim);  // Accumulate residual
                    free_dofs_count++;
                }
            }
        }
    }

    // Debug output
    std::cout << "Free DOFs count: " << free_dofs_count << ", Residual norm in assembleResidual: " << R.norm() << std::endl;
    return R;
}

// Assemble the global stiffness matrix
inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& points, int DOFs, int PD) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(DOFs * 10);  // Reserve space for non-zero entries

    // Debug count
    int stiffness_entries = 0;

    // For each point, add its contribution to the global stiffness matrix
    for (size_t i = 0; i < points.size(); ++i) {
        const Points& point_i = points[i];

        // Add point i's diagonal contribution to itself (self-stiffness)
        for (int dim_i = 0; dim_i < PD; ++dim_i) {
            if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
            int row = point_i.DOF[dim_i] - 1;  // Convert 1-based DOF to 0-based index

            for (int dim_j = 0; dim_j < PD; ++dim_j) {
                if (point_i.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                int col = point_i.DOF[dim_j] - 1;

                if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                    triplets.emplace_back(row, col, point_i.Kab_sum(dim_i, dim_j));
                    stiffness_entries++;
                }
            }
        }

        // 1N interactions
        for (const auto& neighbor : point_i.neighbour_list_1N) {
            size_t j = neighbor[0];
            if (j >= points.size()) continue;

            const Points& point_j = points[j];

            // Calculate off-diagonal contribution (i to j)
            for (int dim_i = 0; dim_i < PD; ++dim_i) {
                if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                int row = point_i.DOF[dim_i] - 1;

                for (int dim_j = 0; dim_j < PD; ++dim_j) {
                    if (point_j.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                    int col = point_j.DOF[dim_j] - 1;

                    if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                        // For 1N neighbors, use the negative of Kab_1
                        triplets.emplace_back(row, col, -point_i.Kab_1(dim_i, dim_j));
                        stiffness_entries++;
                    }
                }
            }
        }

        // 2N interactions (only for PD >= 2)
        if (PD >= 2) {
            for (const auto& neighbor_pair : point_i.neighbour_list_2N) {
                size_t j1 = neighbor_pair[0];
                size_t j2 = neighbor_pair[1];

                if (j1 >= points.size() || j2 >= points.size()) continue;

                const Points& point_j1 = points[j1];
                const Points& point_j2 = points[j2];

                // Calculate off-diagonal contribution for j1
                for (int dim_i = 0; dim_i < PD; ++dim_i) {
                    if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                    int row = point_i.DOF[dim_i] - 1;

                    for (int dim_j = 0; dim_j < PD; ++dim_j) {
                        if (point_j1.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                        int col = point_j1.DOF[dim_j] - 1;

                        if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                            // For 2N neighbors, use the negative of Kab_2
                            triplets.emplace_back(row, col, -point_i.Kab_2(dim_i, dim_j));
                            stiffness_entries++;
                        }
                    }
                }

                // Calculate off-diagonal contribution for j2
                for (int dim_i = 0; dim_i < PD; ++dim_i) {
                    if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                    int row = point_i.DOF[dim_i] - 1;

                    for (int dim_j = 0; dim_j < PD; ++dim_j) {
                        if (point_j2.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                        int col = point_j2.DOF[dim_j] - 1;

                        if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                            // For 2N neighbors, use the negative of Kab_2
                            triplets.emplace_back(row, col, -point_i.Kab_2(dim_i, dim_j));
                            stiffness_entries++;
                        }
                    }
                }
            }
        }

        // 3N interactions (only for PD == 3)
        if (PD == 3) {
            for (const auto& neighbor_triplet : point_i.neighbour_list_3N) {
                size_t j1 = neighbor_triplet[0];
                size_t j2 = neighbor_triplet[1];
                size_t j3 = neighbor_triplet[2];

                if (j1 >= points.size() || j2 >= points.size() || j3 >= points.size()) continue;

                const Points& point_j1 = points[j1];
                const Points& point_j2 = points[j2];
                const Points& point_j3 = points[j3];

                // Calculate off-diagonal contribution for j1
                for (int dim_i = 0; dim_i < PD; ++dim_i) {
                    if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                    int row = point_i.DOF[dim_i] - 1;

                    for (int dim_j = 0; dim_j < PD; ++dim_j) {
                        if (point_j1.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                        int col = point_j1.DOF[dim_j] - 1;

                        if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                            // For 3N neighbors, use the negative of Kab_3
                            triplets.emplace_back(row, col, -point_i.Kab_3(dim_i, dim_j));
                            stiffness_entries++;
                        }
                    }
                }

                // Calculate off-diagonal contribution for j2
                for (int dim_i = 0; dim_i < PD; ++dim_i) {
                    if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                    int row = point_i.DOF[dim_i] - 1;

                    for (int dim_j = 0; dim_j < PD; ++dim_j) {
                        if (point_j2.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                        int col = point_j2.DOF[dim_j] - 1;

                        if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                            // For 3N neighbors, use the negative of Kab_3
                            triplets.emplace_back(row, col, -point_i.Kab_3(dim_i, dim_j));
                            stiffness_entries++;
                        }
                    }
                }

                // Calculate off-diagonal contribution for j3
                for (int dim_i = 0; dim_i < PD; ++dim_i) {
                    if (point_i.BC[dim_i] != 1) continue;  // Skip constrained DOFs
                    int row = point_i.DOF[dim_i] - 1;

                    for (int dim_j = 0; dim_j < PD; ++dim_j) {
                        if (point_j3.BC[dim_j] != 1) continue;  // Skip constrained DOFs
                        int col = point_j3.DOF[dim_j] - 1;

                        if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                            // For 3N neighbors, use the negative of Kab_3
                            triplets.emplace_back(row, col, -point_i.Kab_3(dim_i, dim_j));
                            stiffness_entries++;
                        }
                    }
                }
            }
        }
    }

    // Debug output
    std::cout << "Added " << stiffness_entries << " entries to stiffness matrix" << std::endl;

    // Build the sparse stiffness matrix
    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();

    return K;
}

#endif // ASSEMBLE_H