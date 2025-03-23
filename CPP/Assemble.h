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

    return global_energy;
}

// Assemble the global residual vector
inline Eigen::VectorXd assembleResidual(const std::vector<Points>& points, int DOFs, int PD) {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(DOFs);

    for (const auto& p : points) {
        for (int dim = 0; dim < PD; ++dim) {
            if (p.BC[dim] == 1) {  // Free DOF
                int idx = p.DOF[dim] - 1;  // Convert 1-based DOF to 0-based index
                if (idx >= 0 && idx < DOFs) {
                    R(idx) = p.Ra_sum(dim);  // Note: using = instead of += could be an issue!
                }
            }
        }
    }

    std::cout << "Global residual norm: " << R.norm() << std::endl;
    return R;
}

// Assemble the global stiffness matrix
inline Eigen::SparseMatrix<double> assembleStiffness(const std::vector<Points>& points, int DOFs, int PD) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(DOFs * 10);  // Reserve space for efficiency

    // First pass: collect all triplets
    for (size_t p_idx = 0; p_idx < points.size(); ++p_idx) {
        const auto& point = points[p_idx];

        // Skip points that are fully constrained
        bool has_free_dof = false;
        for (int i = 0; i < PD; ++i) {
            if (point.BC[i] == 1) {
                has_free_dof = true;
                break;
            }
        }
        if (!has_free_dof) continue;

        // Self-contribution (diagonal block)
        for (int i = 0; i < PD; ++i) {
            if (point.BC[i] != 1) continue;
            int row = point.DOF[i] - 1;

            for (int j = 0; j < PD; ++j) {
                if (point.BC[j] != 1) continue;
                int col = point.DOF[j] - 1;

                if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                    tripletList.push_back(T(row, col, point.Kab_sum(i, j)));
                }
            }
        }

        // Process 1-neighbor interactions
        for (const auto& neighbor_indices : point.neighbour_list_1N) {
            int neighbor_idx = neighbor_indices[0];
            if (neighbor_idx < 0 || neighbor_idx >= points.size()) continue;

            const auto& neighbor = points[neighbor_idx];

            for (int i = 0; i < PD; ++i) {
                if (point.BC[i] != 1) continue;
                int row = point.DOF[i] - 1;

                for (int j = 0; j < PD; ++j) {
                    if (neighbor.BC[j] != 1) continue;
                    int col = neighbor.DOF[j] - 1;

                    if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                        tripletList.push_back(T(row, col, -point.Kab_1(i, j)));
                    }
                }
            }
        }

        // Process 2-neighbor interactions if applicable
        if (PD >= 2) {
            for (const auto& neighbor_indices : point.neighbour_list_2N) {
                if (neighbor_indices.size() < 2) continue;

                for (int n = 0; n < 2; ++n) {
                    int neighbor_idx = neighbor_indices[n];
                    if (neighbor_idx < 0 || neighbor_idx >= points.size()) continue;

                    const auto& neighbor = points[neighbor_idx];

                    for (int i = 0; i < PD; ++i) {
                        if (point.BC[i] != 1) continue;
                        int row = point.DOF[i] - 1;

                        for (int j = 0; j < PD; ++j) {
                            if (neighbor.BC[j] != 1) continue;
                            int col = neighbor.DOF[j] - 1;

                            if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                                tripletList.push_back(T(row, col, -point.Kab_2(i, j)));
                            }
                        }
                    }
                }
            }
        }

        // Process 3-neighbor interactions if applicable
        if (PD == 3) {
            for (const auto& neighbor_indices : point.neighbour_list_3N) {
                if (neighbor_indices.size() < 3) continue;

                for (int n = 0; n < 3; ++n) {
                    int neighbor_idx = neighbor_indices[n];
                    if (neighbor_idx < 0 || neighbor_idx >= points.size()) continue;

                    const auto& neighbor = points[neighbor_idx];

                    for (int i = 0; i < PD; ++i) {
                        if (point.BC[i] != 1) continue;
                        int row = point.DOF[i] - 1;

                        for (int j = 0; j < PD; ++j) {
                            if (neighbor.BC[j] != 1) continue;
                            int col = neighbor.DOF[j] - 1;

                            if (row >= 0 && col >= 0 && row < DOFs && col < DOFs) {
                                tripletList.push_back(T(row, col, -point.Kab_3(i, j)));
                            }
                        }
                    }
                }
            }
        }
    }

    // Create sparse matrix from triplets
    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << "Assembled Stiffness Matrix (Sparse Format):\n";
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            std::cout << "K(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
        }
    }


    return K;
}

#endif // ASSEMBLE_H