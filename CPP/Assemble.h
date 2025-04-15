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
                    R(idx) += p.Ra_sum(dim);  // FIXED: Using += instead of = to accumulate residuals
                }
            }
        }
    }

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

        // Process 1-neighbor interactions (simplified for 1D)
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
    }

    // Create sparse matrix from triplets
    Eigen::SparseMatrix<double> K(DOFs, DOFs);
    K.setFromTriplets(tripletList.begin(), tripletList.end());

    // Verify matrix is not empty
    if (K.nonZeros() == 0) {
        std::cerr << "WARNING: Empty stiffness matrix created!" << std::endl;
    }

    return K;
}

#endif // ASSEMBLE_H