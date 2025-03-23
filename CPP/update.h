#ifndef UPDATE_H
#define UPDATE_H

#include <vector>
#include <Eigen/Dense>
#include "Points.h"
#include "Neighbour.h"

// === APPLY PRESCRIBED DISPLACEMENT === //
std::vector<Points> update_prescribed(const std::vector<Points>& points, double LF, int PD, double delta) {
    std::vector<Points> updated_points = points;
    Eigen::Matrix3d FF = Compute_FF(PD, LF * 0.25, "EXP");  // Use 0.25 from your d parameter

    for (auto& point : updated_points) {
        if (point.Flag == "Patch") {
            // Left patches stay fixed
            continue;
        }
        else if (point.Flag == "RightPatch") {
            // Right patches move according to deformation gradient
            point.x = FF * point.X;
        }
        else if (point.Flag == "Point") {
            // Only apply prescribed displacement for constrained DOFs
            if (PD == 1) {
                if (point.BC[0] == 0) {
                    point.x = FF * point.X;
                }
            }
        }
    }
    return updated_points;
}

// === APPLY SOLVER DISPLACEMENT === //
std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx, int PD, double delta) {
    // Store previous positions for all points before updating
    for (auto& point : point_list) {
        point.x_prev = point.x;
    }

    // Apply displacements only to unconstrained points
    for (auto& point : point_list) {
        if (point.Flag == "Patch" || point.Flag == "RightPatch") {
            // Patches (both left and right) don't move from solver
            continue;
        }

        // Regular points update from solver - only for free DOFs (BC==1)
        for (int dim = 0; dim < PD; ++dim) {
            if (point.BC[dim] == 1) {
                int dof_idx = point.DOF[dim] - 1;
                if (dof_idx >= 0 && dof_idx < dx.size()) {
                    point.x[dim] += dx(dof_idx);
                }
            }
        }
    }

    return point_list;
}

#endif // UPDATE_H