#ifndef UPDATE_H
#define UPDATE_H

#include <vector>
#include <Eigen/Dense>
#include "Points.h"
#include "Neighbour.h"

// Apply prescribed displacement (simplified for 1D)
inline std::vector<Points> update_prescribed(const std::vector<Points>& points, double LF, int PD, double d) {
    std::vector<Points> updated_points = points;

    // For 1D, we'll use a simple scaling factor
    double scale_factor = 1.0 + (LF * d); // Based on your d=0.25

    for (auto& point : updated_points) {
        if (point.Flag == "Patch") {
            // Left patches stay fixed
            continue;
        }
        else if (point.Flag == "RightPatch") {
            // Right patches move according to prescribed displacement
            point.x[0] = scale_factor * point.X[0];
        }
    }

    return updated_points;
}

// Apply displacement increment from solver
inline std::vector<Points> update_displaced(std::vector<Points>& points, const Eigen::VectorXd& dx, int PD, double delta) {
    // Store previous positions
    for (auto& point : points) {
        point.x_prev = point.x;
    }

    // Apply displacements only to points with free DOFs
    for (auto& point : points) {
        // Skip patches
        if (point.Flag == "Patch" || point.Flag == "RightPatch") {
            continue;
        }

        // For 1D, we only need to update the x-coordinate for points with BC[0]=1
        if (point.BC[0] == 1) {
            int dof_idx = point.DOF[0] - 1;
            if (dof_idx >= 0 && dof_idx < dx.size()) {
                point.x[0] += dx(dof_idx);
            }
        }
    }

    return points;
}

#endif // UPDATE_H