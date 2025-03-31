#ifndef UPDATE_1D_H
#define UPDATE_1D_H

#include <vector>
#include <Eigen/Dense>
#include "Points_1D.h"

// Update positions based on prescribed BCs (1D)
std::vector<Points> update_prescribed(const std::vector<Points>& points, double LF) {
    std::vector<Points> updated_points = points;
    int free_dof_count = 0;
    int fixed_dof_count = 0;

    for (auto& point : updated_points) {
        if (point.BC == 0) {  // Fixed DOF
            point.x = point.X + LF * point.BCval;
            fixed_dof_count++;
            std::cout << "BC x = " << point.x << std::endl;
        } else {  // Free DOF
            point.x += 1e-4;  // Small perturbation
            free_dof_count++;
            std::cout << "Free x = " << point.x << std::endl;
        }
    }

    std::cout << "Free DOFs: " << free_dof_count << "\n";
    std::cout << "Fixed DOFs: " << fixed_dof_count << "\n";

    if (free_dof_count == 0) {
        std::cout << "WARNING: No free DOFs found. Simulation will not progress!\n";
    }

    return updated_points;
}

// Update positions using displacement vector dx (1D)
std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx) {
    if (dx.size() == 0) {
        std::cout << "WARNING: Empty displacement vector in update_displaced!\n";
        return point_list;
    }

    int updated_count = 0;
    double max_disp = 0.0;

    for (auto& point : point_list) {
        if (point.BC == 1) {  // Free DOF
            int dof_idx = point.DOF - 1;
            if (dof_idx >= 0 && dof_idx < dx.size()) {
                double disp = dx(dof_idx);
                point.x += disp;
                max_disp = std::max(max_disp, std::abs(disp));
                updated_count++;
            }
        }
    }

    std::cout << "Updated " << updated_count << " DOFs. Max displacement: " << max_disp << "\n";
    return point_list;
}

#endif // UPDATE_1D_H
