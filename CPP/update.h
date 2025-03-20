#ifndef UPDATE_H
#define UPDATE_H

#include <vector>
#include <Eigen/Dense>
#include "Points.h"
#include "Neighbour.h"

std::vector<Points> update_prescribed(const std::vector<Points>& points, double LF, int PD, double delta) {
    std::vector<Points> updated_points = points;
    int free_dof_count = 0;
    int fixed_dof_count = 0;

    for (auto& point : updated_points) {
        // Update prescribed boundary conditions
        for (int dim = 0; dim < PD; ++dim) {
            if (point.BC[dim] == 0) { // Fixed DOF
                point.x[dim] = point.X[dim] + LF * point.BCval[dim];
                fixed_dof_count++;
            } else { // Free DOF
                free_dof_count++;
            }
        }
    }

    std::cout << "Number of free DOFs after BC update: " << free_dof_count << std::endl;
    std::cout << "Number of fixed DOFs after BC update: " << fixed_dof_count << std::endl;

    // Check if we have any free DOFs
    if (free_dof_count == 0) {
        std::cout << "WARNING: No free DOFs found. Simulation will not progress!" << std::endl;
    }

    return updated_points;
}

std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx, int PD, double delta) {
    // Safety check
    if (dx.size() == 0) {
        std::cout << "WARNING: Empty displacement vector in update_displaced!" << std::endl;
        return point_list;
    }

    int updated_count = 0;
    double max_disp = 0.0;

    for (auto& point : point_list) {
        for (int dim = 0; dim < PD; ++dim) {
            if (point.BC[dim] == 1) { // Free DOF
                int dof_idx = point.DOF[dim] - 1;
                if (dof_idx >= 0 && dof_idx < dx.size()) {
                    double disp = dx(dof_idx);
                    point.x[dim] += disp;
                    max_disp = std::max(max_disp, std::abs(disp));
                    updated_count++;
                }
            }
        }
    }

    std::cout << "Updated " << updated_count << " DOFs with max displacement: " << max_disp << std::endl;

    // Regenerate neighbor lists after updating positions
    // generate_neighbour_list(PD, point_list, delta);

    return point_list;
}


#endif // UPDATE_H