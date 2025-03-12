#ifndef UPDATE_H
#define UPDATE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"
#include "Neighbour.h"

std::vector<Points> update_prescribed(std::vector<Points>& point_list, const double LF, const int PD, const double delta) {
    for (auto& i : point_list) {
        for (int j = 0; j < PD; j++) {
            if (j < i.x.size() && j < i.X.size() && j < i.BCval.size()) {
                i.x[j] = i.X[j] + LF * i.BCval[j];
            } else {
                std::cerr << "Error: Vector index out of range in update_prescribed" << std::endl;
            }
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}

std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx, const int PD, const double delta) {
    for (auto& i : point_list) {
        for (int j = 0; j < PD; j++) {
            if (i.BC(j) == 1) {  // Only update if this is an active DOF
                int dof_idx = static_cast<int>(i.DOF(j));
                if (dof_idx > 0 && dof_idx <= dx.size()) {
                    // Use the displacement increment to update the current position
                    i.x[j] += dx.coeff(dof_idx - 1);  // Not replacing but incrementing
                }
            }
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}

#endif // UPDATE_H