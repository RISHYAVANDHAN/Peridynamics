//
// Created by srini on 28/02/2025.
//

#ifndef UPDATE_H
#define UPDATE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Points.h"
#include "Neighbour.h"

std::vector<Points> update_prescribed(std::vector<Points>& point_list, const int LF, const int PD, const double delta) {
    for (auto& i : point_list) {
        for (int j = 0; j < PD; j++) {
            i.x[j] = i.X[j] + LF * i.BCval[j];
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}

std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx, const int PD, const double delta) {
    for (auto& i : point_list) {
        for (int j = 0; j < PD; j++) {
            // Make sure the DOF index is within valid range
            int dof_idx = static_cast<int>(i.DOF(j));
            if (dof_idx >= 0 && dof_idx < dx.size()) {
                i.x[j] = i.X[j] + dx.coeff(dof_idx);
            }
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}

#endif //UPDATE_H