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

std::vector<Points> update_prescribed(std::vector<Points>& point_list, const int LF, const int PD, const double delta)
{
    for (auto& i : point_list)
    {
        for (int j = 0; j < PD; j++)  // Fix: Start from 0 and use < instead of <=
        {
            i.x[j] = i.X[j] + LF * i.BCval[j];  // Fix: Use 0-based indexing
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}

std::vector<Points> update_displaced(std::vector<Points>& point_list, const Eigen::VectorXd& dx, const int PD, const double delta)
{
    for (auto& i : point_list)
    {
        for (int j = 0; j < PD; j++)  // Fix: Start from 0 and use < instead of <=
        {
            i.x[j] = i.X[j] + dx.coeff(i.DOF(j));  // Fix: Use 0-based indexing
        }
    }
    generate_neighbour_list(PD, point_list, delta);
    return point_list;
}
#endif //UPDATE_H
