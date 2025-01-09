//
// Created by srini on 08/01/2025.
//

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"

// Calculating point-wise residual and tangent
void cal_pw_residual_tangent(std::vector<Points>& point_list, double NN, double C1)
{
    Points points;
    double Xi, xi, L, l, s, Psi;
    for (size_t i = 0; i < point_list.size(); ++i)
    {
        for (size_t j = 0; j < points.neighbour_list[i].size(); ++j)
        {
            Xi = points.X[i] - points.X[j];
            xi = points.x[i] - points.x[j];
            L = std::sqrt(Xi * Xi);
            l = std::sqrt(xi * Xi);
            s = (1.0 / NN) * (std::pow((l/L), NN) - 1);
            //CC!
            Psi = 0.5 * C1 * L * (s*s);
            //Residual
            //Tangent
        }
    }
}




#endif //SOLVER_H
