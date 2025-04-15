#ifndef SOLVER_1D_H
#define SOLVER_1D_H

#include <iostream>
#include <vector>
#include <cmath>
#include "Points_1D.h"
#include "mesh_1D.h"
#include "Neighbour_1D.h"

// 1D residual and tangent stiffness calculation
void cal_pw_residual_tangent(std::vector<Points>& point_list, double NN, double C1, double horizon)
{
    for (size_t i = 0; i < point_list.size(); ++i)
    {
        double Ra_i = 0.0;
        double Kab_i = 0.0;

        for (const auto& neighbor : point_list[i].neighbour_list)
        {
            int j = neighbor[0];  // Only one neighbor index in 1D

            double Xi = point_list[i].X - point_list[j].X;
            double xi = point_list[i].x - point_list[j].x;

            double L = std::abs(Xi);
            double l = std::abs(xi);

            if (L == 0.0 || l == 0.0) continue; // avoid division by zero

            double s = (1.0 / NN) * (std::pow((l / L), NN) - 1.0);
            double C11 = 4.0 * C1 * std::abs(1.0 - std::pow((L / horizon), 2));

            // Residual
            double dL_dXi = (Xi > 0) ? 1.0 : -1.0;
            double dl_dXi = (xi > 0) ? 1.0 : -1.0;

            double dL = dL_dXi;
            double dl = dl_dXi;

            double ds_dXi = (std::pow(l, NN - 1) / std::pow(L, NN)) * dl -
                            (NN * std::pow(l, NN) / std::pow(L, NN + 1)) * dL;

            double res = C11 * (L * s * ds_dXi + 0.5 * s * s * dL);
            Ra_i += res;

            // Tangent stiffness
            double d2L_dXi2 = 0.0;
            double d2s_dXi2 = (NN - 1) * std::pow(l, NN - 2) / std::pow(L, NN) * (dl * dl) -
                              2 * NN * std::pow(l, NN) / std::pow(L, NN + 1) * (dl * dL) -
                              NN * (NN + 1) * std::pow(l, NN) / std::pow(L, NN + 2) * (dL * dL);

            double K = C11 * (ds_dXi * ds_dXi * dL + L * d2s_dXi2 +
                             0.5 * (2 * ds_dXi * ds_dXi * dL + ds_dXi * ds_dXi * d2L_dXi2));
            Kab_i += K;
        }

        // Save result into the point's data
        point_list[i].Ra_1 = Ra_i;
        point_list[i].Ra_sum = Ra_i;

        point_list[i].Kab_1 = Kab_i;
        point_list[i].Kab_sum = Kab_i;
    }
}

#endif // SOLVER_1D_H
