#pragma once
#include <fstream>
#include <memory>
#include "Grid.hpp"
#include "hyperdual.h"
#include <cmath>
#include <cstdlib>
#include "Eigen/Dense"
#include "Eigen/Sparse"
using EMx = Eigen::MatrixXd;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////             Norm Calculation               ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double Norm2(std::array<double, 2> &p)
{
    return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}
double Norm2(double x, double y)
{
    return std::sqrt(x * x + y * y);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////                 Tangent                    ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Cal_K2(double nn, int idx, int deb_idx, double C2, double V2, double Psi_1_x, double Psi_1_y, double Psi_1_z, double Psi_2_x, double Psi_2_y, double Psi_2_z, double psi_1_x, double psi_1_y, double psi_1_z, double psi_2_x, double psi_2_y, double psi_2_z, double *AA2I, double *AA2J)
{
    if (nn != 0)
    {
        hyperdual a, s, psi, a_HD[3];
        double A;
        double dpsi[3];
        double ddpsi[3][3];
        double A_tmp[3] = {0.0, 0.0, Psi_1_x * Psi_2_y - Psi_1_y * Psi_2_x};
        double a_tmp[3] = {0.0, 0.0, psi_1_x * psi_2_y - psi_1_y * psi_2_x};
        a_HD[0] = hyperdual(a_tmp[0], 1.0, 1.0, 0.0);
        a_HD[1] = hyperdual(a_tmp[1], 0.0, 0.0, 0.0);
        a_HD[2] = hyperdual(a_tmp[2], 0.0, 0.0, 0.0);
        A = std::fabs(A_tmp[2]);
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[0][0] = psi.eps1eps2();
        ///////////////////////////////////////////
        a_HD[0].f1 = 0.0;
        a_HD[0].f2 = 0.0;
        a_HD[1].f1 = 1.0;
        a_HD[1].f2 = 1.0;
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        dpsi[1] = psi.eps1();
        ddpsi[1][1] = psi.eps1eps2();
        ///////////////////////////////////////////////
        a_HD[1].f1 = 0.0;
        a_HD[1].f2 = 0.0;
        a_HD[2].f1 = 1.0;
        a_HD[2].f2 = 1.0;
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        dpsi[2] = psi.eps1();
        ddpsi[2][2] = psi.eps1eps2();
        ////////////////////////////////////////////
        a_HD[0].f1 = 1.0;
        a_HD[1].f2 = 1.0;
        a_HD[2].f1 = 0.0;
        a_HD[2].f2 = 0.0;
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[0][1] = psi.eps1eps2();
        ddpsi[1][0] = ddpsi[0][1];
        ////////////////////////////////////////////
        a_HD[1].f2 = 0.0;
        a_HD[2].f2 = 1.0;
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[0][2] = psi.eps1eps2();
        ddpsi[2][0] = ddpsi[0][2];
        ///////////////////////////////////////////////
        a_HD[0].f1 = 0.0;
        a_HD[1].f1 = 1.0;
        a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
        s = (1.0 / nn) * (pow(a / A, nn) - 1);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[1][2] = psi.eps1eps2();
        ddpsi[2][1] = ddpsi[1][2];
        // CPP
        double dpsi_e[3][3] = {{0, dpsi[2], -dpsi[1]}, {-dpsi[2], 0, dpsi[0]}, {dpsi[1], -dpsi[0], 0}};
        double AA1[3][3], AA2[3][3];

        AA1[0][0] = -ddpsi[2][2] * psi_2_y * (2 * psi_1_y - 2 * psi_2_y);
        AA1[0][1] = +ddpsi[2][2] * psi_2_x * (2 * psi_1_y - 2 * psi_2_y);
        AA1[0][2] = -ddpsi[2][1] * psi_2_x * (2 * psi_1_y - 2 * psi_2_y);
        AA1[1][0] = +ddpsi[2][2] * psi_2_y * (2 * psi_1_x - 2 * psi_2_x);
        AA1[1][1] = -ddpsi[2][2] * psi_2_x * (2 * psi_1_x - 2 * psi_2_x);
        AA2[0][0] = +ddpsi[2][2] * psi_1_y * (2 * psi_1_y - 2 * psi_2_y);
        AA2[0][1] = -ddpsi[2][2] * psi_1_x * (2 * psi_1_y - 2 * psi_2_y);
        AA2[0][2] = +ddpsi[2][1] * psi_1_x * (2 * psi_1_y - 2 * psi_2_y);
        AA2[1][0] = -ddpsi[2][2] * psi_1_y * (2 * psi_1_x - 2 * psi_2_x);
        AA2[1][1] = +ddpsi[2][2] * psi_1_x * (2 * psi_1_x - 2 * psi_2_x);
        /////////////////////////////////////////////////////////////////

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                AA2I[i * 2 + j] = V2 * (-2.0 * dpsi_e[i][j] + AA1[i][j]);
                AA2J[i * 2 + j] = V2 * (2.0 * dpsi_e[i][j] + AA2[i][j]);
            }
        }
    }
    else
    {
        hyperdual a, s, psi, a_HD[3];
        double A;
        double dpsi[3];
        double ddpsi[3][3];
        double A_tmp[3] = {0.0, 0.0, Psi_1_x * Psi_2_y - Psi_1_y * Psi_2_x};
        double a_tmp[3] = {0.0, 0.0, psi_1_x * psi_2_y - psi_1_y * psi_2_x};
        a_HD[0] = hyperdual(a_tmp[0], 1.0, 1.0, 0.0);
        a_HD[1] = hyperdual(a_tmp[1], 0.0, 0.0, 0.0);
        a_HD[2] = hyperdual(a_tmp[2], 0.0, 0.0, 0.0);
        A = std::fabs(A_tmp[2]);
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        dpsi[0] = psi.eps1();
        ddpsi[0][0] = psi.eps1eps2();
        //////////////////////////////////////////
        a_HD[0].f1 = 0.0;
        a_HD[0].f2 = 0.0;
        a_HD[1].f1 = 1.0;
        a_HD[1].f2 = 1.0;
        a = sqrt(pow(a_HD[0], 2.0) + pow(a_HD[1], 2.0) + pow(a_HD[2], 2.0));
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        dpsi[1] = psi.eps1();
        ddpsi[1][1] = psi.eps1eps2();
        ///////////////////////////////////////////////
        a_HD[1].f1 = 0.0;
        a_HD[1].f2 = 0.0;
        a_HD[2].f1 = 1.0;
        a_HD[2].f2 = 1.0;
        a = sqrt(pow(a_HD[0], 2.0) + pow(a_HD[1], 2.0) + pow(a_HD[2], 2.0));
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        dpsi[2] = psi.eps1();
        ddpsi[2][2] = psi.eps1eps2();
        ////////////////////////////////////////////
        a_HD[0].f1 = 1.0;
        a_HD[1].f2 = 1.0;
        a_HD[2].f1 = 0.0;
        a_HD[2].f2 = 0.0;
        a = sqrt(pow(a_HD[0], 2.0) + pow(a_HD[1], 2.0) + pow(a_HD[2], 2.0));
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[0][1] = psi.eps1eps2();
        ddpsi[1][0] = ddpsi[0][1];
        ////////////////////////////////////////////
        a_HD[1].f2 = 0.0;
        a_HD[2].f2 = 1.0;
        a = sqrt(pow(a_HD[0], 2.0) + pow(a_HD[1], 2.0) + pow(a_HD[2], 2.0));
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[0][2] = psi.eps1eps2();
        ddpsi[2][0] = ddpsi[0][2];
        ///////////////////////////////////////////////
        a_HD[0].f1 = 0.0;
        a_HD[1].f1 = 1.0;
        a = sqrt(pow(a_HD[0], 2.0) + pow(a_HD[1], 2.0) + pow(a_HD[2], 2.0));
        s = log(a / A);
        psi = 0.5 * C2 * A * s * s;
        ddpsi[1][2] = psi.eps1eps2();
        ddpsi[2][1] = ddpsi[1][2];
        // CPP
        double dpsi_e[3][3] = {{0, dpsi[2], -dpsi[1]}, {-dpsi[2], 0, dpsi[0]}, {dpsi[1], -dpsi[0], 0}};
        double AA1[3][3], AA2[3][3];
        AA1[0][0] = -ddpsi[2][2] * psi_2_y * (2 * psi_1_y - 2 * psi_2_y);
        AA1[0][1] = +ddpsi[2][2] * psi_2_x * (2 * psi_1_y - 2 * psi_2_y);
        AA1[0][2] = -ddpsi[2][1] * psi_2_x * (2 * psi_1_y - 2 * psi_2_y);
        AA1[1][0] = +ddpsi[2][2] * psi_2_y * (2 * psi_1_x - 2 * psi_2_x);
        AA1[1][1] = -ddpsi[2][2] * psi_2_x * (2 * psi_1_x - 2 * psi_2_x);
        AA2[0][0] = +ddpsi[2][2] * psi_1_y * (2 * psi_1_y - 2 * psi_2_y);
        AA2[0][1] = -ddpsi[2][2] * psi_1_x * (2 * psi_1_y - 2 * psi_2_y);
        AA2[0][2] = +ddpsi[2][1] * psi_1_x * (2 * psi_1_y - 2 * psi_2_y);
        AA2[1][0] = -ddpsi[2][2] * psi_1_y * (2 * psi_1_x - 2 * psi_2_x);
        AA2[1][1] = +ddpsi[2][2] * psi_1_x * (2 * psi_1_x - 2 * psi_2_x);

        /////////////////////////////////////////////////////////////////

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                AA2I[i * 2 + j] = V2 * (-2.0 * dpsi_e[i][j] + AA1[i][j]);
                AA2J[i * 2 + j] = V2 * (2.0 * dpsi_e[i][j] + AA2[i][j]);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////                 Residual                    //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Compute_RK(const std::vector<Node> &Nodes, const std::vector<int> &node_idx, const double C1, const double C2, const double C3, const double horizon, const double V_base, const int max_num_neighb, std::vector<std::array<double, 2>> &R_arr, std::vector<std::vector<std::array<double, 4>>> &K_glob, double nn, bool verbose, int deb_idx, double &V1_total, double &V2_total)
{
    ///////////////////////////////////////
     #pragma omp parallel for
    for (int ele = 0; ele < node_idx.size(); ++ele) // iterate over compute nodes
    {

        int idx = node_idx[ele];
        int n1_size = Nodes[idx].neighbours1.size();           // read n1 neighbours
        std::array<double, 2> R_temp = {0.0, 0.0};             // R is neccessary for each Node and will overwrite for each n1 neighbour
        std::vector<std::array<double, 4>> K_loc(n1_size + 1); // K is neccessary for each Node
        double vol_frac = (n1_size + 1.0) / (max_num_neighb + 1);
        double V1 = vol_frac * V_base / n1_size;
        V1_total = V1_total + V1;
        hyperdual psi_1_HD[2];
        hyperdual l, s, psi;

        if (C1 != 0)
        {
            for (int n1 = 0; n1 < n1_size; ++n1)
            {
                // Ra
                double Psi_1_x = Nodes[Nodes[idx].neighbours1[n1]].X[0] - Nodes[idx].X[0];
                double Psi_1_y = Nodes[Nodes[idx].neighbours1[n1]].X[1] - Nodes[idx].X[1];
                double psi_1_x = Nodes[(Nodes[idx].neighbours1[n1])].coordinate[0] - Nodes[idx].coordinate[0];
                double psi_1_y = Nodes[(Nodes[idx].neighbours1[n1])].coordinate[1] - Nodes[idx].coordinate[1];
                psi_1_HD[0] = hyperdual(psi_1_x, 1.0, 0.0, 0.0);
                psi_1_HD[1] = hyperdual(psi_1_y, 0.0, 1.0, 0.0);
                ///////////////////////////////////////////////////////////////
                if (nn != 0)
                {
                    l = sqrt(psi_1_HD[0] * psi_1_HD[0] + psi_1_HD[1] * psi_1_HD[1]);
                    double L = std::sqrt(Psi_1_x * Psi_1_x + Psi_1_y * Psi_1_y);
                    s = (1.0 / nn) * (pow(l / L, nn) - 1);
                    //////////////////////////////////////////////
                    double C11 = 4 * C1 * abs(1 - (L/horizon,2));
                    //////////////////////////////////////////////
                    psi = 0.5 * C11 * L * s * s;
                    R_temp[0] += psi.eps1() * V1;
                    R_temp[1] += psi.eps2() * V1;
                    //// HDN-R ////////////////////
                    K_loc[n1 + 1][1] = V1 * psi.eps1eps2();
                    K_loc[n1 + 1][2] = K_loc[n1 + 1][1];
                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    if (idx == deb_idx && verbose)
                    {
                        std::cout << " L is = " << L << std::endl;
                        std::cout << "for neighbor1 with index = " << Nodes[idx].neighbours1[n1] << " and compute index = " << Nodes[Nodes[idx].neighbours1[n1]].global_idx << " R_temp_x = " << psi.eps1() * V1 << " and R_temp_y = " << psi.eps2() * V1 << std::endl;
                    }
                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    psi_1_HD[0].f2 = 1;
                    psi_1_HD[1].f2 = 0;
                    l = sqrt(psi_1_HD[0] * psi_1_HD[0] + psi_1_HD[1] * psi_1_HD[1]);
                    s = (1.0 / nn) * (pow(l / L, nn) - 1);
                    psi = 0.5 * C11 * L * s * s;
                    K_loc[n1 + 1][0] = V1 * psi.eps1eps2();
                    psi_1_HD[0].f1 = 0.0;
                    psi_1_HD[0].f2 = 0.0;
                    psi_1_HD[1].f1 = 1.0;
                    psi_1_HD[1].f2 = 1.0;
                    l = sqrt(psi_1_HD[0] * psi_1_HD[0] + psi_1_HD[1] * psi_1_HD[1]);
                    s = (1.0 / nn) * (pow(l / L, nn) - 1);
                    psi = 0.5 * C11 * L * s * s;
                    K_loc[n1 + 1][3] = V1 * psi.eps1eps2();

                    //     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    if (idx == deb_idx && verbose)
                    {
                        std::cout << "ps1_x = " << psi_1_x << " , ps1_Y =  " << psi_1_y << std::endl;
                        std::cout << "K1 of neighbor with n1 in neighbor1 list = " << n1 << " is  = " << K_loc[n1 + 1][0] << " and " << K_loc[n1 + 1][1] << " and " << K_loc[n1 + 1][2] << " and " << K_loc[n1 + 1][3] << std::endl;
                        std::cout << " ============================ " << std::endl;
                    }
                    //     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                    K_loc[0][0] -= K_loc[n1 + 1][0];
                    K_loc[0][1] -= K_loc[n1 + 1][1];
                    K_loc[0][2] -= K_loc[n1 + 1][2];
                    K_loc[0][3] -= K_loc[n1 + 1][3];
                }
                else
                {
                    double pp, alpha;
                    hyperdual s_1, s_2;
                    double ddpsi[2][2];
                    l = sqrt(psi_1_HD[0] * psi_1_HD[0] + psi_1_HD[1] * psi_1_HD[1]);
                    double L = std::sqrt(Psi_1_x * Psi_1_x + Psi_1_y * Psi_1_y);
                    pp = 2.0;
                    alpha = 0.5;
                    s = log(l / L);
                    //////////////////////////////////////////////
                    double C11 = 4 * C1 * abs(1 - (L/horizon,2));
                    //////////////////////////////////////////////
                    psi = 0.5 * C11 * L * s * s;
                    R_temp[0] += psi.eps1() * V1;
                    R_temp[1] += psi.eps2() * V1;
                    /////////////////////////////////////////////////////////////////
                    ddpsi[1][0] = psi.eps1eps2();
                    ddpsi[0][1] = ddpsi[1][0];
                    psi_1_HD[0].f1 = 0.0;
                    psi_1_HD[0].f2 = 0.0;
                    psi_1_HD[1].f1 = 1.0;
                    psi_1_HD[1].f2 = 1.0;
                    l = sqrt(pow(psi_1_HD[0], 2.0) + pow(psi_1_HD[1], 2.0));
                    s = log(l / L);
                    psi = 0.5 * C11 * L * s * s;
                    ddpsi[1][1] = psi.eps1eps2();
                    psi_1_HD[1].f1 = 0.0;
                    psi_1_HD[1].f2 = 0.0;
                    psi_1_HD[0].f1 = 1.0;
                    psi_1_HD[0].f2 = 1.0;
                    l = sqrt(pow(psi_1_HD[0], 2.0) + pow(psi_1_HD[1], 2.0));
                    s = log(l / L);
                    psi = 0.5 * C11 * L * s * s;
                    ddpsi[0][0] = psi.eps1eps2();

                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            K_loc[n1 + 1][i * 2 + j] = V1 * ddpsi[j][i];
                        }
                    }
                    if (idx == deb_idx && verbose)
                    {
                        std::cout << "ps1_x = " << psi_1_x << " , ps1_Y =  " << psi_1_y << std::endl;
                        std::cout << "K1 of neighbor with n1 in neighbor1 list = " << n1 << " is  = " << K_loc[n1 + 1][0] << " and " << K_loc[n1 + 1][1] << " and " << K_loc[n1 + 1][2] << " and " << K_loc[n1 + 1][3] << std::endl;
                        std::cout << " ============================ " << std::endl;
                    }
                    K_loc[0][0] -= K_loc[n1 + 1][0];
                    K_loc[0][1] -= K_loc[n1 + 1][1];
                    K_loc[0][2] -= K_loc[n1 + 1][2];
                    K_loc[0][3] -= K_loc[n1 + 1][3];
                    ////////////////////////////////////////////////////////////////
                }
            }
            if (idx == deb_idx && verbose)
            {

                std::cout << "K1_total_aa   is  x= " << K_loc[0][0] << " and " << K_loc[0][1] << " and " << K_loc[0][2] << " and " << K_loc[0][3] << std::endl;
            }

            K_glob[ele] = K_loc;
            R_arr[ele] = R_temp;

            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if (idx == deb_idx && verbose)
            {
                std::cout << "R1 final is  x = " << R_arr[ele][0] << "R1 final  is  y = " << R_arr[ele][1] << std::endl;
            }
        }
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (C2 != 0)
        {
            double V2 = V_base * V1 * V1 / (V_base - V1);
            V2_total = V2 + V2_total;
            R_temp = {0.0, 0.0};
            hyperdual a_HD[3], a;
            double A;
            double dpsi[3];
            double ddpsi[3][3];

            for (int n2 = 0; n2 < Nodes[idx].neighbours2.size(); ++n2) // iterating over pairs
            {

                auto n1_idx_it = std::find(Nodes[idx].neighbours1.begin(), Nodes[idx].neighbours1.end(), Nodes[idx].neighbours2[n2][0]);
                auto n2_idx_it = std::find(Nodes[idx].neighbours1.begin(), Nodes[idx].neighbours1.end(), Nodes[idx].neighbours2[n2][1]);
                int n1_idx = std::distance(Nodes[idx].neighbours1.begin(), n1_idx_it);
                int n2_idx = std::distance(Nodes[idx].neighbours1.begin(), n2_idx_it);

                if (idx == deb_idx && verbose)
                {
                    std::cout << "n1_idx = " << n1_idx << "n2_idx = " << n2_idx << std::endl;
                }

                double Psi_1_x = Nodes[Nodes[idx].neighbours1[n1_idx]].X[0] - Nodes[idx].X[0];
                double Psi_1_y = Nodes[Nodes[idx].neighbours1[n1_idx]].X[1] - Nodes[idx].X[1];
                double Psi_1_z = 0;
                double psi_1_x = Nodes[Nodes[idx].neighbours1[n1_idx]].coordinate[0] - Nodes[idx].coordinate[0];
                double psi_1_y = Nodes[Nodes[idx].neighbours1[n1_idx]].coordinate[1] - Nodes[idx].coordinate[1];
                double psi_1_z = 0;
                double Psi_2_x = Nodes[Nodes[idx].neighbours1[n2_idx]].X[0] - Nodes[idx].X[0];
                double Psi_2_y = Nodes[Nodes[idx].neighbours1[n2_idx]].X[1] - Nodes[idx].X[1];
                double Psi_2_z = 0;
                double psi_2_x = Nodes[Nodes[idx].neighbours1[n2_idx]].coordinate[0] - Nodes[idx].coordinate[0];
                double psi_2_y = Nodes[Nodes[idx].neighbours1[n2_idx]].coordinate[1] - Nodes[idx].coordinate[1];
                double psi_2_z = 0;
                double AA2I[4];
                double AA2J[4];

                ////////////////////////////////////////////////////////////////////////
                if (nn != 0)
                {
                    double A_tmp[3] = {0.0, 0.0, Psi_1_x * Psi_2_y - Psi_1_y * Psi_2_x};
                    double a_tmp[3] = {0.0, 0.0, psi_1_x * psi_2_y - psi_1_y * psi_2_x};
                    a_HD[0] = hyperdual(a_tmp[0], 1.0, 0.0, 0.0);
                    a_HD[1] = hyperdual(a_tmp[1], 0.0, 1.0, 0.0);
                    a_HD[2] = hyperdual(a_tmp[2], 0.0, 0.0, 0.0);
                    A = fabs(A_tmp[2]);
                    a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
                    s = (1.0 / nn) * (pow(a / A, nn) - 1);
                    psi = 0.5 * C2 * A * s * s;
                    dpsi[0] = psi.eps1();
                    dpsi[1] = psi.eps2();
                    /////////////////////////////////////////////////////////
                    a_HD[0].f1 = 0.0;
                    a_HD[1].f2 = 0.0;
                    a_HD[2].f1 = 1.0;

                    a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
                    s = (1.0 / nn) * (pow(a / A, nn) - 1);
                    psi = 0.5 * C2 * A * s * s;
                    dpsi[2] = psi.eps1();
                    ////////////////////////////////////////////////////////////
                    R_temp[0] = V2 * 2.0 * ((psi_2_y * dpsi[2] - psi_2_z * dpsi[1]) - (psi_1_y * dpsi[2] - psi_1_z * dpsi[1]));
                    R_temp[1] = V2 * 2.0 * ((psi_2_z * dpsi[0] - psi_2_x * dpsi[2]) - (psi_1_z * dpsi[0] - psi_1_x * dpsi[2]));

                    if (idx == deb_idx && verbose)
                    {
                        std::cout << " =  = = = = = == = = " << std::endl;
                        std::cout << "for neighbor1 with index = " << Nodes[idx].neighbours2[n2][0] << " and compute index = " << Nodes[Nodes[idx].neighbours2[n2][0]].global_idx << " R2_temp_x = " << R_temp[0] << " and R_temp_y = " << R_temp[1] << std::endl;
                    }
                    //////////////////////////////////////////////////////////////////////////////////////////////
                    // K2
                    if (idx == deb_idx && verbose)
                    {

                        std::cout << "N2_ idx1 = " << n1_idx << " idx2 = " << n2_idx << std::endl;
                        std::cout << "psi1 = " << psi_1_x << " psi1 = " << psi_1_y << " psi1 = " << psi_1_z << std::endl;
                        std::cout << "Psi1 = " << Psi_1_x << " Psi1 = " << Psi_1_y << " Psi1 = " << Psi_1_z << std::endl;

                        std ::cout << "psi2 = " << psi_2_x << " psi2 = " << psi_2_y << " psi2 = " << psi_2_z << std::endl;

                        std ::cout << "Psi2 = " << Psi_2_x << " Psi2 = " << Psi_2_y << " Psi2 = " << Psi_2_z << std::endl;
                    }
                    Cal_K2(nn, idx, deb_idx, C2, V2, Psi_1_x, Psi_1_y, Psi_1_z, Psi_2_x, Psi_2_y, Psi_2_z, psi_1_x, psi_1_y, psi_1_z, psi_2_x, psi_2_y, psi_2_z, AA2I, AA2J);

                    K_loc[n1_idx + 1][0] += AA2I[0];
                    K_loc[n1_idx + 1][1] += AA2I[1];
                    K_loc[n1_idx + 1][2] += AA2I[2];
                    K_loc[n1_idx + 1][3] += AA2I[3];
                    ///
                    K_loc[n2_idx + 1][0] += AA2J[0];
                    K_loc[n2_idx + 1][1] += AA2J[1];
                    K_loc[n2_idx + 1][2] += AA2J[2];
                    K_loc[n2_idx + 1][3] += AA2J[3];
                    ///
                    K_loc[0][0] -= AA2I[0] + AA2J[0];
                    K_loc[0][1] -= AA2I[1] + AA2J[1];
                    K_loc[0][2] -= AA2I[2] + AA2J[2];
                    K_loc[0][3] -= AA2I[3] + AA2J[3];
                    // R_arr[ele][0] += R_temp[0];
                    // R_arr[ele][1] += R_temp[1];
                }
                else
                {
                    double A_tmp[3] = {0.0, 0.0, Psi_1_x * Psi_2_y - Psi_1_y * Psi_2_x};
                    double a_tmp[3] = {0.0, 0.0, psi_1_x * psi_2_y - psi_1_y * psi_2_x};
                    a_HD[0] = hyperdual(a_tmp[0], 1.0, 0.0, 0.0);
                    a_HD[1] = hyperdual(a_tmp[1], 0.0, 1.0, 0.0);
                    a_HD[2] = hyperdual(a_tmp[2], 0.0, 0.0, 0.0);
                    A = fabs(A_tmp[2]);
                    a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
                    s = log(a / A);
                    psi = 0.5 * C2 * A * s * s;
                    dpsi[0] = psi.eps1();
                    dpsi[1] = psi.eps2();
                    a_HD[0].f1 = 0.0;
                    a_HD[1].f2 = 0.0;
                    a_HD[2].f1 = 1.0;
                    a = sqrt(a_HD[0] * a_HD[0] + a_HD[1] * a_HD[1] + a_HD[2] * a_HD[2]);
                    s = log(a / A);
                    psi = 0.5 * C2 * A * s * s;
                    dpsi[2] = psi.eps1();
                    R_temp[0] = V2 * 2.0 * ((psi_2_y * dpsi[2]) - (psi_1_y * dpsi[2]));
                    R_temp[1] = V2 * 2.0 * ((-psi_2_x * dpsi[2]) + (psi_1_x * dpsi[2]));

                    Cal_K2(nn, idx, deb_idx, C2, V2, Psi_1_x, Psi_1_y, Psi_1_z, Psi_2_x, Psi_2_y, Psi_2_z, psi_1_x, psi_1_y, psi_1_z, psi_2_x, psi_2_y, psi_2_z, AA2I, AA2J);
                    K_loc[n1_idx + 1][0] += AA2I[0];
                    K_loc[n1_idx + 1][1] += AA2I[1];
                    K_loc[n1_idx + 1][2] += AA2I[2];
                    K_loc[n1_idx + 1][3] += AA2I[3];
                    K_loc[n2_idx + 1][0] += AA2J[0];
                    K_loc[n2_idx + 1][1] += AA2J[1];
                    K_loc[n2_idx + 1][2] += AA2J[2];
                    K_loc[n2_idx + 1][3] += AA2J[3];
                    K_loc[0][0] -= AA2J[0] + AA2I[0];
                    K_loc[0][1] -= AA2J[1] + AA2I[1];
                    K_loc[0][2] -= AA2J[2] + AA2I[2];
                    K_loc[0][3] -= AA2J[3] + AA2I[3];
                }
                R_arr[ele][0] += R_temp[0];
                R_arr[ele][1] += R_temp[1];
            }
        }
        K_glob[ele] = K_loc;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////                 Save and write R                    /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void save_R(const std::vector<std::array<double, 2>> &R, const std::vector<int> &computeNodesIdx, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < computeNodesIdx.size(); ++i)
        {
            outputFile << " Node: " << i << " Global_idx = " << computeNodesIdx[i] << " R = " << R[i][0] << " , " << R[i][1] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void Write_R(const std::vector<std::array<double, 2>> &R, const std::vector<int> &computeNodesIdx, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < computeNodesIdx.size(); ++i)
        {
            outputFile << R[i][0] << std::endl;
            outputFile << R[i][1] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////                 Save and write K                    /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Write_K(const std::vector<std::vector<std::array<double, 4>>> &K_glob, const std::vector<int> &computeNodesIdx, std::vector<Node> &Nodes, const std::string &filename)
{
    int size = computeNodesIdx.size();
    std::vector<std::array<double, 4>> full_stiffness(size * size, std::array<double, 4>{0.0, 0.0, 0.0, 0.0});
    for (size_t i = 0; i < size; ++i)
    {
        int global_indx = computeNodesIdx[i];

        for (size_t j = 0; j < Nodes[global_indx].neighbours1.size(); ++j)
        {
            int n1_global_idx = Nodes[Nodes[global_indx].neighbours1[j]].global_idx;

            if (n1_global_idx == -1)
            {
                continue;
            }
            for (size_t k = 0; k < 4; ++k)
            {

                full_stiffness[i * size + n1_global_idx][k] = K_glob[i][j + 1][k];
            }
        }
        for (size_t k = 0; k < 4; ++k)
        {
            full_stiffness[i * size + i][k] = K_glob[i][0][k];
        }
    }
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < size; ++j)
            {

                outputFile << full_stiffness[i * size + j][0] << " , " << full_stiffness[i * size + j][1] << " , ";
            }
            outputFile << std::endl;

            for (size_t j = 0; j < size; ++j)
            {

                outputFile << full_stiffness[i * size + j][2] << " , " << full_stiffness[i * size + j][3] << " , ";
            }
            // std::cout << std::endl;

            outputFile << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void save_K(const std::vector<std::vector<std::array<double, 4>>> &K_glob, const std::vector<int> &computeNodesIdx, std::vector<Node> &Nodes, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < computeNodesIdx.size(); ++i)
        {
            // std::cout << "compute node = " << i << std::endl;
            int global_indx = computeNodesIdx[i];
            // std::cout << global_indx << std::endl;
            outputFile << " Node: " << i << " Global_idx = " << global_indx << " ( " << i << ", " << i << ") = " << K_glob[i][0][0] << " , " << K_glob[i][0][1] << " , " << K_glob[i][0][2] << " , " << K_glob[i][0][3];
            // std::cout << " passed" << std::endl;
            // std::cout << "size of n1 = " << Nodes[global_indx].neighbours1.size() << std::endl;
            for (int j = 0; j < Nodes[global_indx].neighbours1.size(); ++j)
            {
                int neigh_idx = Nodes[global_indx].neighbours1[j];
                // std::cout << global_indx << " and " << j << "and " << neigh_idx << std::endl;
                outputFile << " K (" << i << ", " << Nodes[neigh_idx].global_idx << " ) = " << K_glob[i][j + 1][0] << " , " << K_glob[i][j + 1][1] << " , " << K_glob[i][j + 1][2] << " , " << K_glob[i][j + 1][3];
            }
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////                 Counting N1 neighbors                    ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void count_n1(std::vector<Node> &Nodes)
{

    for (int j = 0; j < Nodes.size(); ++j)
    {
        std::cout << "size of neighbors1 = " << Nodes[j].neighbours1.size() << std::endl;
        for (int i = 0; i < Nodes[j].neighbours1.size(); ++i)
        {
            std::cout << "neighbor compute indx is = " << Nodes[Nodes[j].neighbours1[i]].global_idx << std::endl;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Write_dx(double *dx, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << dx[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_ia(int *ia, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << ia[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_ja(int *ja, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << ja[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_data(double *k, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << k[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_ai(const std::unique_ptr<int[]> &ai, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << ai[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_aj(const std::unique_ptr<int[]> &aj, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << aj[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_R(double *R, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << R[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}


double cal_norm2_R(const EMx &R_arr, int size)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+ : sum2)
    for (int i = 0; i < size; ++i)
    {
        sum2 += R_arr(i) * R_arr(i);
    }
    return std::sqrt(sum2);
}

void Write_ai(std::vector<int> ai, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << ai[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
void Write_aj(std::vector<int> aj, int size, const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        for (size_t i = 0; i < size; ++i)
        {
            outputFile << aj[i] << std::endl;
        }
        outputFile.close();
        std::cout << "Output saved to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}