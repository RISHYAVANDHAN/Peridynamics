#ifndef SOLVER_H
#define SOLVER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include "Points.h"

inline double compute_psi(double C, double L, double l) {
    double s = l / L;
    return 0.5 * C * L * std::pow((s - 1), 2);
}

inline Eigen::Matrix3d compute_Kab_1(const Eigen::Vector3d& xi_1, double L, double l, double V_eff, double C1) {
    Eigen::Matrix3d Kab;
    Eigen::Matrix3d xi_outer = xi_1 * xi_1.transpose();
    double l_cubed = std::pow(l, 3);

    Eigen::Matrix3d term1 = (1.0 / l_cubed) * xi_outer;
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    double scalar_term = ((1.0 / L) - (1.0 / l));

    Kab = C1 * V_eff * (term1 + scalar_term * identity);
    return Kab;
}

inline Eigen::Matrix3d compute_Kab_2(const Eigen::Vector3d& xi_1, const Eigen::Vector3d& xi_2, const Eigen::Vector3d& Xi_1, const Eigen::Vector3d& Xi_2, double A, double a, double V_eff, double C2) {
    Eigen::Matrix3d Kab = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    double xi1_dot_xi2 = xi_1.dot(xi_2);
    double xi2_dot_xi2 = xi_2.dot(xi_2);
    double xi1_dot_xi1 = xi_1.dot(xi_1);

    Eigen::Vector3d term_vec1 = xi1_dot_xi2 * xi_2 - xi2_dot_xi2 * xi_1;
    Eigen::Vector3d term_vec2 = xi1_dot_xi2 * xi_1 - xi1_dot_xi1 * xi_2;

    Eigen::Matrix3d term1 = (2.0 * C2 * V_eff / std::pow(a, 3)) * (term_vec1 * term_vec1.transpose());
    double scalar_term = (1.0 / a) - (1.0 / A);
    Eigen::Matrix3d term2 = (2.0 * C2 * V_eff * scalar_term) * (xi_2 * xi_2.transpose() - xi1_dot_xi2 * I);
    Eigen::Matrix3d term3 = (2.0 * C2 * V_eff / std::pow(a, 3)) * (term_vec1 * term_vec2.transpose());
    Eigen::Matrix3d term4 = (2.0 * C2 * V_eff * scalar_term) * (xi_2 * xi_1.transpose() + xi1_dot_xi2 * I - 2.0 * xi_1 * xi_2.transpose());

    Kab = term1 + term2 + term3 + term4;
    return Kab;
}

inline Eigen::Matrix3d compute_Kab_3(const Eigen::Vector3d& xi_1, const Eigen::Vector3d& xi_2, const Eigen::Vector3d& xi_3, const Eigen::Vector3d& Xi_1, const Eigen::Vector3d& Xi_2, const Eigen::Vector3d& Xi_3, double V, double v, double V_eff, double C3) {
    Eigen::Matrix3d Kab = Eigen::Matrix3d::Zero();

    Eigen::Vector3d xi_12 = xi_1.cross(xi_2);
    Eigen::Vector3d xi_23 = xi_2.cross(xi_3);
    Eigen::Vector3d xi_31 = xi_3.cross(xi_1);

    double xi_123 = xi_1.dot(xi_2.cross(xi_3));
    double Xi_123 = Xi_1.dot(Xi_2.cross(Xi_3));

    double scalar_term = 1.0 / xi_123;
    double scalar_term_ref = 1.0 / Xi_123;
    double diff_term = scalar_term - scalar_term_ref;

    Eigen::Matrix3d term1 = (3.0 * C3 * V_eff * scalar_term) * (xi_23 * xi_23.transpose());
    Eigen::Matrix3d epsilon_matrix;
    epsilon_matrix << 0, -xi_3[2], xi_3[1], xi_3[2], 0, -xi_3[0], -xi_3[1], xi_3[0], 0;
    Eigen::Matrix3d term2 = (3.0 * C3 * V_eff * diff_term) * ((xi_12.dot(xi_3)) * epsilon_matrix);
    Eigen::Matrix3d term3 = (3.0 * C3 * V_eff * scalar_term) * (xi_23 * xi_31.transpose());
    Eigen::Matrix3d epsilon_matrix_2;
    epsilon_matrix_2 << 0, -xi_2[2], xi_2[1], xi_2[2], 0, -xi_2[0], -xi_2[1], xi_2[0], 0;
    Eigen::Matrix3d term4 = -(3.0 * C3 * V_eff * diff_term) * ((xi_12.dot(xi_3)) * epsilon_matrix_2);
    Eigen::Matrix3d term5 = (3.0 * C3 * V_eff * scalar_term) * (xi_23 * xi_12.transpose());

    Kab = term1 + term2 + term3 + term4 + term5;
    return Kab;
}

inline void calculate_r(int PD, std::vector<Points>& point_list, double C1, double C2, double C3, double delta) {
    const double pi = 3.14159265358979323846;
    double Vh = (4.0 / 3.0) * pi * (delta * delta * delta);

    for (size_t i = 0; i < point_list.size(); i++) {
        point_list[i].psi = 0.0;
        point_list[i].Ra_1 = Eigen::Vector3d::Zero();
        point_list[i].Ra_2 = Eigen::Vector3d::Zero();
        point_list[i].Ra_3 = Eigen::Vector3d::Zero();
        point_list[i].Kab_1 = Eigen::Matrix3d::Zero();
        point_list[i].Kab_2 = Eigen::Matrix3d::Zero();
        point_list[i].Kab_3 = Eigen::Matrix3d::Zero();

        // 1-Neighbors
        for (size_t j = 0; j < point_list[i].neighbour_list_1N.size(); j++) {
            int neighbor_idx = point_list[i].neighbour_list_1N[j][0];

            Eigen::Vector3d Xi_1 = point_list[i].X - point_list[neighbor_idx].X;
            Eigen::Vector3d xi_1 = point_list[i].x - point_list[neighbor_idx].x;

            double L = Xi_1.norm();
            double l = xi_1.norm();

            point_list[i].psi += compute_psi(C1, L, l);
            point_list[i].V_eff = Vh / point_list[i].n1;

            point_list[i].Ra_1 += xi_1 * C1 * point_list[i].V_eff * ((1 / L) - (1 / l));
            point_list[i].Kab_1 += compute_Kab_1(xi_1, L, l, point_list[i].V_eff, C1);
        }

        // 2-Neighbors
        if (PD >= 2) {
            for (size_t j = 0; j < point_list[i].neighbour_list_2N.size(); j++) {
                int neighbor1_idx = point_list[i].neighbour_list_2N[j][0];
                int neighbor2_idx = point_list[i].neighbour_list_2N[j][1];

                Eigen::Vector3d Xi_1 = point_list[i].X - point_list[neighbor1_idx].X;
                Eigen::Vector3d xi_1 = point_list[i].x - point_list[neighbor1_idx].x;
                Eigen::Vector3d Xi_2 = point_list[i].X - point_list[neighbor2_idx].X;
                Eigen::Vector3d xi_2 = point_list[i].x - point_list[neighbor2_idx].x;

                Eigen::Vector3d Xi = Xi_1.cross(Xi_2);
                Eigen::Vector3d xi = xi_1.cross(xi_2);

                double A = Xi.norm();
                double a = xi.norm();

                point_list[i].psi += compute_psi(C2, A, a);
                point_list[i].V_eff = (Vh * Vh) / (point_list[i].n1 + point_list[i].n2);

                Eigen::Vector3d calc = (xi_2.dot(xi_2) * xi_1 - (xi_1.dot(xi_2) * xi_2));
                point_list[i].Ra_2 += 2 * C2 * ((1 / A) - (1 / a)) * point_list[i].V_eff * calc;
                point_list[i].Kab_2 += compute_Kab_2(xi_1, xi_2, Xi_1, Xi_2, A, a, point_list[i].V_eff, C2);
            }
        }

        // 3-Neighbors
        if (PD == 3) {
            for (size_t j = 0; j < point_list[i].neighbour_list_3N.size(); j++) {
                int neighbor1_idx = point_list[i].neighbour_list_3N[j][0];
                int neighbor2_idx = point_list[i].neighbour_list_3N[j][1];
                int neighbor3_idx = point_list[i].neighbour_list_3N[j][2];

                Eigen::Vector3d Xi_1 = point_list[i].X - point_list[neighbor1_idx].X;
                Eigen::Vector3d xi_1 = point_list[i].x - point_list[neighbor1_idx].x;
                Eigen::Vector3d Xi_2 = point_list[i].X - point_list[neighbor2_idx].X;
                Eigen::Vector3d xi_2 = point_list[i].x - point_list[neighbor2_idx].x;
                Eigen::Vector3d Xi_3 = point_list[i].X - point_list[neighbor3_idx].X;
                Eigen::Vector3d xi_3 = point_list[i].x - point_list[neighbor3_idx].x;

                double Xi = Xi_1.dot(Xi_2.cross(Xi_3));
                double xi = xi_1.dot(xi_2.cross(xi_3));

                double V = Xi;
                double v = xi;

                point_list[i].psi += compute_psi(C3, V, v);
                point_list[i].V_eff = (Vh * Vh * Vh) / (point_list[i].n1 + point_list[i].n2 + point_list[i].n3);

                double calc = xi_2.dot(xi_1.cross(xi_2));
                Eigen::Vector3d calc1 = xi_2.cross(xi_3);

                point_list[i].Ra_3 += 3 * C3 * ((1 / V) - (1 / v)) * point_list[i].V_eff * calc * calc1;
                point_list[i].Kab_3 += compute_Kab_3(xi_1, xi_2, xi_3, Xi_1, Xi_2, Xi_3, V, v, point_list[i].V_eff, C3);
            }
        }

        point_list[i].Ra_sum = point_list[i].Ra_1 + point_list[i].Ra_2 + point_list[i].Ra_3;
        point_list[i].Kab_sum = point_list[i].Kab_1 + point_list[i].Kab_2 + point_list[i].Kab_3;
    }
}

#endif // SOLVER_H