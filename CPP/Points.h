//
// Created by srini on 22/11/2024.
//

#ifndef POINTS_H
#define POINTS_H
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class Points {
public:
    int Nr{};                                               // Nr - position number or number of points
    Eigen::Vector3d X;                                      // X material coordinates
    Eigen::Vector3d x;                                      // x spatial coordinates
    double volume = 1.0;                                        // Volume
    std::vector<std::vector<int>> neighbour_list_1N;        // consolidated neighbour_list containing 1-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_2N;        // consolidated neighbour_list containing 2-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_3N;        // consolidated neighbour_list containing 3-neighbour interaction
    std::vector<std::vector<int>> neighbour_list;           // consolidated neighbour_list containing all possible interactions (1-, 2- and 3- neighbour)
    Eigen::Vector3d Ra_1;                                   // Residual of 1 - neighbour interaction
    Eigen::Vector3d Ra_2;                                   // Residual of 2 - neighbour interaction
    Eigen::Vector3d Ra_3;                                   // Residual of 3 - neighbour interaction
    Eigen::Vector3d Ra_sum;                                 // Sum of residuals
    Eigen::Matrix3d Kab_1;                                  // Tangential Stiffness of 1 - neighbour interaction
    Eigen::Matrix3d Kab_2;                                  // Tangential Stiffness of 2 - neighbour interaction
    Eigen::Matrix3d Kab_3;                                  // Tangential Stiffness of 3 - neighbour interaction
    Eigen::Matrix3d Kab_sum;                                // Sum of tangential stiffness matrices
    std::string Flag = "Patch";                             // Flag (Patches & points)
    Eigen::Vector3d DOF;                                    // Degrees of freedom
    std::vector<double> DOC;                                // Degrees of constraint
    Eigen::Vector3d BC;                                     // Boundary condition
    Eigen::Vector3d BCval;                                  // Boundary condition value
    size_t n1 = 0;                                          // n1 - no. of 1 neighbour interaction
    size_t n2 = 0;                                          // n2 - no. of 2 neighbour interaction
    size_t n3 = 0;                                          // n3 - no. of 3 neighbour interaction
    double V_eff;                                           // Effective volume
    double psi;                                             // Energy density

    // Constructor for initializing Points object with necessary data
    Points(int Nr, const Eigen::Vector3d& X, const Eigen::Vector3d& x, double volume)
        : Nr(Nr), X(X), x(x), volume(volume) {}

    Points() = default;
};


#endif //POINTS_H

