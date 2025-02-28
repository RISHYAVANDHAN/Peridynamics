//
// Created by srini on 22/11/2024.
//

#ifndef POINTS_H
#define POINTS_H
#include <iostream>
#include <vector>
#include <utility> //std::pair - just fyi
#include <tuple>
#include <Eigen/Dense>




class Points {
public:
    int Nr{};                                               // Nr - position number or number of points
    Eigen::VectorXd X;                                      // X material coordinates
    Eigen::VectorXd x;                                      // x spatial coordinates
    double volume{};                                        // Volume
    std::vector<std::vector<int>> point_list;
    std::vector<std::vector<int>> neighbour_list_1N;        // consolidated neighbour_list containing 1-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_2N;        // consolidated neighbour_list containing 2-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_3N;        // consolidated neighbour_list containing 3-neighbour interaction
    Eigen::VectorXd Ra_1;                                   // Residual of 1 - neighbour interaction
    Eigen::VectorXd Ra_2;                                   // Residual of 2 - neighbour interaction
    Eigen::VectorXd Ra_3;                                   // Residual of 3 - neighbour interaction
    Eigen::VectorXd Ra_sum;                                 // Sum of residuals
    Eigen::MatrixXd Kab_1;                                  // Tangential Stiffness of 1 - neighbour interaction
    Eigen::MatrixXd Kab_2;                                  // Tangential Stiffness of 2 - neighbour interaction
    Eigen::MatrixXd Kab_3;                                  // Tangential Stiffness of 3 - neighbour interaction
    std::string Flag = "Patch";                             // Flag (Patches & points)
    //int DOCs = 0;                                         // temporary names DOC and DOF for global free index and global fixed index
    Eigen::VectorXd DOF;
    std::vector<double> DOC;
    Eigen::Vector3d BC;                                     // 0 - patches, 1 - points
    Eigen::Vector3d BCval;                                  // displacement and force calculation
    size_t n1 = 0;                                          // n1 - no. of 1 neighbour interaction
    size_t n2 = 0;                                          // n2 - no. of 2 neighbour interaction
    size_t n3 = 0;                                          // n3 - no. of 3 neighbour interaction
    double V_eff;                                           // this is used in assembly for some reason, that´s the only reason
    double psi;                                             // same reason

    // Constructor for initializing Points object with necessary data
    Points(int Nr, const Eigen::VectorXd& X, const Eigen::VectorXd& x, double volume)
    : Nr(Nr),
      X(Eigen::VectorXd(3, 0.0)),
      x(Eigen::VectorXd(3, 0.0)),
      volume(volume)
    {
        if (X.size() != 3 || x.size() != 3) {
            std::cerr << "Warning: Input vectors must be size 3. Initializing with zeros.\n";
        }
    }

    Points() = default;

};



#endif //POINTS_H

