#ifndef POINTS_H
#define POINTS_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class Points {
public:
    int Nr{};                                               // Nr - position number or number of points
    Eigen::Vector3d X = Eigen::Vector3d::Zero();            // X material coordinates
    Eigen::Vector3d x = Eigen::Vector3d::Zero();            // x spatial coordinates
    Eigen::Vector3d x_prev = Eigen::Vector3d::Zero();       // previous step x spatial coordinates
    double volume{};                                        // Volume
    std::vector<std::vector<int>> neighbour_list_1N;        // consolidated neighbour_list containing 1-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_2N;        // consolidated neighbour_list containing 2-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_3N;        // consolidated neighbour_list containing 3-neighbour interaction
    std::vector<std::vector<int>> neighbour_list;           // consolidated neighbour_list containing all possible interactions (1-, 2- and 3- neighbour)
    Eigen::Vector3d Ra_1 = Eigen::Vector3d::Zero();         // Residual of 1 - neighbour interaction
    Eigen::Vector3d Ra_2 = Eigen::Vector3d::Zero();         // Residual of 2 - neighbour interaction
    Eigen::Vector3d Ra_3 = Eigen::Vector3d::Zero();         // Residual of 3 - neighbour interaction
    Eigen::Vector3d Ra_sum = Eigen::Vector3d::Zero();       // Sum of residuals
    Eigen::Matrix3d Kab_1 = Eigen::Matrix3d::Zero();        // Tangential Stiffness of 1 - neighbour interaction
    Eigen::Matrix3d Kab_2 = Eigen::Matrix3d::Zero();        // Tangential Stiffness of 2 - neighbour interaction
    Eigen::Matrix3d Kab_3 = Eigen::Matrix3d::Zero();        // Tangential Stiffness of 3 - neighbour interaction
    Eigen::Matrix3d Kab_sum = Eigen::Matrix3d::Zero();      // Sum of tangential stiffness matrices
    std::string Flag = "Patch";                             // Flag (Patches & points)
    Eigen::Vector3i DOF = Eigen::Vector3i::Zero();          // Degrees of freedom
    std::vector<double> DOC;                                // Degrees of constraint
    Eigen::Vector3i BC = Eigen::Vector3i::Zero();           // Boundary condition
    Eigen::Vector3d BCval = Eigen::Vector3d::Zero();        // Boundary condition value
    size_t n1 = 0;                                          // n1 - no. of 1 neighbour interaction
    size_t n2 = 0;                                          // n2 - no. of 2 neighbour interaction
    size_t n3 = 0;                                          // n3 - no. of 3 neighbour interaction
    double V_eff = 0.0;                                     // Effective volume
    double psi = 0.0;                                       // Energy density

    // Constructor for initializing Points object with necessary data
    Points(int Nr, const Eigen::Vector3d& X, const Eigen::Vector3d& x, double volume)
        : Nr(Nr), X(X), x(x), volume(volume) {
        // Initialize all member variables

        Ra_1 = Eigen::Vector3d::Zero();
        Ra_2 = Eigen::Vector3d::Zero();
        Ra_3 = Eigen::Vector3d::Zero();
        Ra_sum = Eigen::Vector3d::Zero();
        Kab_1 = Eigen::Matrix3d::Zero();
        Kab_2 = Eigen::Matrix3d::Zero();
        Kab_3 = Eigen::Matrix3d::Zero();
        Kab_sum = Eigen::Matrix3d::Zero();
        DOF = Eigen::Vector3i::Zero();
        BC = Eigen::Vector3i::Zero();
        BCval = Eigen::Vector3d::Zero();
        V_eff = 0.0;
        psi = 0.0;

    }

    // Default constructor
    Points() = default;
};


#endif // POINTS_H