#ifndef POINTS_H
#define POINTS_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class Points {
public:
    int Nr{};                                               // Point number
    Eigen::Vector3d X = Eigen::Vector3d::Zero();            // Material coordinates
    Eigen::Vector3d x = Eigen::Vector3d::Zero();            // Spatial coordinates
    Eigen::Vector3d x_prev = Eigen::Vector3d::Zero();       // Previous spatial coordinates
    double volume{};                                        // Volume
    std::vector<std::vector<int>> neighbour_list_1N;        // 1-neighbor interactions
    std::vector<std::vector<int>> neighbour_list_2N;        // 2-neighbor interactions
    std::vector<std::vector<int>> neighbour_list_3N;        // 3-neighbor interactions
    std::vector<std::vector<int>> neighbour_list;           // All interactions
    Eigen::Vector3d Ra_1 = Eigen::Vector3d::Zero();         // 1-neighbor residual
    Eigen::Vector3d Ra_2 = Eigen::Vector3d::Zero();         // 2-neighbor residual
    Eigen::Vector3d Ra_3 = Eigen::Vector3d::Zero();         // 3-neighbor residual
    Eigen::Vector3d Ra_sum = Eigen::Vector3d::Zero();       // Total residual
    Eigen::Matrix3d Kab_1 = Eigen::Matrix3d::Zero();        // 1-neighbor stiffness
    Eigen::Matrix3d Kab_2 = Eigen::Matrix3d::Zero();        // 2-neighbor stiffness
    Eigen::Matrix3d Kab_3 = Eigen::Matrix3d::Zero();        // 3-neighbor stiffness
    Eigen::Matrix3d Kab_sum = Eigen::Matrix3d::Zero();      // Total stiffness
    std::string Flag = "Patch";                             // Point type flag
    Eigen::Vector3i DOF = Eigen::Vector3i::Zero();          // Degrees of freedom
    Eigen::Vector3i DOC = Eigen::Vector3i::Zero();          // Degrees of constraint
    Eigen::Vector3i BC = Eigen::Vector3i::Zero();           // Boundary condition
    Eigen::Vector3d BCval = Eigen::Vector3d::Zero();        // BC values
    int DOFs = 0;                                           // Number of free DOFs
    int DOCs = 0;                                           // Number of constrained DOFs
    size_t n1 = 0;                                          // 1-neighbor count
    size_t n2 = 0;                                          // 2-neighbor count
    size_t n3 = 0;                                          // 3-neighbor count
    double V_eff = 0.0;                                     // Effective volume
    double psi = 0.0;                                       // Energy density

    // Constructor
    Points(int Nr = 0, const Eigen::Vector3d& X = Eigen::Vector3d::Zero(),
           const Eigen::Vector3d& x = Eigen::Vector3d::Zero(), double volume = 0.0)
        : Nr(Nr), X(X), x(x), volume(volume) {
        // Initialize all member variables
        Ra_1.setZero();
        Ra_2.setZero();
        Ra_3.setZero();
        Ra_sum.setZero();
        Kab_1.setZero();
        Kab_2.setZero();
        Kab_3.setZero();
        Kab_sum.setZero();
        DOF.setZero();
        DOC.setZero();
        BC.setZero();
        BCval.setZero();
    }
};

#endif // POINTS_H