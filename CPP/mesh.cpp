#include "mesh.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Robust DOF assignment with comprehensive validation
void AssignDOF(std::vector<Points>& point_list, int PD, int& DOFs, int& DOCs) {
    DOFs = 0;
    DOCs = 0;

    for (auto& point : point_list) {
        point.DOF.setZero();
        point.DOC.setZero();
        point.BC.setZero();

        // Assign DOF for free points (regular points)
        if (point.Flag == "Point") {
            for (int dim = 0; dim < PD; dim++) {
                point.DOF(dim) = ++DOFs;
                point.BC(dim) = 1;  // For free points, BC should be 1
            }
            point.DOFs = PD;
        }
        // Assign DOC for constrained points (fixed boundary conditions)
        else if (point.Flag == "Patch" || point.Flag == "RightPatch") {
            for (int dim = 0; dim < PD; dim++) {
                point.DOC(dim) = ++DOCs;
                point.BC(dim) = 0;
            }
            point.DOCs = PD;
        }

        // Sanity check for each point
        if (point.DOFs + point.DOCs != PD) {
            throw std::runtime_error("Inconsistent DOF/DOC assignment for point");
        }
    }

    // Final validation
    int total_point_dofs = 0;
    for (const auto& point : point_list) {
        total_point_dofs += point.DOFs;
    }

    if (total_point_dofs != DOFs) {
        throw std::runtime_error("Total DOFs do not match point-wise calculation");
    }
}

// More robust volume assignment
void AssignVolumes(std::vector<Points>& point_list, int PD, double Delta) {
    double volume = Delta;
    if (PD == 2) volume = Delta * Delta;
    else if (PD == 3) volume = Delta * Delta * Delta;

    for (auto& point : point_list) {
        point.volume = volume;
    }
}

// Comprehensive boundary condition verification
void VerifyBoundaryConditions(const std::vector<Points>& points, int PD) {
    int constrained_points = 0;
    int free_points = 0;

    for (const auto& p : points) {
        if (p.Flag != "Point") {
            if (p.BC.head(PD).sum() != 0) {
                throw std::runtime_error("Invalid BC for constrained point");
            }
            constrained_points++;
        } else {
            if (p.BC.head(PD).sum() != PD) {
                throw std::runtime_error("Invalid BC for free point");
            }
            free_points++;
        }
    }

    if (constrained_points + free_points != points.size()) {
        throw std::runtime_error("Point type mismatch");
    }
}

// Robust mesh generation with extensive error checking
std::vector<Points> generate_mesh(int PD, double domain_size, int number_of_patches,
                                 double Delta, int number_of_right_patches,
                                 int& DOFs, int& DOCs, const double& d) {
    std::vector<Points> point_list;
    const int number_of_points = std::floor(domain_size / Delta) + 1 ;
    int total_points = (number_of_patches + number_of_right_patches + number_of_points) ;
    int index = 0;

    switch (PD) {
        case 1: {
            for (int i = 0; i < total_points; i++) {
                Points point;
                point.Nr = index++;
                point.X = Eigen::Vector3d(Delta / 2 + i * Delta, 0, 0);

                if (i < number_of_patches) {
                    point.Flag = "Patch";
                    point.x = point.X;
                }
                else if (i >= number_of_patches + number_of_points) {
                    point.Flag = "RightPatch";
                    point.BCval = Eigen::Vector3d(d, 0, 0);
                }
                else {
                    point.Flag = "Point";
                    point.x = point.X;
                }
                point_list.push_back(point);
            }
            break;
        }

        case 2: {
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < total_points; j++) {
                    Points point;
                    point.Nr = index++;
                    point.X = Eigen::Vector3d(Delta / 2 + j * Delta, Delta / 2 + i * Delta, 0);

                    if (j < number_of_patches || j >= number_of_patches + number_of_points) {
                        point.Flag = j < number_of_patches ? "Patch" : "RightPatch";
                        point.x = point.X;
                    } else {
                        point.Flag = "Point";
                    }
                    point_list.push_back(point);
                }
            }
            break;
        }

        case 3: {
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < number_of_points; j++) {
                    for (int k = 0; k < total_points; k++) {
                        Points point;
                        point.Nr = index++;
                        point.X = Eigen::Vector3d(
                            Delta / 2 + k * Delta,
                            Delta / 2 + j * Delta,
                            Delta / 2 + i * Delta
                        );

                        if (k < number_of_patches || k >= number_of_patches + number_of_points) {
                            point.Flag = k < number_of_patches ? "Patch" : "RightPatch";
                            point.x = point.X;
                        } else {
                            point.Flag = "Point";
                        }
                        point_list.push_back(point);
                    }
                }
            }
            break;
        }

        default:
            throw std::runtime_error("Invalid PD value in mesh generation");
    }

    AssignDOF(point_list, PD, DOFs, DOCs);
    AssignVolumes(point_list, PD, Delta);
    VerifyBoundaryConditions(point_list, PD);

    return point_list;
}