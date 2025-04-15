#ifndef MESH_H
#define MESH_H

#include "mesh_1D.h"
#include <iostream>
#include <vector>
#include "Points_1D.h"
#include <Eigen/Dense>

void VerifyBoundaryConditions(const std::vector<Points>& points) {
    int constrained_points = 0;
    int free_points = 0;

    for (const auto& p : points) {
        if (p.Flag == "Patch") {
            if (p.BC != 0) {
                std::cerr << "Invalid BC for Patch point " << p.Nr << "\n";
            }
        } else {
            if (p.BC <= 0) {
                std::cerr << "Invalid BC for Point " << p.Nr << "\n";
            }
        }
    }

    std::cout << "Boundary Condition Report:\n"
              << "- Constrained points: " << constrained_points << "\n"
              << "- Free points: " << free_points << "\n"
              << "- Total DOFs: " << free_points << "\n";
}

double Compute_FF_1D(double d, const std::string& DEFflag) {
    if (DEFflag == "EXT") return 1 + d;
    if (DEFflag == "EXP") return 1 + d;
    return 1.0;
}

void AssignDOF(std::vector<Points>& point_list, int& DOFs) {
    DOFs = 0;
    for (auto& point : point_list) {
        if (point.BC == 1) {
            DOFs++;
            point.DOF = DOFs;
        } else {
            point.DOF = 0;
        }
    }
}

void AssignVolumes(std::vector<Points>& point_list, double Delta) {
    double volume = Delta;
    for (auto& point : point_list) {
        point.volume = volume;
    }
}

std::vector<Points> generate_mesh(double d, double domain_size, int number_of_points,
                                  int number_of_patches, double Delta, int number_of_right_patches,
                                  const std::string& DEFflag, int& DOFs) {
    double extended_domain_size = domain_size + (number_of_patches + number_of_right_patches) * Delta;
    double FF = Compute_FF_1D(d, DEFflag);
    std::vector<Points> point_list;
    int index = 0;

    const int total_points = number_of_patches + number_of_points + number_of_right_patches;
    for (int i = 0; i < total_points; i++) {
        Points point;
        point.Nr = index++;
        point.X = Delta / 2 + i * Delta;   // X is now a double
        point.x = point.X;

        if (i < number_of_patches || i >= number_of_patches + number_of_points) {
            point.BC = 0;
            point.Flag = "Patch";
            point.x = point.X;

            if (i >= number_of_patches + number_of_points) {
                point.BCval = (1 + d) * point.X - point.X;  // 1D deformation
            }
        } else {
            point.BC = 1;
            point.Flag = "Point";
            point.BCval = (1 + d) * point.X - point.X;  // 1D deformation
        }

        point_list.push_back(point);
    }

    AssignDOF(point_list, DOFs);
    AssignVolumes(point_list, Delta);
    return point_list;
}

#endif // MESH_H
