//
// Created by srini on 24/11/2024.
//

#ifndef MESH_H
#define MESH_H

#include "Points.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <Eigen/Dense>

// Function to compute the deformation gradient (FF) based on the deformation flag
Eigen::MatrixXd Compute_FF(int PD, double d, const std::string& DEFflag) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(PD, PD); // Identity matrix
    Eigen::MatrixXd FF = Eigen::MatrixXd::Zero(PD, PD); // Initialize FF as zero matrix

    if (DEFflag == "EXT") {
        FF = I;
        FF(0, 0) = 1 + d; // Set (0, 0) entry as 1 + d
    } else if (DEFflag == "EXP") {
        FF = (1 + d) * I; // Scale identity matrix by (1 + d)
    } else if (DEFflag == "SHR") {
        FF = I;
        FF(1, 0) = d; // Set (1, 0) entry as d
    }

    return FF;
}

// Function to assign global degrees of freedom (DOFs) to points
void AssignDOF(std::vector<Points>& point_list, int PD, int DOFs) {
    for (size_t i = 0; i < point_list.size(); i++) {
        // Loop through each dimension
        for (int p = 1; p <= PD; ++p) {
            // Check if the boundary condition flag is active
            if (point_list[i].BC(p) == 1) { // Assuming BC is an integer flag
                // Increment the DOFs counter and assign it to the current DOF
                DOFs++;
                point_list[i].DOF[p] = DOFs;
            }
        }
    }
}

// Function to generate a mesh based on input parameters
std::vector<Points> generate_mesh(int PD, double d, double domain_size, int number_of_points, int number_of_patches, double Delta, int number_of_right_patches, const std::string& DEFflag, const int DOFs) {
    // Calculate the extended domain size
    double extended_domain_size = domain_size + (number_of_patches + number_of_right_patches) * Delta;
    std::cout << "Domain size: " << domain_size << " & Extended Domain size: " << extended_domain_size << std::endl;

    // Calculate the total number of points
    int total_points = (number_of_patches + number_of_right_patches + number_of_patches);
    int number_of_points = total_points - (number_of_right_patches + number_of_patches);
    std::cout << "Number of Points: " << number_of_points << std::endl;

    // Compute the deformation gradient (FF)
    Eigen::MatrixXd FF = Compute_FF(PD, d, DEFflag);

    // Initialize the point list
    std::vector<Points> point_list;
    int index = 0;

    // Generate the mesh based on the number of dimensions (PD)
    switch (PD) {
        case 1: // 1D mesh
            for (int i = 0; i < total_points; i++) {
                Points point;
                point.Nr = index;
                point.X = Eigen::Vector3d(Delta / 2 + i * Delta, 0, 0);
                point.x = point.X;
                index += 1;

                // Determine if the point is a patch or a point
                if ((index < number_of_patches) || (index > number_of_patches + number_of_points - 1)) {
                    point.BC = Eigen::Vector3d::Zero();
                    point.Flag = "Patch";
                } else {
                    point.BC = {1, 0, 0};
                    point.Flag = "Point";
                    point.BCval = FF * point.X - point.X;
                }
                point_list.push_back(point);
                AssignDOF(point_list, PD, DOFs);
            }
            break;

        case 2: // 2D mesh
            std::cout << "2D implementation trial" << std::endl;
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < total_points; j++) {
                    Points point;
                    point.Nr = index;
                    point.X = Eigen::Vector3d(Delta / 2 + j * Delta, Delta / 2 + i * Delta, 0);
                    point.x = point.X;

                    // Determine if the point is a patch or a point
                    if (j < number_of_patches || j >= number_of_patches + number_of_points) {
                        point.BC = Eigen::Vector3d::Zero();
                        point.Flag = "Patch";
                    } else {
                        point.BC = {1, 1, 0};
                        point.Flag = "Point";
                        point.BCval = FF * point.X - point.X;
                    }
                    index += 1;
                    point_list.push_back(point);
                    AssignDOF(point_list, PD, DOFs);
                }
            }
            break;

        case 3: // 3D mesh
            std::cout << "3D implementation trial" << std::endl;
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < number_of_points; j++) {
                    for (int k = 0; k < total_points; k++) {
                        Points point;
                        point.Nr = index;
                        point.X = Eigen::Vector3d(Delta / 2 + k * Delta, Delta / 2 + j * Delta, Delta / 2 + i * Delta);
                        point.x = point.X;
                        index += 1;
                        point_list.push_back(point);
                        AssignDOF(point_list, PD, DOFs);
                    }
                }
            }
            break;

        default:
            std::cerr << "Invalid PD value. Mesh generation aborted." << std::endl;
            break;
    }

    return point_list;
}



#endif
//MESH_H