#ifndef MESH_H
#define MESH_H

#include "Points_1D.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

// Function declarations
Eigen::Matrix3d Compute_FF(int PD, double d, const std::string& DEFflag);
void AssignDOF(std::vector<Points>& point_list, int PD, int& DOFs);
void AssignVolumes(std::vector<Points>& point_list, int PD, double Delta);
std::vector<Points> generate_mesh(double d, double domain_size, int number_of_points,
                                  int number_of_patches, double Delta, int number_of_right_patches,
                                  const std::string& DEFflag, int& DOFs);


// Add the function declaration for VerifyBoundaryConditions
void VerifyBoundaryConditions(const std::vector<Points>& points, int PD);

#endif // MESH_H