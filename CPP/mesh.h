#ifndef MESH_H
#define MESH_H

#include "Points.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdexcept>

// Compute deformation gradient matrix
Eigen::Matrix3d Compute_FF(int PD, double d, const std::string& DEFflag);

// Assign degrees of freedom and constraints
void AssignDOF(std::vector<Points>& point_list, int PD, int& DOFs, int& DOCs);

// Assign volumes to points
void AssignVolumes(std::vector<Points>& point_list, int PD, double Delta);

// Generate the computational mesh
std::vector<Points> generate_mesh(int PD, double domain_size, int number_of_patches, double Delta, int number_of_right_patches, int& DOFs, int& DOCs, const double& d);

// Verify boundary conditions
void VerifyBoundaryConditions(const std::vector<Points>& points, int PD);

#endif // MESH_H