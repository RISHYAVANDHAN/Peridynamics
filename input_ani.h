//
// Created by ACER on 27-11-2024.
//

#ifndef INPUT_H
#define INPUT_H

#include <iostream>

// Material constants and input parameters
struct InputParams {
    int Problem_dimension;      // 1D, 2D, 3D
    double domainSize;          // Size of the domain
    double gridSpacing;         // Grid spacing (Delta)
    double delta; 				// Horizon size
    double horizonRatio;        // Horizon size / grid spacing (delta/Delta) Number of patches needed
    double deformationMag; 		// Deformation magnitude
    int numSteps;               // Number of steps for deformation
    double C1;                  // Material constant resistance
    double NN;                  // Material power law
    char flag;                  // Displacement ('d') or force ('f')
    int deformationfunction;    // 0=constant, 1=linear, etc.
};

#endif //INPUT_H
