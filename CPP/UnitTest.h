//
// Created by srini on 19/03/2025.
//

#ifndef UNITTEST_H
#define UNITTEST_H

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "solver.h"
#include "Assemble.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>

// Function declarations
void testEnergyCalculation();
void testResidualCalculation();
void testStiffnessCalculation();
void runUnitTests();

#endif // UNITTEST_H