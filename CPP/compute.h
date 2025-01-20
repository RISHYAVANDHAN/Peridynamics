//
// Created by srini on 18/01/2025.
//

#ifndef COMPUTE_H
#define COMPUTE_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"

double Norm2(std::array<double, 2> &p)
{
    return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}
double Norm2(double x, double y)
{
    return std::sqrt(x * x + y * y);
}

// calculate tangent


#endif //COMPUTE_H
