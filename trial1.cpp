//
// Created by ACER on 25-11-2024.
//
#pragma once
#include <iostream>
#include "input_ani.h"
#include "Point_ani.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <vector>
#include <cstring>
#include <algorithm>
#include <array>
#include <unordered_map>
#include "PD -Ref-Code/Eigen/Dense"
#include "PD -Ref-Code/Eigen/Sparse"
#include "PD -Ref-Code/Eigen/IterativeLinearSolvers"
#include <filesystem>


void Generate_Mesh(const InputParams& inputParams, std::vector<Points>& points){
	double domainSize = inputParams.domainSize;
    double Delta = inputParams.gridSpacing;

	int numPointsperDim = static_cast<int>(ceil(domainSize / Delta));
	int index = 0;
    for (int i = 0; i < numPointsperDim; i++){
      for (int j = 0; j < numPointsperDim; j++){
        Points point;
        point.Nr = index;
        point.X = {Delta/2 + i*Delta, Delta/2 +j*Delta};
        point.x = point.X;
        point.Number_of_Neighbors = numPointsperDim;
        index += 1;
        point.Point_BC = true;
        point.DOF = 0;
        point.DOC = 0;
        points.push_back(point);
    };
  };
};

void Assign_BC(std::vector<Points>& points, double delta ){
  double horizon = delta;

}







