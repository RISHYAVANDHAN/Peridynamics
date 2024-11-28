//
// Created by ACER on 27-11-2024.
//

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <vector>

struct Points{
    int Nr;                       // Node ID
    std::vector<double> X;		// Original Coordinates
    std::vector<double> x;		// Deformed Coordinates
    //double Volume;   			//
    int Number_of_Neighbors; 		// Number of Neighbours for a point
    bool Point_BC; 				// 1: Point within domain, 0: Point outside domain
    int DOF;						// Degree-of-freedom for points within domain
    int DOC;						// Degree-of-constraints for points outside domain
};

#endif //POINT_H
