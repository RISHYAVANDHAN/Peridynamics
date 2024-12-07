//
// Created by srini on 22/11/2024.
//

#ifndef POINTS_H
#define POINTS_H
#include <iostream>
#include <vector>



class Points {
public:
    int Nr{};                                           // Nr - position number or number of points
    std::vector<double> X;                              // X material coordinates
    std::vector<double> x;                              // x spatial coordinates
    double volume{};                                      // Volume
    std::vector<std::vector<int>> point_list;
    std::vector<std::vector<int>> neighbour_list;       // consolidated neighbour_list containing 1/2/3 depending on PD
    int Flag = 0;                                       // Flag (0 for patches, 1 for points)
    int DOC = 0;                                        // temporary names DOC and DOF for global free index and global fixed index
    int DOF = 0;
    int BC = 0;                                         // 0 - patches, 1 - points

    // Constructor for initializing Points object with necessary data
    Points(int Nr, std::vector<double>& X, std::vector<double>& x, double volume )
    : Nr(Nr),
      X(std::move(X)),
      x(x),
      volume(volume)
    {
    }

    Points() = default;

};



#endif //POINTS_H

