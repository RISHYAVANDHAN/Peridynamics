//
// Created by srini on 22/11/2024.
//

#ifndef POINTS_H
#define POINTS_H
#include <iostream>
#include <vector>
#include <array>



class Points {
public:
    int Number_of_points{};                             // Nr - position number or number of points
    std::vector<double> X;                              // X material coordinates
    std::vector<double> x;                              // x spatial coordinates
    double volume{};                                    // Volume
    std::vector<int> point_index;
    int Number_of_neighbours{};                         // Number of neighbours (for peridynamics or mesh)
    std::vector<int> point_list;
    std::vector<int> neighbour_list_1;                  // 1 Neighbour
    std::vector<std::array<int, 2>> neighbour_list_2;   // 2 Neighbours
    std::vector<std::array<int, 3>> neighbour_list_3;   // 3 Neighbours
    std::vector<double> neighbour_index;
    int Flag = 0;                                       // Flag (0 for patches, 1 for points)
    int DOC = 0;                                        // temporary names DOC and DOF for global free index and global fixed index
    int DOF = 0;
    int BC = 0;                                         // 0 - patches, 1 - points

    // Constructor for initializing Points object with necessary data
    explicit Points(const int Number_of_points)
    : Number_of_points(Number_of_points),
      X(Number_of_points, 0.0),
      x(Number_of_points, 0.0),
      volume(1.0),
      point_index(Number_of_points, 0),
      Number_of_neighbours(0),
      point_list(Number_of_points, 0),
      neighbour_list_1(Number_of_points, 0),
      neighbour_list_2(Number_of_points, {0, 0}),
      neighbour_list_3(Number_of_points, {0, 0, 0}),
      neighbour_index(Number_of_points, 0.0)
      {
      }

    Points() = default;

};



#endif //POINTS_H

