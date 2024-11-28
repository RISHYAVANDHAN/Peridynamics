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
    int Number_of_points;                               // Nr - position number or number of points
    std::vector<double> X;                              // X material coordinates
    std::vector<double> x;                              // x spatial coordinates
    double volume;                                      // Volume
    std::vector<int> point_index;
    int Number_of_neighbours;                           // Number of neighbours (for peridynamics or mesh)
    std::vector<int> point_list;
    std::vector<int> neighbour_list_1;                  // 1 Neighbour
    std::vector<std::array<int, 2>> neighbour_list_2;   // 2 Neighbours
    std::vector<std::array<int, 3>> neighbour_list_3;   // 3 Neighbours
    std::vector<double> neighbour_index;
    int Flag = 0;                                       // Flag (0 for patches, 1 for points)
    int DOC = 0;                                        // temporary names DOC and DOF for global free index and global fixed index
    int DOF = 0;

    // Constructor for initializing Points object with necessary data
    Points(int Number_of_points, int Number_of_neighbours, const std::vector<double>& X, const std::vector<double>& x,
           double volume, const std::vector<int>& point_index, const std::vector<std::array<int, 2>>& neighbour_list_2,
           const std::vector<std::array<int, 3>>& neighbour_list_3, const std::vector<int>& neighbour_list_1, const std::vector<double>& neighbour_index)
        : Number_of_points(Number_of_points), X(X), x(x), volume(volume), point_index(point_index), Number_of_neighbours(Number_of_neighbours),
          neighbour_list_1(neighbour_list_1), neighbour_list_2(neighbour_list_2), neighbour_list_3(neighbour_list_3), neighbour_index(neighbour_index)
    {
    }
    Points() = default;

};


#endif //POINTS_H

