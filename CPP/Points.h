//
// Created by srini on 22/11/2024.
//

#ifndef POINTS_H
#define POINTS_H
#include <iostream>
#include <vector>
#include <utility> //std::pair - just fyi
#include <tuple>



class Points {
public:
    int Nr{};                                           // Nr - position number or number of points
    std::vector<double> X;                              // X material coordinates
    std::vector<double> x;                              // x spatial coordinates
    double volume{};                                    // Volume
    std::vector<std::vector<int>> point_list;
    std::vector<std::vector<int>> neighbour_list_1N;    // consolidated neighbour_list containing 1-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_2N;    // consolidated neighbour_list containing 2-neighbour interaction
    std::vector<std::vector<int>> neighbour_list_3N;    // consolidated neighbour_list containing 3-neighbour interaction
    std::string Flag = "Patch";                         // Flag (Patches & points)
    int DOC = 0;                                        // temporary names DOC and DOF for global free index and global fixed index
    int DOF = 0;
    int BC = 0;                                         // 0 - patches, 1 - points
    size_t n1 = 0;                                         // n1 - no. of 1 neighbour interaction
    size_t n2 = 0;                                         // n2 - no. of 2 neighbour interaction
    size_t n3 = 0;                                         // n3 - no. of 3 neighbour interaction

    // Constructor for initializing Points object with necessary data
    Points(int Nr, const std::vector<double>& X, const std::vector<double>& x, double volume)
    : Nr(Nr),
      X(X.size() == 3 ? X : std::vector<double>(3, 0.0)),
      x(x.size() == 3 ? x : std::vector<double>(3, 0.0)),
      volume(volume),
      neighbour_list_1N(),
      neighbour_list_2N(),
      neighbour_list_3N()
    {
        if (X.size() != 3 || x.size() != 3) {
            std::cerr << "Warning: Input vectors must be size 3. Initializing with zeros.\n";
        }
    }

    Points() = default;

};



#endif //POINTS_H

