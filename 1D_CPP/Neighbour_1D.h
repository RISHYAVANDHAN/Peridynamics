#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <vector>
#include <iostream>
#include <cmath>
#include "Points_1D.h"

/// Checks if two points are within 'delta' distance of each other (Euclidean norm)
inline bool NeighbourCheck(const Points& a, const Points& b, double delta) {
    return std::abs(a.X - b.X) < delta;
}


/// Generates only 1-neighbor lists for 1D points
inline void generate_neighbour_list_1D(std::vector<Points>& point_list, double delta) {
    if (point_list.empty()) return;

    for (size_t i = 0; i < point_list.size(); ++i) {
        // Clear existing neighbor lists
        point_list[i].neighbour_list_1N.clear();
        point_list[i].n1 = 0;

        // 1-Neighbors
        for (size_t j = 0; j < point_list.size(); ++j) {
            if (i != j && NeighbourCheck(point_list[i], point_list[j], delta)) {
                point_list[i].neighbour_list_1N.push_back({static_cast<int>(j)});
                point_list[i].n1++;
            }
        }

        // Optional debug print
        std::cout << "Point " << i << " has " << point_list[i].n1 << " 1N neighbors\n";

        // Populate generic neighbor list for consistency
        point_list[i].neighbour_list = point_list[i].neighbour_list_1N;
    }
}

#endif // NEIGHBOUR_H
