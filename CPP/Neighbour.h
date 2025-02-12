//
// Created by srini on 28/11/2024.
//


//Debugging reference for future, the idx for 2d with domain 10.0, Delta 1.0, delta 3.0 is 2676

#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

//if we want to use number_of_points ass the mesh.h header beforehand
//std::vector<int> neighbour_list_1;                  // 1 Neighbour
//std::vector<std::array<int, 2>> neighbour_list_2;   // 2 Neighbours
//std::vector<std::array<int, 3>> neighbour_list_3;   // 3 Neighbours
// 3 separate lists are not needed, 1 is just enough, it can have everything in it and this can be easily parallelized with pragma openmp

inline bool NeighbourCheck(const Points& a, const Points& b, double delta)
{
    double diff_x = a.X[0] - b.X[0];
    double diff_y = a.X[1] - b.X[1];
    double diff_z = a.X[2] - b.X[2];
    return sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z) < delta;
}

inline std::vector<std::vector<int>> generate_neighbour_list(int PD, std::vector<Points>& point_list, int number_of_patches, int number_of_right_patches, int delta) {
    std::vector<std::vector<int>> potential_nbrs;
    std::cout << "Generating neighbour list..." << std::endl;
    Points points;

    int idx = 0; // just another debugging variable

    // comparing the distance of the points from the point in question in 1d, 2d and 3d to find the neighbours and updating the neighbour_list
    for (size_t i = 0; i < point_list.size(); ++i) {
        // Debugging
        std::cout << "Neighbours of " << point_list[i].Nr << " are: " << std::endl << "[ ";
        std::set<std::vector<int>> unique_triplets;  // Set to store unique triplets
        for (size_t j = 0; j < point_list.size(); ++j) {
            if (i != j && NeighbourCheck(point_list[i], point_list[j], delta)) {
                potential_nbrs.push_back({point_list[j].Nr}); // 1 neighbour interaction = also this is it for 1D
                std::cout << "{" << point_list[j].Nr << "}, ";
                idx++;
                for (size_t k = 0; k < point_list.size(); ++k) {
                    if (k != j && k != i && PD != 1 && NeighbourCheck(point_list[k], point_list[j], delta) && NeighbourCheck(point_list[k], point_list[i], delta)) {
                        std::vector<int> pair = {std::min(point_list[k].Nr, point_list[j].Nr), std::max(point_list[k].Nr, point_list[j].Nr)};
                        // Use std::any_of to check for duplicate pairs, and capture pair by reference
                        if (!std::ranges::any_of(potential_nbrs, [&pair](const std::vector<int>& existing_pair) { return pair == existing_pair; })) {
                            potential_nbrs.push_back(pair); // 2 neighbour interaction = also this is it for 2D
                            std::cout << " {" << pair[0] << ", " << pair[1] << "}, ";
                        }

                        for (size_t l = 0; l < point_list.size(); ++l) {
                            if (l != k && l != j && l != i && PD != 2 && NeighbourCheck(point_list[l], point_list[k], delta) && NeighbourCheck(point_list[l], point_list[j], delta) && NeighbourCheck(point_list[l], point_list[i], delta)) {
                                std::vector<int> triplet = {std::min({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr}), std::min({std::max(point_list[k].Nr, point_list[j].Nr), point_list[l].Nr}),
                                                            std::max({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr})};
                                if (!unique_triplets.contains(triplet)) {
                                    unique_triplets.insert(triplet); // Only insert unique triplets
                                    potential_nbrs.push_back(triplet); // 3 neighbour interaction
                                    std::cout << " {" << triplet[0] << ", " << triplet[1] << ", " << triplet[2] << "}, ";
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << " ]" << std::endl;
        points.neighbour_list.insert(points.neighbour_list.end(), potential_nbrs.begin(), potential_nbrs.end());
    }

    std::cout << "idx: " << idx << std::endl;
    return points.neighbour_list;
}

#endif //NEIGHBOUR_H
