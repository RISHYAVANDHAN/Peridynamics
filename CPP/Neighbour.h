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

inline bool NeighbourCheck(const Points& a, const Points& b, double delta)
{
    double diff_x = a.X[0] - b.X[0];
    double diff_y = a.X[1] - b.X[1];
    double diff_z = a.X[2] - b.X[2];
    return sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z) < delta;
}

inline std::vector<std::vector<int>> generate_neighbour_list(int PD, std::vector<Points>& point_list, int delta) {
    //std::cout << "Generating neighbour list..." << std::endl;

    if (point_list.empty()) {
        std::cerr << "Point list is empty. Exiting..." << std::endl;
        return {};
    }

    for (size_t i = 0; i < point_list.size(); ++i) {
        std::set<std::vector<int>> unique_pairs;
        std::set<std::vector<int>> unique_triplets;

        //std::cout << "Neighbours of " << point_list[i].Nr << " are: [ ";
        int n1 = 0, n2 = 0, n3 = 0;

        switch (PD) {
            case 1: {
                // Clear previous neighbors if any
                point_list[i].neighbour_list_1N.clear();

                for (size_t j = 0; j < point_list.size(); ++j) {
                    if (i != j && NeighbourCheck(point_list[i], point_list[j], delta)) {
                        std::vector<int> single = {point_list[j].Nr};
                        point_list[i].neighbour_list_1N.push_back(single);
                        //std::cout << point_list[j].Nr << " ";
                        n1++;
                    }
                }
                break;
            }

            case 2: {
                // Clear previous neighbors if any
                point_list[i].neighbour_list_2N.clear();

                for (size_t j = 0; j < point_list.size(); ++j) {
                    for (size_t k = 0; k < point_list.size(); ++k) {
                        if (k != j && k != i &&
                            NeighbourCheck(point_list[k], point_list[j], delta) &&
                            NeighbourCheck(point_list[k], point_list[i], delta)) {

                            std::vector<int> pair = {
                                std::min(point_list[k].Nr, point_list[j].Nr),
                                std::max(point_list[k].Nr, point_list[j].Nr)
                            };

                            if (unique_pairs.insert(pair).second) {
                                point_list[i].neighbour_list_2N.push_back(pair);
                                //std::cout << "{" << pair[0] << ", " << pair[1] << "} ";
                                n2++;
                            }
                        }
                    }
                }
                break;
            }

            case 3: {
                // Clear previous neighbors if any
                point_list[i].neighbour_list_3N.clear();

                for (size_t j = 0; j < point_list.size(); ++j) {
                    for (size_t k = 0; k < point_list.size(); ++k) {
                        for (size_t l = 0; l < point_list.size(); ++l) {
                            if (l != k && l != j && l != i &&
                                NeighbourCheck(point_list[l], point_list[k], delta) &&
                                NeighbourCheck(point_list[l], point_list[j], delta) &&
                                NeighbourCheck(point_list[l], point_list[i], delta)) {

                                std::vector<int> triplet = {
                                    std::min({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr}),
                                    std::min({std::max(point_list[k].Nr, point_list[j].Nr), point_list[l].Nr}),
                                    std::max({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr})
                                };

                                if (unique_triplets.insert(triplet).second) {
                                    point_list[i].neighbour_list_3N.push_back(triplet);
                                    //std::cout << "{" << triplet[0] << ", " << triplet[1] << ", " << triplet[2] << "} ";
                                    n3++;
                                }
                            }
                        }
                    }
                }
                break;
            }
        }

        //std::cout << "]" << std::endl;
        point_list[i].n1 = n1;
        point_list[i].n2 = n2;
        point_list[i].n3 = n3;
        //std::cout << "N1: " << n1 << " N2: " << n2 << " N3: " << n3 << std::endl;
    }

    return {};  // The function modifies point_list directly
}

#endif //NEIGHBOUR_H



/*
 for (size_t i = 0; i < point_list.size(); ++i) {
        int n1 = 0, n2 = 0, n3 = 0;
        std::cout << "Neighbours of " << point_list[i].Nr << " are: " << std::endl << "[ ";
        std::set<std::vector<int>> unique_triplets;  // Set to store unique triplets
        for (size_t j = 0; j < point_list.size(); ++j) {
            if (i != j && NeighbourCheck(point_list[i], point_list[j], delta)) {
                potential_nbrs.push_back({point_list[j].Nr}); // 1 neighbour interaction = also this is it for 1D
                std::cout << "{" << point_list[j].Nr << "}, ";
                n1++;
                idx++;
                for (size_t k = 0; k < point_list.size(); ++k) {
                    if (k != j && k != i && PD != 1 && NeighbourCheck(point_list[k], point_list[j], delta) && NeighbourCheck(point_list[k], point_list[i], delta)) {
                        std::vector<int> pair = {std::min(point_list[k].Nr, point_list[j].Nr), std::max(point_list[k].Nr, point_list[j].Nr)};
                        // Use std::any_of to check for duplicate pairs, and capture pair by reference
                        if (!std::ranges::any_of(potential_nbrs, [&pair](const std::vector<int>& existing_pair) { return pair == existing_pair; })) {
                            potential_nbrs.push_back(pair); // 2 neighbour interaction = also this is it for 2D
                            std::cout << " {" << pair[0] << ", " << pair[1] << "}, ";
                            n2++;
                        }

                        for (size_t l = 0; l < point_list.size(); ++l) {
                            if (l != k && l != j && l != i && PD != 2 && NeighbourCheck(point_list[l], point_list[k], delta) && NeighbourCheck(point_list[l], point_list[j], delta) && NeighbourCheck(point_list[l], point_list[i], delta)) {
                                std::vector<int> triplet = {std::min({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr}), std::min({std::max(point_list[k].Nr, point_list[j].Nr), point_list[l].Nr}),
                                                            std::max({point_list[k].Nr, point_list[j].Nr, point_list[l].Nr})};
                                if (!unique_triplets.contains(triplet)) {
                                    unique_triplets.insert(triplet); // Only insert unique triplets
                                    potential_nbrs.push_back(triplet); // 3 neighbour interaction
                                    std::cout << " {" << triplet[0] << ", " << triplet[1] << ", " << triplet[2] << "}, ";
                                    n3++;
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << " ]" << std::endl;
        point_list[i].n1 = n1;
        point_list[i].n2 = n2;
        point_list[i].n3 = n3;
        point_list[i].neighbour_list.insert(point_list[i].neighbour_list.end(), potential_nbrs.begin(), potential_nbrs.end());
    }

 */

