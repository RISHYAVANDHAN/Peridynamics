//
// Created by srini on 28/11/2024.

#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <iostream>
#include <memory>
#include <vector>


//std::vector<int> neighbour_list_1;                  // 1 Neighbour
//std::vector<std::array<int, 2>> neighbour_list_2;   // 2 Neighbours
//std::vector<std::array<int, 3>> neighbour_list_3;   // 3 Neighbours
// 3 seperate lists are not needed, 1 is just enough, it can have everything in it and this can be easily parallelized with pragma openmp
inline std::vector<Points> working_list;
inline std::vector<std::vector<int>> potential_nbrs;// Temporary list to hold potential neighbours

inline void generate_neighbour_list(int PD, int number_of_points, std::vector<Points>& point_list, int number_of_patches, int number_of_right_patches, int delta) {
    std::cout << "Generating neighbour list..." << std::endl;

    int total_points = (number_of_patches + number_of_right_patches + number_of_points);
    for (auto& i : point_list) {
        for(int j = 0; j < (total_points*total_points); j++){
            if (i.Nr > ((j * total_points) + (number_of_patches - 1)) && (i.Nr < (number_of_patches + number_of_points + (j * total_points)))) {
                std::cout << i.Nr << std::endl;
                working_list.push_back(i);
            }
        }
    }

    // Debugging to check if working list is what we want
    std::cout << "Working list: " << std::endl;
    for (auto& i : working_list) {
        std::cout << "Nr: " << i.Nr << ", X: [";
        for (const auto& val : i.X) {
            std::cout << val << " ";
        }
        std::cout << "], x: [";
        for (const auto& val : i.x) {
            std::cout << val << " ";
        }
        std::cout << "], Volume: " << i.volume << std::endl;
    }


    switch (PD) {
        case 1:

            for (auto &i : working_list) {
                for (auto &j : working_list) {
                    if ((i.Nr != j.Nr) && ((std::abs(i.X[0] - j.X[0]) < static_cast<double>(delta)))) {
                        potential_nbrs.push_back({j.Nr});
                    }
                }
            }
        break;
        case 2:
            for (auto &i : working_list) {
                for (auto &j : working_list) {
                    if ((i.Nr != j.Nr) && ((std::abs(i.X[0] - j.X[0]) < static_cast<double>(delta))) && ((std::abs(i.X[1] - j.X[1]) < static_cast<double>(delta)))) {
                        potential_nbrs.push_back({j.Nr});
                        for (auto &k : working_list) {
                            if ((k.Nr != j.Nr) && ((std::abs(k.X[0] - j.X[0]) < static_cast<double>(delta))) && ((std::abs(k.X[1] - j.X[1]) < static_cast<double>(delta))) && ((std::abs(k.X[0] - i.X[0]) < static_cast<double>(delta))) && ((std::abs(k.X[1] - i.X[1]) < static_cast<double>(delta))) ) {
                                potential_nbrs.push_back({k.Nr, j.Nr});
                            }
                        }
                    }
                }
            }
        break;
        case 3:
            for (auto &i : working_list) {
                for (auto &j : working_list) {
                    if ((i.Nr != j.Nr) && ((std::abs(i.X[0] - j.X[0]) < static_cast<double>(delta))) && ((std::abs(i.X[1] - j.X[1]) < static_cast<double>(delta))) && ((std::abs(i.X[2] - j.X[2]) < static_cast<double>(delta)))) {
                        potential_nbrs.push_back({j.Nr});
                        for (auto &k : working_list) {
                            if ((k.Nr != j.Nr) &&
                                ((std::abs(k.X[0] - j.X[0]) < static_cast<double>(delta))) && ((std::abs(k.X[1] - j.X[1]) < static_cast<double>(delta))) && ((std::abs(k.X[2] - i.X[2]) < static_cast<double>(delta))) &&
                                ((std::abs(k.X[0] - i.X[0]) < static_cast<double>(delta))) && ((std::abs(k.X[1] - i.X[1]) < static_cast<double>(delta))) && ((std::abs(k.X[2] - j.X[2]) < static_cast<double>(delta))) ) {
                                potential_nbrs.push_back({k.Nr, j.Nr});
                                for (auto &l : working_list) {
                                    if ((l.Nr != k.Nr) &&
                                    ((std::abs(l.X[0] - k.X[0]) < static_cast<double>(delta))) && ((std::abs(l.X[1] - k.X[1]) < static_cast<double>(delta))) && ((std::abs(l.X[2] - k.X[2]) < static_cast<double>(delta))) &&
                                    ((std::abs(l.X[0] - j.X[0]) < static_cast<double>(delta))) && ((std::abs(l.X[1] - j.X[1]) < static_cast<double>(delta))) && ((std::abs(l.X[2] - j.X[2]) < static_cast<double>(delta))) &&
                                    ((std::abs(l.X[0] - i.X[0]) < static_cast<double>(delta))) && ((std::abs(l.X[1] - i.X[1]) < static_cast<double>(delta))) && ((std::abs(l.X[2] - i.X[2]) < static_cast<double>(delta)))){
                                        potential_nbrs.push_back({l.Nr, k.Nr, j.Nr});
                                    }
                                }
                            }
                        }
                    }
                }
            }
        break;
        default:
            std::cout<<"Potential neighbour list is empty"<<std::endl;
        break;
    }


    // Debugging to check potential neighbours
    int idx = 0;
    for (const auto& nbrs : potential_nbrs) {
        std::cout << "{";
        for (const auto& nbr : nbrs) {
            std::cout << nbr << ",  ";
        }
        idx+=1;
        std::cout << "} ";
    }
    std::cout << std::endl;
    std::cout<<"idx = "<<idx << std::endl;
}


#endif //NEIGHBOUR_H
