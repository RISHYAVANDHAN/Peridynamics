//
// Created by srini on 23/02/2025.
//

#ifndef DEBUG_H
#define DEBUG_H

#include <vector>
#include <iostream>

void debug_it(int PD, std::vector<Points>& point_list) {
    for (const auto& i : point_list) {
        std::cout << "Nr: " << i.Nr << ", X: [";
        for (const auto& val : i.X) {
            std::cout << val << " ";
        }
        std::cout << "], x: [";
        for (const auto& val : i.x) {
            std::cout << val << " ";
        }
        std::cout << "], Volume: " << i.volume << std::endl;
        std::cout << "BC: " << i.BC << " Flag: " << i.Flag << std::endl;

        std::cout << "Neighbours of " << i.Nr << " are: [";
        if (PD == 1) {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\nNumber of neighbours for point " << i.Nr << ": " << i.n1;
        } else if (PD == 2) {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << " {\n";
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << ", ";
                std::cout << "} ";
            }
            std::cout << "}\n";
            std::cout << "Number of 1-neighbours: " << i.n1 << ", Number of 2-neighbours: " << i.n2;
        } else if (PD == 3) {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << " {\n";
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << ", ";
                std::cout << "} ";
            }
            std::cout << "}\n{\n";
            for (const auto& n : i.neighbour_list_3N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << ", ";
                std::cout << "} ";
            }
            std::cout << "}\n";
            std::cout << "Number of 1-neighbours: " << i.n1 << ", Number of 2-neighbours: " << i.n2
                      << ", Number of 3-neighbours: " << i.n3;
        }
        std::cout << "]\n";
        std::cout << std::endl;
    }
}



#endif //DEBUG_H