//
// Created by srini on 23/02/2025.
//

#ifndef DEBUG_H
#define DEBUG_H

#include <vector>
#include <iostream>

// Forward declaration or include the Points class definition
class Points;

inline void debug_it(int PD, std::vector<Points>& point_list)
{
    // Mesh Debugging
    for (const auto& i : point_list)
    {
        // Changed to const reference for better performance
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

        std::cout << "Neighbours of " << i.Nr << " are: [";  // Removed extra space

        if (PD == 1) {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\nNumber of neighbours for point " << i.Nr << ": " << i.n1;
        }
        else if (PD == 2) {
            // Display 1-neighbor interactions
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }

            // Display 2-neighbor interactions
            std::cout << " {\n";  // Consistent formatting
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << ", ";  // Consistent comma spacing
                std::cout << "} ";
            }
            std::cout << "}\n";
            std::cout << "Number of 1-neighbours: " << i.n1 << ", Number of 2-neighbours: " << i.n2;
        }
        else if (PD == 3) {
            // Display 1-neighbor interactions
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }

            // Display 2-neighbor interactions
            std::cout << " {\n";
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{ ";
                for (const int val : n) std::cout << val << ", ";
                std::cout << "} ";
            }

            // Display 3-neighbor interactions
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
        std::cout << "]\n";  // Consistent newline usage
        std::cout << std::endl;
    }
}

#endif //DEBUG_H