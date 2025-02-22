//
// Created by srini on 23/02/2025.
//

#ifndef DEBUG_H
#define DEBUG_H

inline void debug_it(int PD, std::vector<Points>& point_list)
{
    // Mesh Debugging
    for (auto & i : point_list) {
        std::cout << "Nr: " << i.Nr << ", X: [";
        for (const auto& val : i.X) {
            std::cout << val << " ";
        }
        std::cout << "], x: [";
        for (const auto& val : i.x) {
            std::cout << val << " ";
        }
        std::cout << "], Volume: " << i.volume << std::endl;
        std::cout<< "BC: " << i.BC << " Flag: " << i.Flag << std::endl;
    }

    // Neighbour_list debugging
    for (auto & i : point_list)
    {
        std::cout << "Stored neighbours for point " << i.Nr << ": ";

        if(PD == 1)
        {
            for (const auto& n : i.neighbour_list_1N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << std::endl;
            std::cout << "number of neighbours for point " << i.Nr << ": "<< i.n1;
        }
        if(PD == 2)
        {
            for (const auto& n : i.neighbour_list_2N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << std::endl;
            std::cout << "number of neighbours for point " << i.Nr << ": "<< i.n2;
        }
        if(PD == 3)
        {
            for (const auto& n : i.neighbour_list_3N) {
                std::cout << "{";
                for (const int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << std::endl;
            std::cout << "number of neighbours for point " << i.Nr << ": " << i.n3;
        }
        std::cout << std::endl;

    }
}

#endif //DEBUG_H
