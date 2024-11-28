//
// Created by srini on 24/11/2024.
//

#ifndef MESH_H
#define MESH_H

#include "Points.h"
#include <iostream>
#include <vector>
#include <array>


inline std::tuple<std::vector<int>, std::vector<int>, std::vector<std::array<int, 2>>, std::vector<std::array<int, 3>>>
generate_mesh(int PD, int number_of_points, int Partition, int degree, double domain_size) {
    std::cout << "Generating mesh..." << std::endl;

    int free_points = Partition;
    double dx = domain_size / (degree * free_points);
    double X = 0.0;
    double x = 0.0;
    std::vector point_list(number_of_points, 0);
    std::vector<int> neighbour_list_1(number_of_points);
    std::vector<std::array<int, 2>> neighbour_list_2(number_of_points);
    std::vector<std::array<int, 3>> neighbour_list_3(number_of_points);

    std::vector<Points> point;

    switch (PD) {
        case 1:

            for (int i = 0; i < number_of_points; i++) {
                point.push_back(Points(number_of_points, Number_of_neighbours, X, x, volume, point_index, neighbour_list_2, neighbour_list_3, neighbour_list_1, neighbour_index));
                X = X + (Delta/2 + (i * Delta)) ;
                x = X;
                point_index+=1;




            }
        break;

        default:
            std::cout << "Unsupported PD type" << std::endl;
        break;
    }

    // Return the necessary data: point_list, neighbour_list_1, neighbour_list_2, neighbour_list_3
    return std::make_tuple(point_list, neighbour_list_1, neighbour_list_2, neighbour_list_3);
}







#endif //MESH_H
