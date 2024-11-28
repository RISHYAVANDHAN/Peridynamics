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
generate_mesh(int PD, int number_of_points, int Partition, int degree, double domain_size, int number_of_patches, int BC, int Delta, int number_of_neighbours) {
    std::cout << "Generating mesh..." << std::endl;

    int free_points = Partition;
    double dx = domain_size / (degree * free_points);
    double X = 0.0;
    double x = 0.0;
    std::vector point_list(number_of_points, 0);
    std::vector<int> neighbour_list_1(number_of_points, 0);
    std::vector<std::array<int, 2>> neighbour_list_2(number_of_points, {0,0});
    std::vector<std::array<int, 3>> neighbour_list_3(number_of_points, {0,0,0});

    std::vector<Points> points;
    Points point(number_of_points);
    int total_points = ((number_of_patches * 2) + number_of_points); // patches + points
    int index = 0;

    switch (PD) {
        case 1:

            for (int i = 0; i < total_points; i++) {

                if ((i <= number_of_patches) || (i >= (total_points - number_of_patches))) {
                    BC = 1;
                }
                else {
                    BC = 0;
                }
                for (int j = 0; j < total_points; j++) {
                    point.Nr = index;
                    point.X = {X + Delta/2 + i*Delta, Delta/2 +j*Delta};
                    point.x = point.X;
                    index += 1;
                    points.push_back(point);
                }


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
