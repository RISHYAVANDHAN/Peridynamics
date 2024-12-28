//
// Created by srini on 24/11/2024.
//

#ifndef MESH_H
#define MESH_H

#include "Points.h"
#include <iostream>
#include <vector>

inline int number_of_points;
std::vector<Points> generate_mesh(int PD, int Partition, int degree, double domain_size, int number_of_patches, double Delta, int number_of_right_patches) {
    std::cout << "Generating mesh..." << std::endl;

    //int free_points = Partition;
    //double dx = domain_size / (degree * free_points);
    std::vector<double> X (3, 0.0) ;
    std::vector<double> x (3, 0.0);
    double volume = 1.0;
    double extended_domain_size = domain_size + (number_of_patches + number_of_right_patches) * Delta ;
    std::cout << "Domain size: " << domain_size << " & Extended Domain size: " << extended_domain_size << std::endl;
    int total_points = extended_domain_size * Delta; // patches + points

    std::vector<Points> point_list;
    number_of_points = total_points - (number_of_right_patches + number_of_patches);
    std::cout << "Number of Points:" << number_of_points << std::endl;
    Points point(number_of_points, X, x , volume);


    int index = 0;

    switch (PD) {
        case 1:
            for (int i = 0; i < total_points; i++) {
                point.Nr = index;
                point.X = { Delta/2 + i*Delta,0 , 0};
                point.x = point.X;
                index += 1;
                point_list.push_back(point);
                if ((index < (number_of_patches)) || (index > (number_of_patches + number_of_points - 1))) {
                    point.BC = 0;
                }
                else {
                    point.BC = 1;
                }
            }

        break;
        case 2:
            std::cout << "2d implementation trial" << std::endl;
            for (int i = 0; i < total_points; i++) {
                for (int j = 0; j < total_points; j++) {
                    point.Nr = index;
                    point.X = { Delta/2 + j*Delta, Delta/2 + i * Delta , 0};
                    point.x = point.X;
                    index += 1;
                    point_list.push_back(point);
                }
            }
            for (auto& i : point_list) {
                for(int j = 0; j < (total_points*total_points); j++){
                    if (i.Nr > ((j * total_points) + (number_of_patches - 1)) || (i.Nr < (number_of_patches + number_of_points + (j * total_points)))) {
                        point.BC = 0;
                        //std::cout << " , BC : " << i.BC <<std::endl;
                    }
                    else {
                        point.BC = 1;
                        //std::cout << " , BC : " << i.BC <<std::endl;
                    }
                }
            }

        break;
        case 3:
            std::cout << "3d implementation trial" << std::endl;
            for (int i = 0; i < total_points; i++) {
                for (int j = 0; j < total_points; j++) {
                    for (int k = 0; k < total_points; k++) {
                        point.Nr = index;
                        point.X = { Delta/2 + k * Delta, Delta/2 + j * Delta , Delta/2 + i * Delta};
                        point.x = point.X;
                        index += 1;
                        point_list.push_back(point);
                    }
                }
            }
        break;

        default:

        break;

    }

    // Debugging

    for (auto & i : point_list) {
        std::cout << "Nr: " << i.Nr << ", X: [";

        // Print X vector
        for (const auto& val : i.X) {
            std::cout << val << " ";
        }
        std::cout << "], x: [";

        // Print x vector
        for (const auto& val : i.x) {
            std::cout << val << " ";
        }
        std::cout << "], Volume: " << i.volume << std::endl;

    }
    return point_list;
}



#endif //MESH_H
