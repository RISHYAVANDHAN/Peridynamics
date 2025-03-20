#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <vector>
#include <set>
#include <cmath>
#include "Points.h"

inline bool NeighbourCheck(const Points& a, const Points& b, double delta) {
    Eigen::Vector3d diff = a.X - b.X;
    return diff.norm() < delta;
}

inline void generate_neighbour_list(int PD, std::vector<Points>& point_list, double delta) {
    if (point_list.empty()) return;

    for (size_t i = 0; i < point_list.size(); ++i) {
        // Clear existing neighbor lists
        point_list[i].neighbour_list_1N.clear();
        point_list[i].neighbour_list_2N.clear();
        point_list[i].neighbour_list_3N.clear();
        point_list[i].n1 = point_list[i].n2 = point_list[i].n3 = 0;

        // 1-Neighbors
        for (size_t j = 0; j < point_list.size(); ++j) {
            if (i != j && NeighbourCheck(point_list[i], point_list[j], delta)) {
                point_list[i].neighbour_list_1N.push_back({static_cast<int>(j)});
                point_list[i].n1++;
            }
        }

        // 2-Neighbors (non-collinear pairs)
        if (PD >= 2) {
            for (size_t j = 0; j < point_list.size(); ++j) {
                for (size_t k = j + 1; k < point_list.size(); ++k) {
                    if (i != j && i != k && j != k &&
                        NeighbourCheck(point_list[i], point_list[j], delta) &&
                        NeighbourCheck(point_list[i], point_list[k], delta) &&
                        NeighbourCheck(point_list[j], point_list[k], delta)) {

                        Eigen::Vector3d vec1 = point_list[j].X - point_list[i].X;
                        Eigen::Vector3d vec2 = point_list[k].X - point_list[i].X;
                        if (vec1.cross(vec2).norm() > 1e-12) {  // Non-collinear
                            point_list[i].neighbour_list_2N.push_back({static_cast<int>(j), static_cast<int>(k)});
                            point_list[i].n2++;
                        }
                    }
                }
            }
        }

        // 3-Neighbors (non-coplanar triplets)
        if (PD == 3) {
            for (size_t j = 0; j < point_list.size(); ++j) {
                for (size_t k = j + 1; k < point_list.size(); ++k) {
                    for (size_t l = k + 1; l < point_list.size(); ++l) {
                        if (i != j && i != k && i != l && j != k && j != l && k != l &&
                            NeighbourCheck(point_list[i], point_list[j], delta) &&
                            NeighbourCheck(point_list[i], point_list[k], delta) &&
                            NeighbourCheck(point_list[i], point_list[l], delta) &&
                            NeighbourCheck(point_list[j], point_list[k], delta) &&
                            NeighbourCheck(point_list[j], point_list[l], delta) &&
                            NeighbourCheck(point_list[k], point_list[l], delta)) {

                            Eigen::Vector3d vec1 = point_list[j].X - point_list[i].X;
                            Eigen::Vector3d vec2 = point_list[k].X - point_list[i].X;
                            Eigen::Vector3d vec3 = point_list[l].X - point_list[i].X;
                            double vol = std::abs(vec1.dot(vec2.cross(vec3)));
                            if (vol > 1e-12) {  // Non-coplanar
                                point_list[i].neighbour_list_3N.push_back({static_cast<int>(j), static_cast<int>(k), static_cast<int>(l)});
                                point_list[i].n3++;
                            }
                        }
                    }
                }
            }
        }

        // Debug: Print neighbor counts
        std::cout << "Point " << i << " 1N neighbors: " << point_list[i].n1 << "\n";
        std::cout << "Point " << i << " 2N neighbors: " << point_list[i].n2 << "\n";
        std::cout << "Point " << i << " 3N neighbors: " << point_list[i].n3 << "\n";
        point_list[i].neighbour_list = point_list[i].neighbour_list_1N;
        point_list[i].neighbour_list.insert(point_list[i].neighbour_list.end(),
                                            point_list[i].neighbour_list_2N.begin(),
                                            point_list[i].neighbour_list_2N.end());
        point_list[i].neighbour_list.insert(point_list[i].neighbour_list.end(),
                                            point_list[i].neighbour_list_3N.begin(),
                                            point_list[i].neighbour_list_3N.end());
    }
}


#endif // NEIGHBOUR_H