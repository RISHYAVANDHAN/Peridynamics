#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include "Points.h"
#include <vector>
#include <cmath>
#include <algorithm>

inline bool NeighbourCheck(const Points& a, const Points& b, double delta, int PD) {
    double distance_sq = 0.0;
    for (int dim = 0; dim < PD; dim++) {
        double diff = a.X(dim) - b.X(dim);
        distance_sq += diff * diff;
    }
    double distance = std::sqrt(distance_sq);

    // Add these critical checks:
    if (distance > delta * 100) {  // Prevent extreme stretching
        std::cerr << "ERROR: Invalid neighbor distance " << distance
                  << " between points " << a.Nr << " and " << b.Nr
                  << " (delta = " << delta << ")" << std::endl;
        return false;
    }
    if (distance < 1e-12) {
        return false;  // Skip self-interaction
    }
    return (distance < delta);
}

inline void generate_neighbour_list(int PD, std::vector<Points>& points, double delta) {

    for (size_t i = 0; i < points.size(); ++i) {
        auto& point = points[i];
        point.neighbour_list.clear();
        point.neighbour_list_1N.clear();
        point.neighbour_list_2N.clear();
        point.neighbour_list_3N.clear();
        point.n1 = point.n2 = point.n3 = 0;

        // 1-Neighbors
        for (size_t j = 0; j < points.size(); ++j) {
            if (i != j && NeighbourCheck(point, points[j], delta, PD)) {
                point.neighbour_list_1N.push_back({static_cast<int>(j)});
                point.n1++;
            }
        }

        // 2-Neighbors (for 2D and 3D)
        if (PD >= 2) {
            for (size_t j = 0; j < points.size(); ++j) {
                for (size_t k = j + 1; k < points.size(); ++k) {
                    if (i != j && i != k && j != k &&
                        NeighbourCheck(point, points[j], delta, PD) &&
                        NeighbourCheck(point, points[k], delta, PD) &&
                        NeighbourCheck(points[j], points[k], delta, PD)) {

                        Eigen::Vector3d vec1 = points[j].x - point.x;
                        Eigen::Vector3d vec2 = points[k].x - point.x;
                        if (vec1.cross(vec2).norm() > 1e-12) {
                            point.neighbour_list_2N.push_back({static_cast<int>(j), static_cast<int>(k)});
                            point.n2++;
                        }
                    }
                }
            }
        }

        // 3-Neighbors (for 3D only)
        if (PD == 3) {
            for (size_t j = 0; j < points.size(); ++j) {
                for (size_t k = j + 1; k < points.size(); ++k) {
                    for (size_t l = k + 1; l < points.size(); ++l) {
                        if (i != j && i != k && i != l && j != k && j != l && k != l &&
                            NeighbourCheck(point, points[j], delta, PD) &&
                            NeighbourCheck(point, points[k], delta, PD) &&
                            NeighbourCheck(point, points[l], delta, PD) &&
                            NeighbourCheck(points[j], points[k], delta, PD) &&
                            NeighbourCheck(points[j], points[l], delta, PD) &&
                            NeighbourCheck(points[k], points[l], delta, PD)) {

                            Eigen::Vector3d vec1 = points[j].X - point.X;
                            Eigen::Vector3d vec2 = points[k].X - point.X;
                            Eigen::Vector3d vec3 = points[l].X - point.X;
                            if (std::abs(vec1.dot(vec2.cross(vec3))) > 1e-12) {
                                point.neighbour_list_3N.push_back(
                                    {static_cast<int>(j), static_cast<int>(k), static_cast<int>(l)});
                                point.n3++;
                            }
                        }
                    }
                }
            }
        }

        // Consolidate all neighbor lists
        point.neighbour_list = point.neighbour_list_1N;
        if (PD >= 2) {
            point.neighbour_list.insert(point.neighbour_list.end(), point.neighbour_list_2N.begin(),point.neighbour_list_2N.end());
        }
        if (PD == 3) {
            point.neighbour_list.insert(point.neighbour_list.end(), point.neighbour_list_3N.begin(), point.neighbour_list_3N.end());
        }
    }
}

#endif // NEIGHBOUR_H