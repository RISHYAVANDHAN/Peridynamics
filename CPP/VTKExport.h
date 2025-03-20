#ifndef VTKEXPORT_H
#define VTKEXPORT_H

#include <vector>
#include <fstream>
#include "Points.h"

void write_vtk(const std::vector<Points>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK file for point data\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << points.size() << " double\n";
    for (const auto& point : points) {
        file << point.x[0] << " " << point.x[1] << " " << point.x[2] << "\n";
    }

    file << "CELLS " << points.size() << " " << 2 * points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1 " << i << "\n";
    }

    file << "CELL_TYPES " << points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1\n";
    }

    file << "POINT_DATA " << points.size() << "\n";
    file << "SCALARS volume double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        file << point.volume << "\n";
    }

    file.close();
}

#endif // VTKEXPORT_H