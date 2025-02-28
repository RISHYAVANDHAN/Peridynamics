//
// Created by srini on 28/02/2025.
//

#ifndef VTKEXPORT_H
#define VTKEXPORT_H

#include "Points.h"
#include "mesh.h"
#include "Neighbour.h"
#include "compute.h"

void write_vtk(const std::vector<Points>& point_list, const std::string& filename) {
    std::ofstream vtk_file;
    vtk_file.open(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "Failed to open VTK file for writing: " << filename << std::endl;
        return;
    }

    vtk_file << "# vtk DataFile Version 4.2" << std::endl;
    vtk_file << "Generated Mesh Data" << std::endl;
    vtk_file << "ASCII" << std::endl;
    vtk_file << "DATASET POLYDATA" << std::endl;
    vtk_file << "POINTS " << point_list.size() << " float" << std::endl;
    for (const auto& point : point_list) {
        vtk_file << std::fixed << std::setprecision(6);
        vtk_file << point.X[0] << " " << point.X[1] << " " << point.X[2] << std::endl;
    }
    vtk_file << "POINT_DATA " << point_list.size() << std::endl;
    vtk_file << "SCALARS BC int 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& point : point_list) {
        vtk_file << point.BC << std::endl;
    }
    vtk_file << "SCALARS Color float 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& point : point_list) {
        if (point.Flag != "Patch") {
            vtk_file << "1.0\n";  // Points (colored)
        } else {
            vtk_file << "0.0\n";  // Patches (gray)
        }
    }

    vtk_file.close();
    std::cout << "VTK file written to " << filename << std::endl;
}


#endif //VTKEXPORT_H
