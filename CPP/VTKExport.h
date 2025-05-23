#ifndef VTKEXPORT_H
#define VTKEXPORT_H

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include "Points.h"

void write_vtk(const std::vector<Points>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK file for peridynamics mesh visualization\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Write point coordinates - horizontal 1D layout
    file << "POINTS " << points.size() << " double\n";
    for (const auto& point : points) {
        // Using X[0] for proper original spacing, y=0, z=0 for clean 1D view
        file << point.X[0] << " 0.0 0.0\n";
    }

    // Write cell data (each point is a vertex cell)
    file << "CELLS " << points.size() << " " << 2 * points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1 " << i << "\n";
    }

    // Cell types (vertex = 1)
    file << "CELL_TYPES " << points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1\n";
    }

    // Start point data section
    file << "POINT_DATA " << points.size() << "\n";

    // Output point type (0=patch, 1=point)
    file << "SCALARS point_type int 1\n";
    file << "LOOKUP_TABLE point_type_table\n";
    for (const auto& point : points) {
        file << (point.Flag == "Point" ? 1 : 0) << "\n";
    }

    // Color lookup table (RED patches, BLUE points)
    file << "LOOKUP_TABLE point_type_table 2\n";
    file << "1.0 0.0 0.0 1.0\n";  // RED for patches (0)
    file << "0.0 0.0 1.0 1.0\n";  // BLUE for points (1)

    // Output volume
    file << "SCALARS volume double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        file << point.volume << "\n";
    }

    // Output displacement magnitude
    file << "SCALARS displacement_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        double dx = point.x[0] - point.X[0];
        double dy = point.x[1] - point.X[1];
        double dz = point.x[2] - point.X[2];
        file << sqrt(dx*dx + dy*dy + dz*dz) << "\n";
    }

    // Output boundary conditions
    file << "VECTORS boundary_condition double\n";
    for (const auto& point : points) {
        file << point.BC[0] << " " << point.BC[1] << " " << point.BC[2] << "\n";
    }

    // Output strain energy
    file << "SCALARS strain_energy double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        file << point.psi << "\n";
    }

    // Output displacements
    file << "VECTORS displacement double\n";
    for (const auto& point : points) {
        file << (point.x[0]-point.X[0]) << " "
             << (point.x[1]-point.X[1]) << " "
             << (point.x[2]-point.X[2]) << "\n";
    }

    file.close();
    std::cout << "VTK file written to: " << filename << std::endl;
}

void write_residual_vtk(const std::vector<Points>& points, const std::string& filename, int iteration) {
    std::string iter_filename = filename + "_iter_" + std::to_string(iteration) + ".vtk";
    std::ofstream file(iter_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << iter_filename << " for writing.\n";
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK file for residual forces in Newton-Raphson iteration " << iteration << "\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Horizontal 1D point layout
    file << "POINTS " << points.size() << " double\n";
    for (const auto& point : points) {
        file << point.X[0] << " 0.0 0.0\n";
    }

    // Cell data
    file << "CELLS " << points.size() << " " << 2 * points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1 " << i << "\n";
    }

    file << "CELL_TYPES " << points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1\n";
    }

    file << "POINT_DATA " << points.size() << "\n";

    // Point type coloring
    file << "SCALARS point_type int 1\n";
    file << "LOOKUP_TABLE point_type_table\n";
    for (const auto& point : points) {
        file << (point.Flag == "Point" ? 1 : 0) << "\n";
    }

    file << "LOOKUP_TABLE point_type_table 2\n";
    file << "1.0 0.0 0.0 1.0\n";  // RED patches
    file << "0.0 0.0 1.0 1.0\n";  // BLUE points

    // Residual magnitude
    file << "SCALARS residual_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        file << point.Ra_sum.norm() << "\n";
    }

    // Residual vectors
    file << "VECTORS residual_force double\n";
    for (const auto& point : points) {
        file << point.Ra_sum[0] << " " << point.Ra_sum[1] << " " << point.Ra_sum[2] << "\n";
    }

    file.close();
    std::cout << "Residual VTK file written to: " << iter_filename << std::endl;
}

#endif // VTKEXPORT_H