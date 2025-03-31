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

    // Write point coordinates (current deformed configuration)
    file << "POINTS " << points.size() << " double\n";
    for (const auto& point : points) {
        file << point.x[0] << " " << point.x[1] << " " << point.x[2] << "\n";
    }

    // Write cell data (each point is a cell)
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

    // Output volume as a scalar
    file << "SCALARS volume double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        file << point.volume << "\n";
    }

    // Output displacement magnitude as a scalar
    file << "SCALARS displacement_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        // Calculate displacement magnitude
        double dx = point.x[0] - point.X[0];
        double dy = point.x[1] - point.X[1];
        double dz = point.x[2] - point.X[2];
        double disp_mag = sqrt(dx*dx + dy*dy + dz*dz);
        file << disp_mag << "\n";
    }

    // Output point type (Patch=0, Point=1) for visualization
    file << "SCALARS point_type int 1\n";
    file << "LOOKUP_TABLE point_type_table\n";
    for (const auto& point : points) {
        int type_value = (point.Flag == "Patch") ? 0 : 1;
        file << type_value << "\n";
    }

    // Create a custom lookup table for point types
    file << "LOOKUP_TABLE point_type_table 2\n";
    file << "1.0 0.0 0.0 1.0\n";  // Red for Patch points
    file << "0.0 0.0 1.0 1.0\n";  // Blue for regular Points

    // Output boundary condition flags as a vector
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

    // If there are displacements, output them as vectors
    file << "VECTORS displacement double\n";
    for (const auto& point : points) {
        double dx = point.x[0] - point.X[0];
        double dy = point.x[1] - point.X[1];
        double dz = point.x[2] - point.X[2];
        file << dx << " " << dy << " " << dz << "\n";
    }

    file.close();
    std::cout << "VTK file written to: " << filename << std::endl;
}

// Optional: Helper function to visualize residual forces during Newton-Raphson iterations
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

    // Write point coordinates (current configuration)
    file << "POINTS " << points.size() << " double\n";
    for (const auto& point : points) {
        file << point.x[0] << " " << point.x[1] << " " << point.x[2] << "\n";
    }

    // Write cell data
    file << "CELLS " << points.size() << " " << 2 * points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1 " << i << "\n";
    }

    // Cell types
    file << "CELL_TYPES " << points.size() << "\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << "1\n";
    }

    // Start point data section
    file << "POINT_DATA " << points.size() << "\n";

    // Output residual force magnitudes as a scalar
    file << "SCALARS residual_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : points) {
        double res_mag = point.Ra_sum.norm();
        file << res_mag << "\n";
    }

    // Output residual forces as vectors
    file << "VECTORS residual_force double\n";
    for (const auto& point : points) {
        file << point.Ra_sum[0] << " " << point.Ra_sum[1] << " " << point.Ra_sum[2] << "\n";
    }

    // Output point type (to maintain color distinction)
    file << "SCALARS point_type int 1\n";
    file << "LOOKUP_TABLE point_type_table\n";
    for (const auto& point : points) {
        int type_value = (point.Flag == "Patch") ? 0 : 1;
        file << type_value << "\n";
    }

    file << "LOOKUP_TABLE point_type_table 2\n";
    file << "1.0 0.0 0.0 1.0\n";  // Red for Patch points
    file << "0.0 0.0 1.0 1.0\n";  // Blue for regular Points

    file.close();
    std::cout << "Residual VTK file written to: " << iter_filename << std::endl;
}

#endif // VTKEXPORT_H