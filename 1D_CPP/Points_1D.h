#ifndef POINTS_H
#define POINTS_H

#include <iostream>
#include <vector>
#include <string>

class Points {
public:
    int Nr{};                         // Point index
    double X{};                       // Material coordinate (1D)
    double x{};                       // Spatial coordinate (1D)
    double volume{};

    // 1D neighbor list
    std::vector<std::vector<int>> neighbour_list_1N;
    std::vector<std::vector<int>> neighbour_list;

    // Residuals and stiffness values (1D scalars)
    double Ra_1 = 0.0;
    double Ra_sum = 0.0;
    double Kab_1 = 0.0;
    double Kab_sum = 0.0;

    std::string Flag = "Patch";

    int DOF = 0;                      // Degree of freedom (1D)
    std::vector<double> DOC;         // Degree of constraint (if used)
    int BC = 0;                       // Boundary condition (1D)
    double BCval = 0.0;              // BC value (1D)
    size_t n1 = 0;                   // Number of 1-neighbors

    double V_eff = 0.0;              // Effective volume
    double psi = 0.0;                // Energy density

    // Constructor
    Points(int Nr, double X, double x, double volume)
        : Nr(Nr), X(X), x(x), volume(volume) {
        Ra_1 = Ra_sum = 0.0;
        Kab_1 = Kab_sum = 0.0;
        DOF = 0;
        BC = 0;
        BCval = 0.0;
        V_eff = psi = 0.0;
    }

    // Default constructor
    Points() = default;
};

#endif // POINTS_H
