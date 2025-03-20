#include "UnitTest.h"

void testEnergyCalculation() {
    std::cout << "=== Testing Energy Calculation ===\n";

    // Create points with proper initialization
    Points p1;
    Points p2;

    // Set coordinates
    p1.X = Eigen::Vector3d(0.0, 0.0, 0.0);
    p1.x = Eigen::Vector3d(0.0, 0.0, 0.0);
    p2.X = Eigen::Vector3d(1.0, 0.0, 0.0);
    p2.x = Eigen::Vector3d(2.0, 0.0, 0.0);

    // Set neighbor lists
    p1.neighbour_list_1N.push_back({1});
    p2.neighbour_list_1N.push_back({0});

    // Set initial counts
    p1.n1 = 1;
    p2.n1 = 1;

    // Constants
    double C1 = 1.0, C2 = 1.0, C3 = 1.0;
    double delta = 2.0;
    int PD = 1;

    // Create vector of points
    std::vector<Points> points = {p1, p2};

    // Calculate energy, residuals, and stiffness
    calculate_r(PD, points, C1, C2, C3, delta);

    // Manual calculation for comparison
    double L = (p2.X - p1.X).norm();
    double l = (p2.x - p1.x).norm();
    double expected_energy = 0.5 * C1 * L * std::pow((l / L - 1), 2);
    std::cout << "L = " << L << ", l = " << l << std::endl;

    // Check energy calculation
    if (std::abs(points[0].psi - expected_energy) < 1e-6) {
        std::cout << "Energy Calculation Test PASSED!\n";
    } else {
        std::cout << "Energy Calculation Test FAILED!\n";
        std::cout << "Expected Energy: " << expected_energy << "\n";
        std::cout << "Calculated Energy: " << points[0].psi << "\n";
    }
}

void testResidualCalculation() {
    std::cout << "\n=== Testing Residual Calculation ===\n";

    // Create points with proper initialization
    Points p1;
    Points p2;

    // Set coordinates
    p1.X = Eigen::Vector3d(0.0, 0.0, 0.0);
    p1.x = Eigen::Vector3d(0.0, 0.0, 0.0);
    p2.X = Eigen::Vector3d(1.0, 0.0, 0.0);
    p2.x = Eigen::Vector3d(2.0, 0.0, 0.0);

    // Set neighbor lists
    p1.neighbour_list_1N.push_back({1});
    p2.neighbour_list_1N.push_back({0});

    // Set initial counts
    p1.n1 = 1;
    p2.n1 = 1;

    // Constants
    double C1 = 1.0, C2 = 1.0, C3 = 1.0;
    double delta = 2.0;
    int PD = 1;

    // Create vector of points
    std::vector<Points> points = {p1, p2};

    // Calculate energy, residuals, and stiffness
    calculate_r(PD, points, C1, C2, C3, delta);

    // Manual calculation for comparison
    double L = (p2.X - p1.X).norm();
    double l = (p2.x - p1.x).norm();
    const double pi = 3.14159265358979323846;
    double Vh = (4.0 / 3.0) * pi * (delta * delta * delta);
    double V_eff = Vh / points[0].n1;
    std::cout << "V_eff = " << V_eff << std::endl;
    std::cout << "L = " << L << ", l = " << l << std::endl;

    Eigen::Vector3d xi_1 = p2.x - p1.x;
    Eigen::Vector3d expected_residual = xi_1 * C1 * V_eff * ((1 / L) - (1 / l));
    std::cout << "Expected Residual: " << expected_residual.transpose() << "\n";
    std::cout << "Calculated Residual: " << points[0].Ra_1.transpose() << "\n";
    std::cout << "Difference Norm: " << (points[0].Ra_1 - expected_residual).norm() << "\n";

    // Check residual calculation
    if ((points[0].Ra_1 - expected_residual).norm() < 1e-6) {
        std::cout << "Residual Calculation Test PASSED!\n";
    } else {
        std::cout << "Residual Calculation Test FAILED!\n";
        std::cout << "Expected Residual: " << expected_residual.transpose() << "\n";
        std::cout << "Calculated Residual: " << points[0].Ra_1.transpose() << "\n";
        std::cout << "Difference Norm: " << (points[0].Ra_1 - expected_residual).norm() << "\n";
    }
}

void testStiffnessCalculation() {
    std::cout << "\n=== Testing Stiffness Matrix Calculation ===\n";

    // Create points with proper initialization
    Points p1;
    Points p2;

    // Set coordinates
    p1.X = Eigen::Vector3d(0.0, 0.0, 0.0);
    p1.x = Eigen::Vector3d(0.0, 0.0, 0.0);
    p2.X = Eigen::Vector3d(1.0, 0.0, 0.0);
    p2.x = Eigen::Vector3d(2.0, 0.0, 0.0);

    // Set neighbor lists
    p1.neighbour_list_1N.push_back({1});
    p2.neighbour_list_1N.push_back({0});

    // Set initial counts
    p1.n1 = 1;
    p2.n1 = 1;

    // Constants
    double C1 = 1.0, C2 = 1.0, C3 = 1.0;
    double delta = 2.0;
    int PD = 1;

    // Create vector of points
    std::vector<Points> points = {p1, p2};

    // Calculate energy, residuals, and stiffness
    calculate_r(PD, points, C1, C2, C3, delta);

    // Manual calculation for comparison
    double L = (p2.X - p1.X).norm();
    double l = (p2.x - p1.x).norm();
    const double pi = 3.14159265358979323846;
    double Vh = (4.0 / 3.0) * pi * (delta * delta * delta);
    double V_eff = Vh / points[0].n1;
    std::cout << "V_eff = " << V_eff << std::endl;
    std::cout << "L = " << L << ", l = " << l << std::endl;

    Eigen::Vector3d xi_1 = p2.x - p1.x;
    Eigen::Matrix3d expected_stiffness = compute_Kab_1(xi_1, L, l, V_eff, C1);
    std::cout << "Expected Stiffness:\n" << expected_stiffness << "\n";
    std::cout << "Calculated Stiffness:\n" << points[0].Kab_1 << "\n";
    std::cout << "Difference Norm: " << (points[0].Kab_1 - expected_stiffness).norm() << "\n";

    // Check stiffness calculation
    if ((points[0].Kab_1 - expected_stiffness).norm() < 1e-6) {
        std::cout << "Stiffness Calculation Test PASSED!\n";
    } else {
        std::cout << "Stiffness Calculation Test FAILED!\n";
        std::cout << "Expected Stiffness:\n" << expected_stiffness << "\n";
        std::cout << "Calculated Stiffness:\n" << points[0].Kab_1 << "\n";
        std::cout << "Difference Norm: " << (points[0].Kab_1 - expected_stiffness).norm() << "\n";
    }
}

void runUnitTests() {
    testEnergyCalculation();
    testResidualCalculation();
    testStiffnessCalculation();

}
