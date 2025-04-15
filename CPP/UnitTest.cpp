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
    double expected_energy = 0.5 * C1 * L * std::pow((l / L - 1), 2) ;
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

    // Print neighbor lists for debugging
    //std::cout << "DEBUG: Point 1 neighbors: ";
    for (const auto& neighbor : p1.neighbour_list_1N) {
        //std::cout << neighbor[0] << " ";
    }
    std::cout << std::endl;

    //std::cout << "DEBUG: Point 2 neighbors: ";
    for (const auto& neighbor : p2.neighbour_list_1N) {
        //std::cout << neighbor[0] << " ";
    }
    std::cout << std::endl;

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

    Eigen::Vector3d xi_1 = p2.x - p1.x;
    Eigen::Vector3d expected_residual = -xi_1 * C1 * V_eff * ((1 / L) - (1 / l));

    // Check residual calculation - try both points
    bool point1_passed = (points[0].Ra_1 - expected_residual).norm() < 1e-6;
    bool point2_passed = (points[1].Ra_1 - (-expected_residual)).norm() < 1e-6;

    if (point1_passed) {
        std::cout << "Residual Calculation Test for Point 1 PASSED!\n";
    } else if (point2_passed) {
        std::cout << "Residual Calculation Test for Point 2 (with negative expected) PASSED!\n";
    } else {
        std::cout << "Residual Calculation Test FAILED for both points!\n";
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

void testTwoNeighborInteraction() {
    std::cout << "\n=== Testing Two-Neighbor Interaction ===\n";

    // Create points with proper initialization
    Points p1, p2, p3;

    // Set coordinates
    p1.X = Eigen::Vector3d(0.0, 0.0, 0.0);
    p1.x = Eigen::Vector3d(0.0, 0.0, 0.0);
    p2.X = Eigen::Vector3d(1.0, 0.0, 0.0);
    p2.x = Eigen::Vector3d(2.0, 0.0, 0.0);
    p3.X = Eigen::Vector3d(0.0, 1.0, 0.0);
    p3.x = Eigen::Vector3d(0.0, 2.0, 0.0);

    // Set neighbor lists
    p1.neighbour_list_2N.push_back({1, 2});
    p2.neighbour_list_2N.push_back({0, 2});
    p3.neighbour_list_2N.push_back({0, 1});

    // Set initial counts
    p1.n2 = 1;
    p2.n2 = 1;
    p3.n2 = 1;

    // Constants
    double C2 = 1.0;
    double delta = 2.0;
    int PD = 2;

    // Create vector of points
    std::vector<Points> points = {p1, p2, p3};

    // Calculate energy, residuals, and stiffness
    calculate_r(PD, points, 0.0, C2, 0.0, delta);

    // Manual calculation for comparison
    Eigen::Vector3d Xi_1 = p2.X - p1.X;
    Eigen::Vector3d Xi_2 = p3.X - p1.X;
    Eigen::Vector3d xi_1 = p2.x - p1.x;
    Eigen::Vector3d xi_2 = p3.x - p1.x;

    Eigen::Vector3d Xi = Xi_1.cross(Xi_2);
    Eigen::Vector3d xi = xi_1.cross(xi_2);

    double A = Xi.norm();
    double a = xi.norm();

    double expected_energy = 0.5 * C2 * A * std::pow((a / A - 1), 2);
    std::cout << "A = " << A << ", a = " << a << std::endl;

    // Check energy calculation
    if (std::abs(points[0].psi - expected_energy) < 1e-6) {
        std::cout << "Two-Neighbor Energy Calculation Test PASSED!\n";
    } else {
        std::cout << "Two-Neighbor Energy Calculation Test FAILED!\n";
        std::cout << "Expected Energy: " << expected_energy << "\n";
        std::cout << "Calculated Energy: " << points[0].psi << "\n";
    }
}

void testThreeNeighborInteraction() {
    std::cout << "\n=== Testing Three-Neighbor Interaction ===\n";

    // Create points with proper initialization
    Points p1, p2, p3, p4;

    // Set coordinates
    p1.X = Eigen::Vector3d(0.0, 0.0, 0.0);
    p1.x = Eigen::Vector3d(0.0, 0.0, 0.0);
    p2.X = Eigen::Vector3d(1.0, 0.0, 0.0);
    p2.x = Eigen::Vector3d(2.0, 0.0, 0.0);
    p3.X = Eigen::Vector3d(0.0, 1.0, 0.0);
    p3.x = Eigen::Vector3d(0.0, 2.0, 0.0);
    p4.X = Eigen::Vector3d(0.0, 0.0, 1.0);
    p4.x = Eigen::Vector3d(0.0, 0.0, 2.0);

    // Set neighbor lists
    p1.neighbour_list_3N.push_back({1, 2, 3});
    p2.neighbour_list_3N.push_back({0, 2, 3});
    p3.neighbour_list_3N.push_back({0, 1, 3});
    p4.neighbour_list_3N.push_back({0, 1, 2});

    // Set initial counts
    p1.n3 = 1;
    p2.n3 = 1;
    p3.n3 = 1;
    p4.n3 = 1;

    // Constants
    double C3 = 1.0;
    double delta = 2.0;
    int PD = 3;

    // Create vector of points
    std::vector<Points> points = {p1, p2, p3, p4};

    // Calculate energy, residuals, and stiffness
    calculate_r(PD, points, 0.0, 0.0, C3, delta);

    // Manual calculation for comparison
    Eigen::Vector3d Xi_1 = p2.X - p1.X;
    Eigen::Vector3d Xi_2 = p3.X - p1.X;
    Eigen::Vector3d Xi_3 = p4.X - p1.X;
    Eigen::Vector3d xi_1 = p2.x - p1.x;
    Eigen::Vector3d xi_2 = p3.x - p1.x;
    Eigen::Vector3d xi_3 = p4.x - p1.x;

    double V = Xi_1.dot(Xi_2.cross(Xi_3));
    double v = xi_1.dot(xi_2.cross(xi_3));

    double expected_energy = 0.5 * C3 * V * std::pow((v / V - 1), 2);
    std::cout << "V = " << V << ", v = " << v << std::endl;

    // Check energy calculation
    if (std::abs(points[0].psi - expected_energy) < 1e-6) {
        std::cout << "Three-Neighbor Energy Calculation Test PASSED!\n";
    } else {
        std::cout << "Three-Neighbor Energy Calculation Test FAILED!\n";
        std::cout << "Expected Energy: " << expected_energy << "\n";
        std::cout << "Calculated Energy: " << points[0].psi << "\n";
    }
}

void runUnitTests(int PD) {
    testEnergyCalculation();
    testResidualCalculation();
    testStiffnessCalculation();
    if(PD == 2)
        testTwoNeighborInteraction();
    if(PD == 3)
        testThreeNeighborInteraction();
}