#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include "Points.h"

// Debug verbosity levels
enum class DebugLevel {
    NONE = 0,   // No debugging
    LOW = 1,    // Minimal output
    MEDIUM = 2, // Moderate detail
    HIGH = 3,   // Comprehensive output
    EXTREME = 4 // Extremely detailed
};

class DebugUtils {
public:
    // Declare the static member with extern
    static DebugLevel currentLevel;

    // Set debug level globally
    static void setDebugLevel(DebugLevel level) {
        currentLevel = level;
    }

    // Conditional debug print for force calculations
    static void printForceDebug(int pointNr, double L, double l, double stretch,
                                const Eigen::Vector3d& Ra_1,
                                const Eigen::Vector3d& Ra_2,
                                const Eigen::Vector3d& Ra_3) {
        if (currentLevel >= DebugLevel::MEDIUM) {
            std::cout << "Force Debug [Point " << pointNr << "]:\n"
                      << "  Reference Length (L): " << L << "\n"
                      << "  Current Length (l): " << l << "\n"
                      << "  Stretch: " << stretch << "\n"
                      << "  Ra_1: " << Ra_1.transpose() << "\n"
                      << "  Ra_2: " << Ra_2.transpose() << "\n"
                      << "  Ra_3: " << Ra_3.transpose() << "\n";
        }
    }

    // Detailed neighbor list inspection
    static void printNeighborDetails(const Points& point, int PD) {
        if (currentLevel >= DebugLevel::MEDIUM) {
            std::cout << "Neighbor Details [Point " << point.Nr << "]:\n";

            std::cout << "  1-Neighbors (" << point.n1 << "): ";
            for (const auto& n : point.neighbour_list_1N) {
                std::cout << "{ ";
                for (int val : n) std::cout << val << " ";
                std::cout << "} ";
            }
            std::cout << "\n";

            if (PD >= 2) {
                std::cout << "  2-Neighbors (" << point.n2 << "): ";
                for (const auto& n : point.neighbour_list_2N) {
                    std::cout << "{ ";
                    for (int val : n) std::cout << val << " ";
                    std::cout << "} ";
                }
                std::cout << "\n";
            }

            if (PD == 3) {
                std::cout << "  3-Neighbors (" << point.n3 << "): ";
                for (const auto& n : point.neighbour_list_3N) {
                    std::cout << "{ ";
                    for (int val : n) std::cout << val << " ";
                    std::cout << "} ";
                }
                std::cout << "\n";
            }
        }
    }

    // Boundary condition inspection
    static void printBoundaryConditions(const Points& point) {
        if (currentLevel >= DebugLevel::LOW) {
            std::cout << "Boundary Conditions [Point " << point.Nr << "]:\n"
                      << "  Flag: " << point.Flag << "\n"
                      << "  BC: " << point.BC.transpose() << "\n"
                      << "  DOF: " << point.DOF.transpose() << "\n"
                      << "  DOC: " << point.DOC.transpose() << "\n";
        }
    }

    // Configuration comparison
    static void compareConfigurations(const Points& point) {
        if (currentLevel >= DebugLevel::LOW) {
            std::cout << "Configuration Comparison [Point " << point.Nr << "]:\n"
                      << "  Reference (X): " << point.X.transpose() << "\n"
                      << "  Current (x):   " << point.x.transpose() << "\n"
                      << "  Difference:    " << (point.x - point.X).transpose() << "\n";
        }
    }

    // Energy and Volume Diagnostics
    static void printEnergyDiagnostics(const Points& point) {
        if (currentLevel >= DebugLevel::MEDIUM) {
            std::cout << "Energy Diagnostics [Point " << point.Nr << "]:\n"
                      << "  Effective Volume: " << point.V_eff << "\n"
                      << "  Energy Density (psi): " << point.psi << "\n";
        }
    }
};

// Define the static member with inline to avoid multiple definition errors
inline DebugLevel DebugUtils::currentLevel = DebugLevel::NONE;

#endif // DEBUG_UTILS_H
