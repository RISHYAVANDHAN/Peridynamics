
#include "mesh.h"

void VerifyBoundaryConditions(const std::vector<Points>& points, int PD) {
    int constrained_points = 0;
    int free_points = 0;

    for (const auto& p : points) {
        if (p.Flag == "Patch") {
            if (p.BC != Eigen::Vector3i::Zero()) {
                std::cerr << "Invalid BC for Patch point " << p.Nr << "\n";
            }
        } else {
            Eigen::Vector3d BC_double = p.BC.cast<double>();
            if (BC_double.head(PD).minCoeff() <= 0) {
                std::cerr << "Invalid BC for Point " << p.Nr << "\n";
            }
        }
    }

    std::cout << "Boundary Condition Report:\n"
              << "- Constrained points: " << constrained_points << "\n"
              << "- Free points: " << free_points << "\n"
              << "- Total DOFs: " << free_points * PD << "\n";
}

Eigen::Matrix3d Compute_FF(int PD, double d, const std::string& DEFflag) {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d FF = Eigen::Matrix3d::Zero();

    if (DEFflag == "EXT") {
        FF = I;
        FF(0, 0) = 1 + d;
    } else if (DEFflag == "EXP") {
        FF = (1 + d) * I;
    } else if (DEFflag == "SHR") {
        FF = I;
        FF(1, 0) = d;
    }
    return FF;
}

void AssignDOF(std::vector<Points>& point_list, int PD, int& DOFs) {
    DOFs = 0;
    for (auto& point : point_list) {
        for (int dim = 0; dim < PD; dim++) {
            if (point.BC(dim) == 1) {
                DOFs++;
                point.DOF(dim) = DOFs;
            } else {
                point.DOF(dim) = 0;
            }
        }
    }
}

void AssignVolumes(std::vector<Points>& point_list, int PD, double Delta) {
    double volume = Delta;
    if (PD == 2) volume = Delta * Delta;
    else if (PD == 3) volume = Delta * Delta * Delta;

    for (auto& point : point_list) {
        point.volume = volume;
    }
}

std::vector<Points> generate_mesh(int PD, double d, double domain_size, int number_of_points,
                                  int number_of_patches, double Delta, int number_of_right_patches,
                                  const std::string& DEFflag, int& DOFs) {
    double extended_domain_size = domain_size + (number_of_patches + number_of_right_patches) * Delta;
    Eigen::Matrix3d FF = Compute_FF(PD, d, DEFflag);
    std::vector<Points> point_list;
    int index = 0;

    switch (PD) {
        case 1: {
            const int total_points = number_of_patches + number_of_points + number_of_right_patches;
            for (int i = 0; i < total_points; i++) {
                Points point;
                point.Nr = index++;
                point.X = Eigen::Vector3d(Delta / 2 + i * Delta, 0, 0);
                point.x = point.X;

                if (i < number_of_patches || i >= number_of_patches + number_of_points) {
                    point.BC = Eigen::Vector3i(0, 0, 0);
                    point.Flag = "Patch";
                    point.x = point.X;
                } else {
                    point.BC = Eigen::Vector3i(1, 0, 0);
                    point.Flag = "Point";
                    point.BCval = FF * point.X - point.X;
                }
                point_list.push_back(point);
            }
            break;
        }

        case 2: {
            const int total_points = number_of_patches + number_of_points + number_of_right_patches;
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < total_points; j++) {
                    Points point;
                    point.Nr = index++;
                    point.X = Eigen::Vector3d(Delta / 2 + j * Delta, Delta / 2 + i * Delta, 0);
                    point.x = point.X;

                    if (j < number_of_patches || j >= number_of_patches + number_of_points) {
                        point.BC = Eigen::Vector3i(0, 0, 0);
                        point.Flag = "Patch";
                        point.x = point.X;
                    } else {
                        point.BC = Eigen::Vector3i(1, 1, 0);
                        point.Flag = "Point";
                        point.BCval = FF * point.X - point.X;
                    }
                    point_list.push_back(point);
                }
            }
            break;
        }

        case 3: {
            const int total_points = number_of_patches + number_of_points + number_of_right_patches;
            for (int i = 0; i < number_of_points; i++) {
                for (int j = 0; j < number_of_points; j++) {
                    for (int k = 0; k < total_points; k++) {
                        Points point;
                        point.Nr = index++;
                        point.X = Eigen::Vector3d(Delta / 2 + k * Delta, Delta / 2 + j * Delta, Delta / 2 + i * Delta);
                        point.x = point.X;

                        if (k < number_of_patches || k >= number_of_patches + number_of_points) {
                            point.BC = Eigen::Vector3i(0, 0, 0);
                            point.Flag = "Patch";
                            point.x = point.X;
                        } else {
                            point.BC = Eigen::Vector3i(1, 1, 1);
                            point.Flag = "Point";
                            point.BCval = FF * point.X - point.X;
                        }
                        point_list.push_back(point);
                    }
                }
            }
            break;
        }

        default:
            std::cerr << "Invalid PD value. Mesh generation aborted." << std::endl;
            break;
    }

    AssignDOF(point_list, PD, DOFs);
    AssignVolumes(point_list, PD, Delta);
    return point_list;
}
