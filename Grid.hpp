#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <vector>
#include <cstring>
#include <algorithm>
#include <array>
#include <unordered_map>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include <filesystem>
using Triplet = Eigen::Triplet<double>;
using EMx = Eigen::MatrixXd;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////             Node               ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Node
{
    double coordinate[2];
    double X[2];
    double Vol;
    int Flag;
    int global_idx;
    std::vector<int> neighbours1;
    std::vector<std::array<int, 2>> neighbours2; // index of neighb2 {24,12}
    Node(double x, double y, int flag, double vol) : Flag(flag), Vol(vol)
    {
        coordinate[0] = x;
        coordinate[1] = y;
        X[0] = x;
        X[1] = y;
        global_idx = -1;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////             Grid               ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Grid
{
    int max_Num_Neighbors = 0; // max node in horizon (max neighb + self node(1))
    int Nx_tot;                // number of nodes in x + patch nodes
    int Ny_tot;
    int Empty_Nodes_num;
    int Patch_num; // number of patch idx (basing on horizon)
    double L;      // assumed equidistance grids
    double Corners[4][2];
    double extended_Corners[4][2];
    double Horizon;
    double V_Base;
    std::string Patch_mode;
    std::vector<Node> Nodes;
    std::vector<int> Compute_Nodes_idx;
    // constructor
    Grid(double(Corners_i)[4][2], double l, double horizon, std::string &patch_mode) : L(l), Patch_mode(patch_mode), Horizon(horizon)
    {
        V_Base = horizon * horizon * 3.14159265;
        Patch_num = floor(horizon / l); //
        memcpy(Corners, Corners_i, sizeof(Corners));
        if (patch_mode == "vertical_patch") // we extend the whole domain but only create valid nodes (flag version)
        {
            extended_Corners[0][0] = Corners[0][0] - floor(horizon / l) * l;
            extended_Corners[1][0] = Corners[1][0] + floor(horizon / l) * l;
            extended_Corners[2][0] = Corners[2][0] + floor(horizon / l) * l;
            extended_Corners[3][0] = Corners[3][0] - floor(horizon / l) * l;
            extended_Corners[0][1] = Corners[0][1] - floor(horizon / l) * l;
            extended_Corners[1][1] = Corners[1][1] - floor(horizon / l) * l;
            extended_Corners[2][1] = Corners[2][1] + floor(horizon / l) * l;
            extended_Corners[3][1] = Corners[3][1] + floor(horizon / l) * l;

            // extended_Corners[0][0] = Corners[0][0] - floor(horizon / l) * l;
            // extended_Corners[1][0] = Corners[1][0] + floor(horizon / l) * l;
            // extended_Corners[2][0] = Corners[2][0] + floor(horizon / l) * l;
            // extended_Corners[3][0] = Corners[3][0] - floor(horizon / l) * l;
            // extended_Corners[0][1] = Corners[0][1];
            // extended_Corners[1][1] = Corners[1][1];
            // extended_Corners[2][1] = Corners[2][1];
            // extended_Corners[3][1] = Corners[3][1];

        }

        Nx_tot = (extended_Corners[1][0] - extended_Corners[0][0]) / L + 1; /// number of nodes + patch (valid for inner nodes)
        Ny_tot = (extended_Corners[3][1] - extended_Corners[0][1]) / L + 1;
        std::cout << "Patch mode = " << patch_mode << std::endl;
        std::cout << "Nx_tot = " << Nx_tot << std::endl;
        std::cout << "Ny_tot = " << Ny_tot << std::endl;
        Nodes.reserve(Nx_tot * Ny_tot); // very important
        for (size_t j = 0; j < Ny_tot; ++j)
        {
            double alpha = 1;
            for (size_t i = 0; i < Nx_tot; ++i)
            {
                double x = extended_Corners[0][0] + (i * L);
                double y = extended_Corners[0][1] + (j * L);
                Nodes.emplace_back(x, y, 0, alpha * L * L);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////



    /////////////////////////////////////////////////////////////////////////////
    /////////////////////       Access operators      ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    Node &operator()(size_t i, size_t j)
    {
        return Nodes[Nx_tot * j + i];
    }
    Node operator()(size_t i, size_t j) const { return Nodes[Nx_tot * j + i]; }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    ///////////////////       Finding the distance      /////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    double p_distance2(Node &a, Node &b)
    {
        double diff_x = a.coordinate[0] - b.coordinate[0];
        double diff_y = a.coordinate[1] - b.coordinate[1];
        double distance2 = diff_x * diff_x + diff_y * diff_y;
        return distance2;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    ///////////////////          Post Process           /////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    void writeVTK(int id)
    {
        // Vtk output
        // std::filesystem::create_directory("vtk");
        std::ofstream vtkFile("vtk/paraview_" + std::to_string(id) + ".vtk");

        vtkFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkFile << "Node coordinates with flag values" << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        vtkFile << "POINTS " << (Nodes.size() - Empty_Nodes_num) << " double" << std::endl;

        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.coordinate[0] << " " << node.coordinate[1] << " 0.0" << std::endl;
            }
        }

        vtkFile << std::endl;
        vtkFile << "POINT_DATA " << (Nodes.size() - Empty_Nodes_num) << std::endl;
        vtkFile << "SCALARS flag int" << std::endl; // change data type base of debugg value
        vtkFile << "LOOKUP_TABLE default" << std::endl;

        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.Flag << std::endl; // change for debugg
            }
        }

        vtkFile << "SCALARS G_idx int" << std::endl; // change data type base of debugg value        
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.global_idx << std::endl; // change for debugg
            }
        }

        vtkFile << "SCALARS NNbr_1 int" << std::endl; // change data type base of debugg value
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.neighbours1.size() << std::endl; // change for debugg
            }
        }

        vtkFile << "SCALARS NNbr_2 int" << std::endl; // change data type base of debugg value
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.neighbours2.size() << std::endl; // change for debugg
            }
        }

        vtkFile << "VECTORS X double" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.X[0] << " " << node.X[1] << " 0.0" << std::endl;
            }
        }

        vtkFile << "VECTORS x double" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.coordinate[0] << " " << node.coordinate[1] << " 0.0" << std::endl;
            }
        }

        vtkFile << "VECTORS disp double" << std::endl;
        for (const auto &node : Nodes)
        {
            if (node.Flag != 3)
            {
                vtkFile << node.coordinate[0]-node.X[0] << " " << node.coordinate[1]-node.X[1] << " 0.0" << std::endl;
            }
        }

        vtkFile.close();
        // std::cout << "VTK file 'grid_2d.vtk' has been written." << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////          Patch nodes and Free Nodes         ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    // patch node = 1;
    // Compute node = 0;

    void Set_Patch_flags_Vertical() // left_patch_flag = 1              right_patch_flag=2
    {

        // #pragma omp parallel for // made slow
        for (int j = 0; j < Ny_tot; ++j)
        {
            for (int i = 0; i < Patch_num; ++i)
            {

                operator()(i, j).Flag = 1;
                // Patch_Nodes_left.push_back(Nx_tot * j + i);
            }
        }
        // #pragma omp parallel for

        for (int j = 0; j < Ny_tot; ++j)
        {
            for (int i = Nx_tot - Patch_num; i < Nx_tot; ++i)
            {

                operator()(i, j).Flag = 2;
                // Patch_Nodes_right.push_back(Nx_tot * j + i);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////                    Volume                   ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void Assign_Vol()
    {
        // This Function assign volume factor Alpha ( 0.5 for compute nodes on axis - 0.25 for nodes on corners)
        for (int i = Patch_num + 1; i < Nx_tot - Patch_num - 1; ++i)
        {
            operator()(i, Patch_num).Vol = 0.5 * L * L;
            operator()(i, Ny_tot - Patch_num - 1).Vol = 0.5 * L * L;
        }
        for (int j = Patch_num + 1; j < Ny_tot - Patch_num - 1; ++j)
        {
            operator()(Patch_num, j).Vol = 0.5 * L * L;
            operator()(Nx_tot - Patch_num - 1, j).Vol = 0.5 * L * L;
        }
        operator()(Patch_num, Patch_num).Vol = 0.5 * L * L;
        operator()(Patch_num, Ny_tot - Patch_num - 1).Vol = 0.5 * L * L;
        operator()(Nx_tot - Patch_num - 1, Ny_tot - Patch_num - 1).Vol = 0.5 * L * L;
        operator()(Nx_tot - Patch_num - 1, Patch_num).Vol = 0.5 * L * L;



        // for (int i = Patch_num ; i < Nx_tot - Patch_num ; ++i)
        // {
        //     operator()(i, Patch_num).Vol = 0.5 * L * L;
        //     operator()(i, Ny_tot - Patch_num ).Vol = 0.5 * L * L;
        // }
        // for (int j = Patch_num ; j < Ny_tot - Patch_num - 1; ++j)
        // {
        //     operator()(Patch_num, j).Vol = 0.5 * L * L;
        //     operator()(Nx_tot - Patch_num , j).Vol = 0.5 * L * L;
        // }
        // operator()(Patch_num, Patch_num).Vol = 0.25 * L * L;
        // operator()(Patch_num, Ny_tot - Patch_num ).Vol = 0.25 * L * L;
        // operator()(Nx_tot - Patch_num , Ny_tot - Patch_num - 1).Vol = 0.25 * L * L;
        // operator()(Nx_tot - Patch_num , Patch_num).Vol = 0.25 * L * L;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////                  1 Neighbor                 ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void find_Neighbors1() // works for every types of patches (base on flag number) find the neghibor1 and calculate X1 simultaneously
    {
        double eps = 1e-8;

#pragma omp parallel for collapse(2) // optimum num threads = 3
        for (int j = 0; j < Ny_tot; ++j)
        {
            for (int i = 0; i < Nx_tot; ++i)
            {
                if (operator()(i, j).Flag == 0)
                {
                    int maxcounter = 1;

                    int start_search_idx_j = std::max(j - Patch_num, 0);
                    int end_search_idx_j = std::min(j + Patch_num + 1, Ny_tot);
                    int start_search_idx_i = std::max(i - Patch_num, 0);
                    int end_search_idx_i = std::min(i + Patch_num + 1, Nx_tot);
                    for (int jj = start_search_idx_j; jj < end_search_idx_j; ++jj)
                    {
                        for (int ii = start_search_idx_i; ii < end_search_idx_i; ++ii)
                        {

                            if ((ii == i && jj == j) || (operator()(ii, jj).Flag == 3))
                            {
                                continue;
                            }

                            double distance = sqrt(p_distance2(operator()(i, j), operator()(ii, jj)));
                            if (distance <= Horizon + eps)
                            {
                                int neighbour_idx = Nx_tot * jj + ii;
                                operator()(i, j).neighbours1.push_back(neighbour_idx);
                                ++maxcounter;
                            }
                        }
                    }
                    if (maxcounter > max_Num_Neighbors)
                    {
                        max_Num_Neighbors = maxcounter;
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////                  2 Neighbor                 ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void find_Neighbors2()
    {
        double eps = 1e-8;
#pragma omp parallel for

        for (int n = 0; n < Nodes.size(); ++n)
        {
            if (Nodes[n].Flag == 0)
            {
                for (size_t i = 0; i < Nodes[n].neighbours1.size(); ++i)
                {
                    for (size_t j = i + 1; j < Nodes[n].neighbours1.size(); ++j)
                    {
                        double diff_x = Nodes[Nodes[n].neighbours1[i]].coordinate[0] - Nodes[Nodes[n].neighbours1[j]].coordinate[0];
                        double diff_y = Nodes[Nodes[n].neighbours1[i]].coordinate[1] - Nodes[Nodes[n].neighbours1[j]].coordinate[1];
                        double distance2 = diff_x * diff_x + diff_y * diff_y;
                        if ((distance2 <= Horizon * Horizon + eps))
                        {
                            int tmp1 = Nodes[n].neighbours1[i];
                            int tmp2 = Nodes[n].neighbours1[j];
                            double X1_x = Nodes[Nodes[n].neighbours1[i]].coordinate[0] - Nodes[n].coordinate[0];
                            double X1_y = Nodes[Nodes[n].neighbours1[i]].coordinate[1] - Nodes[n].coordinate[1];
                            double X2_x = Nodes[Nodes[n].neighbours1[j]].coordinate[0] - Nodes[n].coordinate[0];
                            double X2_y = Nodes[Nodes[n].neighbours1[j]].coordinate[1] - Nodes[n].coordinate[1];
                            double cross_colinear = std::abs(X1_x * X2_y - X2_x * X1_y);
                            if (cross_colinear > eps)
                            {
                                Nodes[n].neighbours2.push_back({tmp1, tmp2});
                            }
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    //////         Flag for patch nodes and free nodes                ///////////
    /////////////////////////////////////////////////////////////////////////////

    void Set_flags_empty() /// Empty nodes =3        vertical patch                         ** Remember : we extend domain for each types of patch but we modify setflag empty for each type
    {
        int counter = 0;
        // #pragma omp parallel
        // #pragma omp for collapse(2)

        for (int j = 0; j < Patch_num; ++j)
            for (int i = Patch_num; i < Nx_tot - Patch_num; ++i)
            {
                {
                    ++counter;
                    operator()(i, j).Flag = 3;
                }
            }
        // #pragma omp parallel for collapse(2)
        // #pragma omp for collapse(2)

        for (int j = Ny_tot - Patch_num; j < Ny_tot; ++j)
        {
            for (int i = Patch_num; i < Nx_tot - Patch_num; ++i)
            {
                {
                    ++counter;

                    operator()(i, j).Flag = 3;
                }
            }
        }
        Empty_Nodes_num = counter;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////                      BC                     ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void Apply_BC_rigid(double load, double step)
    {
        // #pragma omp parallel for collapse(2)

        for (int j = 0; j < Ny_tot; ++j)
        {
            for (int i = 0; i < Patch_num; ++i)
            {
                {
                    operator()(i, j).coordinate[0] -= load + (Patch_num - i) * step;
                }
            }
        }
        int counter = 0;
        // #pragma omp parallel for collapse(2)

        for (int j = 0; j < Ny_tot; ++j)
        {
            for (int i = Nx_tot - Patch_num; i < Nx_tot; ++i)
            {
                {
                    operator()(i, j).coordinate[0] += load + counter * step;
                    ++counter;
                }
            }
            counter = 0;
        }
    }

    void Apply_BC_move_mirror(double load)
    {
        for (int j = 0; j < Ny_tot; ++j){
            operator()(Patch_num, j).coordinate[0] -= load;
            operator()(Nx_tot - Patch_num - 1, j).coordinate[0] += load;
        }
    }

    void Apply_BC_mirror_patch()
    {
        //---------- mid patches
        for (int j = Patch_num; j < Ny_tot - Patch_num; ++j){
            for (int i = 0; i < Patch_num; ++i){
                operator()(i, j).coordinate[0] = operator()(i, j).X[0] - (
                2*(operator()(Patch_num, j).X[0] - operator()(Patch_num, j).coordinate[0]) - 
                (operator()(2*Patch_num-i, j).X[0] - operator()(2*Patch_num-i, j).coordinate[0]));
            }

            for (int i = Nx_tot - Patch_num ; i < Nx_tot; ++i){
                operator()(i, j).coordinate[0] = operator()(i, j).X[0] - (
                2*(operator()(Nx_tot - Patch_num - 1, j).X[0] - operator()(Nx_tot - Patch_num - 1, j).coordinate[0]) - 
                (operator()(2*(Nx_tot - Patch_num - 1) - i, j).X[0] - operator()(2*(Nx_tot - Patch_num - 1) - i, j).coordinate[0]));
            }
        }

        //------------ top bottom patches
        for (int j = 0; j < Patch_num; ++j){
            for (int i = 0; i < Patch_num; ++i){
                operator()(i, j).coordinate[0] = operator()(i, Patch_num).coordinate[0] ;
            }
            for (int i = Nx_tot - Patch_num ; i < Nx_tot; ++i){
                operator()(i, j).coordinate[0] = operator()(i, Patch_num).coordinate[0] ;
            }
        }
        for (int j = Ny_tot - Patch_num; j < Ny_tot; ++j){
            for (int i = 0; i < Patch_num; ++i){
                operator()(i, j).coordinate[0] = operator()(i, Ny_tot - Patch_num - 1).coordinate[0] ;
            }
            for (int i = Nx_tot - Patch_num ; i < Nx_tot; ++i){
                operator()(i, j).coordinate[0] = operator()(i, Ny_tot - Patch_num - 1).coordinate[0] ;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////               Find free nodes               ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void Find_Compute_Nodes()
    {
        int count = 0;

        for (int i = 0; i < Nodes.size(); ++i)
        {
            if (Nodes[i].Flag == 0)
            {
                Compute_Nodes_idx.push_back(i);
                Nodes[i].global_idx = count;
                ++count;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////                 Assign DOF                  ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    void update_compute_nodes(const EMx &dx_arr)
    {
        #pragma omp parallel for
        for (int i = 0; i < Compute_Nodes_idx.size(); ++i)
        {
            Nodes[Compute_Nodes_idx[i]].coordinate[0] += dx_arr(2 * i);
            Nodes[Compute_Nodes_idx[i]].coordinate[1] += dx_arr(2 * i + 1);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////             Solver Utility                  ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    int create_sparse_CSR(std::vector<int> &ia_vec, std::vector<int> &ja_vec)
    {

        int nnzr = 0;
        int ia_counter = 0;
        int ia_counter_glob = 0;
        ia_vec.push_back(0);
        for (int node = 0; node < Compute_Nodes_idx.size(); ++node)
        {
            // std::cout << "compute nodes size = " << Compute_Nodes_idx.size() << std::endl;
            int node_compute_idx = node;
            int node_global_idx = Compute_Nodes_idx[node];
            bool is_D_written = false;
            for (int j = 0; j < Nodes[node_global_idx].neighbours1.size(); ++j)
            {
                int neighb_global_idx = Nodes[node_global_idx].neighbours1[j];
                int n1_compute_idx = Nodes[neighb_global_idx].global_idx;

                if (n1_compute_idx != -1)
                {
                    if (n1_compute_idx < node_compute_idx)
                    {
                        ++nnzr;
                        ++ia_counter;
                        ja_vec.push_back(2 * n1_compute_idx);
                        ja_vec.push_back(2 * n1_compute_idx + 1);
                    }
                    else
                    {
                        if (!is_D_written)
                        {

                            ++nnzr;
                            ++ia_counter;
                            ja_vec.push_back(2 * node);
                            ja_vec.push_back(2 * node + 1);
                            is_D_written = true;
                            ++nnzr;
                            ++ia_counter;
                            ja_vec.push_back(2 * n1_compute_idx);
                            ja_vec.push_back(2 * n1_compute_idx + 1);
                        }
                        else
                        {

                            ++nnzr;
                            ++ia_counter;
                            ja_vec.push_back(2 * n1_compute_idx);
                            ja_vec.push_back(2 * n1_compute_idx + 1);
                        }
                    }
                }
            }
            if (!is_D_written)
            {

                ++nnzr;
                ++ia_counter;
                ja_vec.push_back(2 * node);
                ja_vec.push_back(2 * node + 1);
            }
            ////////////////////////////
            int ja_size = ja_vec.size();
            for (int jj = (ja_size - 2 * ia_counter); jj < ja_size; ++jj)
            {
                // std::cout << "jj = " << std::endl;
                ja_vec.push_back(ja_vec[jj]);
            }
            for (int it = 0; it < 2; ++it)
            {
                ia_counter_glob += 2 * (ia_counter);
                ia_vec.push_back(ia_counter_glob);
            }

            ia_counter = 0; // reset local counter
        }
        return 4 * nnzr;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    int create_sparse_index(std::vector<int> &ia_vec, std::vector<int> &ja_vec)
    {

        int nnzr = 0;
        int ia_counter = 0;
        for (int node = 0; node < Compute_Nodes_idx.size(); ++node)
        {
            // std::cout << "compute nodes size = " << Compute_Nodes_idx.size() << std::endl;
            int node_compute_idx = node;
            int node_global_idx = Compute_Nodes_idx[node];
            bool is_D_written = false;
            for (int j = 0; j < Nodes[node_global_idx].neighbours1.size(); ++j)
            {
                int neighb_global_idx = Nodes[node_global_idx].neighbours1[j];
                int n1_compute_idx = Nodes[neighb_global_idx].global_idx;

                if (n1_compute_idx != -1)
                {
                    if (n1_compute_idx < node_compute_idx)
                    {
                        ++nnzr;
                        ++ia_counter;
                        ia_vec.push_back(2 * node);
                        ja_vec.push_back(2 * n1_compute_idx);
                        ja_vec.push_back(2 * n1_compute_idx + 1);
                    }
                    else
                    {
                        if (!is_D_written)
                        {

                            ++nnzr;
                            ia_vec.push_back(2 * node);
                            ++ia_counter;

                            ja_vec.push_back(2 * node);
                            ja_vec.push_back(2 * node + 1);
                            is_D_written = true;
                            ++nnzr;
                            ia_vec.push_back(2 * node);
                            ++ia_counter;

                            ja_vec.push_back(2 * n1_compute_idx);
                            ja_vec.push_back(2 * n1_compute_idx + 1);
                        }
                        else
                        {

                            ++nnzr;
                            ia_vec.push_back(2 * node);
                            ++ia_counter;

                            ja_vec.push_back(2 * n1_compute_idx);
                            ja_vec.push_back(2 * n1_compute_idx + 1);
                        }
                    }
                }
            }
            if (!is_D_written)
            {

                ++nnzr;
                ia_vec.push_back(2 * node);
                ++ia_counter;

                ja_vec.push_back(2 * node);
                ja_vec.push_back(2 * node + 1);
            }
            ////////////////////////////
            int ja_size = ja_vec.size();
            for (int jj = (ja_size - 2 * ia_counter); jj < ja_size; ++jj)
            {
                // std::cout << "jj = " << std::endl;
                ja_vec.push_back(ja_vec[jj]);
            }
            for (int it = 0; it < ia_counter; ++it)
            {
                ia_vec.push_back(node * 2 + 1);
            }

            ia_counter = 0; // reset local counter
        }
        return 4 * nnzr;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////