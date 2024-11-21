#pragma once
#include "Grid.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
// #include <mkl.h>
#include <vector>
#include <array>
#include <numeric>
void rearrange_K_global(std::vector<std::vector<std::array<double, 4>>> &K_glob, std::vector<std::vector<std::array<double, 4>>> &K_glob_rearanged, const std::vector<int> &computeNodesIdx, std::vector<Node> &Nodes)
{
    #pragma omp parallel for
    for (int node = 0; node < computeNodesIdx.size(); ++node)
    {
        int node_compute_idx = node;
        int node_global_idx = computeNodesIdx[node];
        bool is_D_written = false;
        for (int j = 0; j < Nodes[node_global_idx].neighbours1.size(); ++j)
        {
            int neighb_global_idx = Nodes[node_global_idx].neighbours1[j];
            int n1_compute_idx = Nodes[neighb_global_idx].global_idx;

            if (n1_compute_idx != -1)
            {
                if (n1_compute_idx < node_compute_idx)
                {
                    K_glob_rearanged[node].push_back(K_glob[node][j + 1]);
                }
                else
                {
                    if (!is_D_written)
                    {
                        K_glob_rearanged[node].push_back(K_glob[node][0]);
                        is_D_written = true;
                        K_glob_rearanged[node].push_back(K_glob[node][j + 1]);
                    }
                    else
                    {
                        K_glob_rearanged[node].push_back(K_glob[node][j + 1]);
                    }
                }
            }
        }
        if (!is_D_written)
        {
            K_glob_rearanged[node].push_back(K_glob[node][0]);
        }
    }
    K_glob.clear();
    K_glob.shrink_to_fit();
}

void vec_to_eigen(const double *vec, EMx &eig, int N){
    for(int i=0; i<N; ++i) eig(i) = vec[i];
}

void eigen_to_vec(const EMx &eig, double *vec, int N){
    for(int i=0; i<N; ++i) vec[i] = eig(i);
}

void crs_to_eigen(int *ia_arr, int *ja_arr, std::vector<Triplet> &tripletList,
                Eigen::SparseMatrix<double, Eigen::RowMajor> &Ktot, int DOF){
    size_t trc=0;
    for(int row=0; row<DOF; ++row){
        int colbeg = ia_arr[row];
        int colend = ia_arr[row+1];
        for(int ja_idx=colbeg; ja_idx<colend; ++ja_idx){
            int col = ja_arr[ja_idx];
            tripletList[trc++] = Triplet(row, col, 1);
        }
    }
    Ktot.setFromTriplets(tripletList.begin(), tripletList.begin() + trc);
    Ktot.makeCompressed();
}

void Solve_System(const std::vector<std::vector<std::array<double, 4>>> &data, std::vector<std::array<double, 2>> residual,
                    EMx &dx, EMx &Rtot, 
                    Eigen::SparseMatrix<double, Eigen::RowMajor> &Ktot, 
                    Eigen::BiCGSTAB< Eigen::SparseMatrix<double, Eigen::RowMajor>> &solver)
{
    int k = 0;
    double* values = Ktot.valuePtr();
    for (int i = 0; i < data.size(); ++i)
    {
        for (int j = 0; j < data[i].size(); ++j)
        {
            values[k++] = -data[i][j][0];
            values[k++] = -data[i][j][1];
        }
        for (int j = 0; j < data[i].size(); ++j)
        {
            values[k++] = -data[i][j][2];
            values[k++] = -data[i][j][3];
        }
    }

    int r = 0;
    for (int i = 0; i < residual.size(); ++i)
    {
        Rtot(r++) = residual[i][0];
        Rtot(r++) = residual[i][1];
    }
    

    solver.compute(Ktot);
    dx = solver.solve(Rtot);

}