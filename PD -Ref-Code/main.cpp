#include <iostream>
#include <iomanip> 
#include <vector>
#include <chrono>
#include <string>
#include "Grid.hpp"
#include "compute.hpp"
#include "solver.hpp"
#include <fstream>
#include <chrono>
#include <memory>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Compute_Corners(const double SiZe, const std::string BCflag, double (&Corners)[4][2])
{
    if (BCflag == "COOK")
    {
        double tmpx = SiZe * 24.0;
        double tmpy = SiZe * 30.0;
        // left south corner
        Corners[0][0] = -tmpx;
        Corners[0][1] = -tmpy;
        // right south corner
        Corners[1][0] = tmpx;
        Corners[1][1] = -tmpy;
        // right north corner
        Corners[2][0] = tmpx;
        Corners[2][1] = tmpy;
        // left north corner
        Corners[3][0] = -tmpx;
        Corners[3][1] = tmpy;
    }
    else
    {
        double tmp = SiZe * 0.5;
        // left south corner
        Corners[0][0] = -tmp;
        Corners[0][1] = -tmp;
        // right south corner
        Corners[1][0] = tmp;
        Corners[1][1] = -tmp;
        // right north corner
        Corners[2][0] = tmp;
        Corners[2][1] = tmp;
        // left north corner
        Corners[3][0] = -tmp;
        Corners[3][1] = tmp;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Compute_PatchFlag(const std::string &BCflag, const std::string &DEFflag, std::string &PatchFlag)
{
    if (BCflag == "DBC")
    {
        PatchFlag = "fullpatch"; // "fullpatch"  or  "horzpatch"   or  "vertpatch"
    }
    else if (BCflag == "STD")
    {
        PatchFlag = "horzpatch"; // "fullpatch"  or  "horzpatch"   or  "vertpatch"
    }
    else if (BCflag == "XTD")
    {
        PatchFlag = "horzpatch"; // "fullpatch"  or  "horzpatch"   or  "vertpatch"
    }
    else if (BCflag == "COOK")
    {
        PatchFlag = "horzpatch"; // "fullpatch"  or  "horzpatch"   or  "vertpatch"
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////      Taking the inputs      ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    if (argc != 11)
    {
        std::cerr << "Usage: " << argv[0] << " SiZe L horizon load loadsteps C1 C2  C3 nn patch_BC_type" << std::endl;
        return 1;
    }
    double SiZe = std::strtod(argv[1], nullptr);
    double L = std::strtod(argv[2], nullptr);
    double Delta_C = std::strtod(argv[3], nullptr);
    double load = std::strtod(argv[4], nullptr);
    int steps = std::strtod(argv[5], nullptr);
    double C1 = std::strtod(argv[6], nullptr);
    double C2 = std::strtod(argv[7], nullptr);
    double C3 = std::strtod(argv[8], nullptr);
    double nn = std::strtod(argv[9], nullptr);
    int patch_BC_type = std::strtod(argv[10], nullptr);    // 0 := Rigid, 1 := Mirror

    double Corners[4][2];
    std::string BCflag = "Normal";
    double horizon = Delta_C * L;
    double load_steps = load / steps;

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////      Grid Specification     ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    std::string patch_mod = "vertical_patch";
    double load_step_patches = 0; // 1; // 0;
    Compute_Corners(SiZe, BCflag, Corners);
    Grid grid(Corners, L, horizon, patch_mod);

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    ///////////////////      Topology Initialization     ////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    double V_Base = grid.V_Base;
    grid.Set_Patch_flags_Vertical();
    grid.Set_flags_empty();
    grid.Find_Compute_Nodes();
    grid.find_Neighbors1();
    grid.find_Neighbors2();
    std::cout << "Number of Compute Nodes = " << grid.Compute_Nodes_idx.size() << std::endl;
    //std::cout << grid.Empty_Nodes_num << std::endl;
    //std::cout << grid.Nodes.size() << std::endl;

    std::cout << " ======== Topology Init Done ! ======== " << std::endl;

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    ///////////////////              Debugger            ////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    bool debug = false;
    int deb_idx = 40;
    double norm2_R = 0;

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////
    ///////////////////         Create sparse info       ////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    int DOF = 2 * grid.Compute_Nodes_idx.size();
    std::vector<int> ia;
    std::vector<int> ja;
    // //////////////////////////////////////////////////////////////////////////
    int nnz = grid.create_sparse_CSR(ia, ja);
    int size_of_ia = ia.size();
    int size_of_ja = ja.size();
    ia.shrink_to_fit(); ja.shrink_to_fit();

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    //////////////              Solver Preparation            ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    
    //// ------------------ dx, Rtot, Ktot, tripletList -------------------------
        Eigen::SparseMatrix<double, Eigen::RowMajor> Ktot(DOF,DOF);
        std::vector<Triplet> tripletList;
        tripletList.reserve(nnz);

        crs_to_eigen(ia.data(), ja.data(), tripletList, Ktot, DOF);

        ia.clear(); ia.shrink_to_fit();
        ja.clear(); ja.shrink_to_fit();
        tripletList.clear(); tripletList.shrink_to_fit();

        EMx dx(DOF,1);
        EMx Rtot(DOF,1);
    //// ------------------ solver ---------------------------
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
        solver.setTolerance(1e-10);
        // if(DOF>4000) solver.setMaxIterations(4000);
        Eigen::setNbThreads(8);

    //================================
    grid.writeVTK(0);

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////
    ///////////////////              Solution            ////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    auto start1 = std::chrono::high_resolution_clock::now();
    for (int loadit = 0; loadit < steps; ++loadit)
    {
        
        auto start2 = std::chrono::high_resolution_clock::now();
        if(patch_BC_type==0) grid.Apply_BC_rigid(load_steps, load_step_patches);
        if(patch_BC_type==1) grid.Apply_BC_move_mirror(load_steps);
        
        int relaxIter = 1;  // for mirror two-way coupling
        if(patch_BC_type==1) relaxIter=3; 
        
        for (int r=0; r<relaxIter; ++r){
            if(patch_BC_type==1) grid.Apply_BC_mirror_patch();

            std::cout << std::endl;
            std::cout << "====================== Load Applied ============================= = " << loadit * load_steps << std::endl;

            auto end2 = std::chrono::high_resolution_clock::now();
            auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();

            double normnull;
            double normR;
            int error_counter = 1;
            bool isNotAccurate = true;
            double tol = 1e-10;

            auto start3 = std::chrono::high_resolution_clock::now();

            while (isNotAccurate)
            {
                std::vector<std::array<double, 2>> R(DOF / 2);
                std::vector<std::vector<std::array<double, 4>>> K_glob(DOF / 2);
                std::vector<std::vector<std::array<double, 4>>> K_glob_rearanged(DOF / 2);
                ///////// Computing R K x ////////////////////
                double V1_total = 0;
                double V2_total = 0;
                Compute_RK(grid.Nodes, grid.Compute_Nodes_idx, C1, C2, C3, horizon, V_Base, grid.max_Num_Neighbors, R, K_glob, nn, debug, deb_idx, V1_total, V2_total);
                rearrange_K_global(K_glob, K_glob_rearanged, grid.Compute_Nodes_idx, grid.Nodes);
                // std::cout << "V1_total = " << V1_total << std::endl;
                // std::cout << "V2_total = " << V2_total << std::endl;
                Solve_System(K_glob_rearanged, R, dx, Rtot, Ktot, solver);
                norm2_R = cal_norm2_R(Rtot, DOF);
                grid.update_compute_nodes(dx);
                if (error_counter == 1)
                    normnull = cal_norm2_R(Rtot, DOF);
                normR = cal_norm2_R(Rtot, DOF);
                std::cout << "Residual Norm @ Increment " << loadit + 1 << " @ Iteration " << error_counter << " with Residual :" << scientific << setprecision(2) << normR << ",  normalized :" << normR / normnull << std::endl;
                ++error_counter;

                if (error_counter > 1 && normR / normnull < tol)
                {
                    isNotAccurate = false;
                }
            }

            auto end3 = std::chrono::high_resolution_clock::now();
            auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();

            std::cout << "Iteration Runtime = " << duration3 << " milliseconds" << std::endl;

        } // end relax
        grid.writeVTK(loadit+1);
    }

    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
    std::cout <<  std::endl;
    std::cout << "Overall Iteration Runtime = " << duration1 << " milliseconds" << std::endl;
    std::cout <<  std::endl;

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////