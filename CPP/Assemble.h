//
// Created by srini on 26/02/2025.
//

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

struct Point {
    vector<double> volumes;
    double Vol;
    double energy;
    vector<double> residual;
    vector<int> BCflg;
    vector<int> DOF;
    vector<int> neighbors;
    MatrixXd stiffness;
};

MatrixXd Assembly(const vector<Point>& PL, double Delta, int DOFs, double LF, const string& flag) {
    int NoPs = PL.size();
    int PD = PL[0].BCflg.size();
    MatrixXd out;

    if (flag == "point-volume") {
        out = MatrixXd::Zero(NoPs, 3);
        for (int p = 0; p < NoPs; p++) {
            out.row(p) = Map<const RowVectorXd>(PL[p].volumes.data(), 3);
        }
    } else if (flag == "volume") {
        double Vol = 0;
        for (int p = 0; p < NoPs; p++) {
            Vol += PL[p].Vol;
        }
        out = MatrixXd::Constant(1, 1, Vol);
    } else if (flag == "energy") {
        double Psi = 0;
        for (int p = 0; p < NoPs; p++) {
            double Psi_P = PL[p].energy;
            double Vol_P = PL[p].Vol;
            Psi += Vol_P * Psi_P;
        }
        out = MatrixXd::Constant(1, 1, Psi);
    } else if (flag == "residual") {
        VectorXd R = VectorXd::Zero(DOFs);
        for (int i = 0; i < NoPs; i++) {
            VectorXd R_P = Map<const VectorXd>(PL[i].residual.data(), PD);
            vector<int> BCflg = PL[i].BCflg;
            vector<int> DOF = PL[i].DOF;
            for (int ii = 0; ii < PD; ii++) {
                if (BCflg[ii] == 1) {
                    R(DOF[ii]) += R_P(ii);
                }
            }
        }
        out = R;
    } else if (flag == "stiffness") {
        vector<int> ii, jj;
        vector<double> vv;
        int sprC = 0;

        for (int p = 0; p < NoPs; p++) {
            MatrixXd K_P = PL[p].stiffness;
            vector<int> BCflg_p = PL[p].BCflg;
            vector<int> DOF_p = PL[p].DOF;

            for (int pp = 0; pp < PD; pp++) {
                if (BCflg_p[pp] == 1) {
                    vector<int> nbrL = PL[p].neighbors;
                    nbrL.push_back(p); // append the point to its neighbors!

                    for (size_t q = 0; q < nbrL.size(); q++) {
                        VectorXd K_PQtmp = K_P.col(q);
                        MatrixXd K_PQ;

                        if (PD == 2) {
                            K_PQ.resize(2, 2);
                            K_PQ << K_PQtmp(0), K_PQtmp(2),
                                     K_PQtmp(1), K_PQtmp(3);
                        } else if (PD == 3) {
                            K_PQ.resize(3, 3);
                            K_PQ << K_PQtmp(0), K_PQtmp(3), K_PQtmp(6),
                                     K_PQtmp(1), K_PQtmp(4), K_PQtmp(7),
                                     K_PQtmp(2), K_PQtmp(5), K_PQtmp(8);
                        }

                        vector<int> BCflg_q = PL[nbrL[q]].BCflg;
                        vector<int> DOF_q = PL[nbrL[q]].DOF;

                        for (int qq = 0; qq < PD; qq++) {
                            if (BCflg_q[qq] == 1) {
                                sprC++;
                                ii.push_back(DOF_p[pp]);
                                jj.push_back(DOF_q[qq]);
                                vv.push_back(K_PQ(pp, qq));
                            }
                        }
                    }
                }
            }
        }

        SparseMatrix<double> K(DOFs, DOFs);
        K.setFromTriplets(ii.begin(), jj.begin(), vv.begin());
        out = K;
    }

    return out;
}



#endif //ASSEMBLE_H
