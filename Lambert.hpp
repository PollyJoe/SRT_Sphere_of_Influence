//
//  Lambert.hpp
//  Lambert's Problem (Self-Version)
//
//  Created by 周懿 on 2021/11/22.
//

#ifndef _LAMBERT_H_
#define _LAMBERT_H_

//Subfunction that evaluates the time of flight as a function of x
double x2tof(double x, double s, double c, int lw, int N);

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
    double tf, const double* unith, int& flag, double mu, int way = 0, int N = 0, int branch = 0,
    int Maxiter = 200, double tol = 1.0e-11);

#endif
