//
//  NormalGrids3D.cpp
//  
//
//  Created by Tao Han on 4/30/14.
//
//

#include "NormalGrids3D.h"
#include <cassert>
#include <cmath>
#include <iostream>

// constructor
NormalGrids3D::NormalGrids3D(double ro, int NR, int NTheta, int NPhi) {
    ro_ = ro; // outer radius limit
    // step
    dR = ro_ / NR;
    dTheta = (thetaMax-thetaMin) / NTheta;
    dPhi = (phiMax-phiMin) / NPhi;
    // vector
    rVal.resize(NR, 0);
    thetaVal.resize(NTheta+1, 0);
    phiVal.resize(NPhi+1, 0);
    for (int i = 0; i < NR; i++) rVal[i] = dR + i*dR;
    for (int i = 0; i <= NTheta; i++) thetaVal[i] = thetaMin + i*dTheta;
    for (int i = 0; i <= NPhi; i++) phiVal[i] = phiMin + i*dPhi;
    
    const int n1 = rVal.size();
    const int n2 = thetaVal.size();
    const int n3 = phiVal.size();
    // (x, y, z) in its own coordinates
    x0.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    y0.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    z0.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    // (x, y, z) in its opposite coordinates
    x1.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    y1.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    z1.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    
    
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                x0[k][i][j] = rVal[k] * sin(thetaVal[i]) * cos(phiVal[j]);
                y0[k][i][j] = rVal[k] * sin(thetaVal[i]) * sin(phiVal[j]);
                z0[k][i][j] = rVal[k] * cos(thetaVal[i]);
            }
        }
    }
    
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                x1[k][i][j] = -x0[k][i][j];
                y1[k][i][j] = z0[k][i][j];
                z1[k][i][j] = y0[k][i][j];
            }
        }
    }
}


double NormalGrids3D::rad(int k) const {
    assert(k >= 0 && k < (int)rVal.size());
    return rVal[k];
}

double NormalGrids3D::theta(int i) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    return thetaVal[i];
}

double NormalGrids3D::phi(int i) const {
    assert(i >= 0 && i < (int)phiVal.size());
    return phiVal[i];
}


double NormalGrids3D::x(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x0[k][i][j];
}

double NormalGrids3D::y(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y0[k][i][j];
}

double NormalGrids3D::z(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z0[k][i][j];
}


double NormalGrids3D::xprime(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x1[k][i][j];
}

double NormalGrids3D::yprime(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y1[k][i][j];
}

double NormalGrids3D::zprime(int k, int i, int j) const {
    assert(k >= 0 && k < (int)rVal.size());
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z1[k][i][j];
}

int NormalGrids3D::radiusToIndex(double r) const {
    return (int)round((r-rVal[0])/dR);
}

int NormalGrids3D::thetaToIndex(double angle) const {
    return (int)floor((angle-thetaMin)/dTheta);
}

int NormalGrids3D::phiToIndex(double angle) const {
    return (int)floor((angle-phiMin)/dPhi);
}