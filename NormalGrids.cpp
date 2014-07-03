//
//  NormalGrids.cpp
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#include "NormalGrids.h"
#include <cassert>
#include <cmath>
#include <iostream>

// constructor
NormalGrids::NormalGrids(int NTheta, int NPhi, double r) {
    // step
    dTheta = (thetaMax-thetaMin) / NTheta;
    dPhi = (phiMax-phiMin) / NPhi;
    radius = r;
    // vector
    thetaVal.resize(NTheta+1, 0);
    phiVal.resize(NPhi+1, 0);
    for (int i = 0; i <= NTheta; i++) thetaVal[i] = thetaMin + i*dTheta;
    for (int i = 0; i <= NPhi; i++) phiVal[i] = phiMin + i*dPhi;
    
    const int nrol = thetaVal.size();
    const int ncol = phiVal.size();
    // (x, y, z) in its own coordinates
    x0.resize(nrol, vector<double>(ncol, 0));
    y0.resize(nrol, vector<double>(ncol, 0));
    z0.resize(nrol, vector<double>(ncol, 0));
    // (x, y, z) in its opposite coordinates
    x1.resize(nrol, vector<double>(ncol, 0));
    y1.resize(nrol, vector<double>(ncol, 0));
    z1.resize(nrol, vector<double>(ncol, 0));
    for (int i = 0; i < nrol; i++) {
        for (int j = 0; j < ncol; j++) {
            x0[i][j] = radius * sin(thetaVal[i]) * cos(phiVal[j]);
            y0[i][j] = radius * sin(thetaVal[i]) * sin(phiVal[j]);
            z0[i][j] = radius * cos(thetaVal[i]);
        }
    }
    for (int i = 0; i < nrol; i++) {
        for (int j = 0; j < ncol; j++) {
            x1[i][j] = -x0[i][j];
            y1[i][j] = z0[i][j];
            z1[i][j] = y0[i][j];
        }
    }
}

double NormalGrids::theta(int i) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    return thetaVal[i];
}


double NormalGrids::phi(int i) const {
    assert(i >= 0 && i < (int)phiVal.size());
    return phiVal[i];
}

double NormalGrids::x(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x0[i][j];
}

double NormalGrids::y(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y0[i][j];
}

double NormalGrids::z(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z0[i][j];
}

double NormalGrids::xprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x1[i][j];
}

double NormalGrids::yprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y1[i][j];
}

double NormalGrids::zprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z1[i][j];
}


int NormalGrids::thetaToIndex(double angle) const {
    return (int)floor((angle-thetaMin)/dTheta);
}

int NormalGrids::phiToIndex(double angle) const {
    return (int)floor((angle-phiMin)/dPhi);
}