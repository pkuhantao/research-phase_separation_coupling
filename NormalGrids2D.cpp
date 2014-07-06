//
//  NormalGrids2D.cpp
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#include "NormalGrids2D.h"
#include <cassert>
#include <cmath>

// constructor
NormalGrids2D::NormalGrids2D(int NTheta, int NPhi, double r) {
    // step
    dTheta = (thetaMax-thetaMin) / NTheta;
    dPhi = (phiMax-phiMin) / NPhi;
    // radius
    radius = r;
    // theta and phi value
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
    // calculate (x, y, z) in its own coordinates
    for (int i = 0; i < nrol; i++) {
        for (int j = 0; j < ncol; j++) {
            x0[i][j] = radius * sin(thetaVal[i]) * cos(phiVal[j]);
            y0[i][j] = radius * sin(thetaVal[i]) * sin(phiVal[j]);
            z0[i][j] = radius * cos(thetaVal[i]);
        }
    }
    // calculate (x, y, z) in its opposite coordinates
    for (int i = 0; i < nrol; i++) {
        for (int j = 0; j < ncol; j++) {
            x1[i][j] = -x0[i][j];
            y1[i][j] = z0[i][j];
            z1[i][j] = y0[i][j];
        }
    }
}

// theta
double NormalGrids2D::theta(int i) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    return thetaVal[i];
}

// phi
double NormalGrids2D::phi(int j) const {
    assert(j >= 0 && j < (int)phiVal.size());
    return phiVal[j];
}

// x
double NormalGrids2D::x(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x0[i][j];
}

// y
double NormalGrids2D::y(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y0[i][j];
}

// z
double NormalGrids2D::z(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z0[i][j];
}

// x'
double NormalGrids2D::xprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return x1[i][j];
}

// y'
double NormalGrids2D::yprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return y1[i][j];
}

// z'
double NormalGrids2D::zprime(int i, int j) const {
    assert(i >= 0 && i < (int)thetaVal.size());
    assert(j >= 0 && j < (int)phiVal.size());
    return z1[i][j];
}

// given theta angle, return the index for that angle(round downwards)
int NormalGrids2D::thetaToIndex(double angle) const {
    return (int)floor((angle-thetaMin)/dTheta + 0.000001);
}

// given phi angle, return the index for that angle(round downwards)
int NormalGrids2D::phiToIndex(double angle) const {
    return (int)floor((angle-phiMin)/dPhi + 0.000001);
}





