//
//  GhostFrame2D.cpp
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#include "GhostFrame2D.h"
#include "GhostBC2D.h"
#include "NormalGrids2D.h"

#include <vector>


GhostFrame2D::GhostFrame2D(const NormalGrids2D &currGrids, const NormalGrids2D &oppoGrids) {
    const int sizeTheta = currGrids.sizeInTheta();
    const int sizePhi = currGrids.sizeInPhi();
    const double dtheta = currGrids.dTHETA();
    const double dphi = currGrids.dPHI();
    
    vector<double> thetaVal(sizeTheta, 0);
    for (int i = 0; i < sizeTheta; i++) thetaVal[i] = currGrids.theta(i);
    vector<double> phiVal(sizePhi, 0);
    for (int i = 0; i < sizePhi; i++) phiVal[i] = currGrids.phi(i);
    // create four ghost boundaries
    theta_lo = new GhostBC2D(vector<double>(sizePhi, thetaVal[0]-dtheta), phiVal, oppoGrids);
    theta_hi = new GhostBC2D(vector<double>(sizePhi, thetaVal[sizeTheta-1]+dtheta), phiVal, oppoGrids);
    phi_lo = new GhostBC2D(thetaVal, vector<double>(sizeTheta, phiVal[0]-dphi), oppoGrids);
    phi_hi = new GhostBC2D(thetaVal, vector<double>(sizeTheta, phiVal[sizePhi-1]+dphi), oppoGrids);
}

GhostFrame2D::~GhostFrame2D() {
    delete theta_lo;
    delete theta_hi;
    delete phi_lo;
    delete phi_hi;
}

// from f in complemental grid, interpolate values in four boundaries of current grid
void GhostFrame2D::interpolate(vector<double> &f0_the_lo, vector<double> &f0_the_hi, vector<double> &f0_phi_lo, vector<double> &f0_phi_hi, const vector<vector<double> > &f) {
    theta_lo->interpolation(f0_the_lo, f);
    theta_hi->interpolation(f0_the_hi, f);
    phi_lo->interpolation(f0_phi_lo, f);
    phi_hi->interpolation(f0_phi_hi, f);
}
