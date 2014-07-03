//
//  GhostBC.cpp
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#include "GhostBC.h"
#include "NormalGrids.h"

#include <cmath>
#include <cassert>
#include <utility>
#include <vector>

GhostBC::GhostBC(const vector<double> &theta, const vector<double> &phi, const NormalGrids &oppoGrids) {
    const int N = theta.size();
    angles.resize(N, pair<double, double>());
    anglesPrim.resize(N, pair<double, double>());
    nbs.resize(N, vector<pair<int, int> >(numNbs, pair<int, int>())); // for each cell, 4 neighbors
    weights.resize(N, vector<double>(numNbs, 0));
    // make a copy of angles
    for (int i = 0; i < N; i++) {
        angles[i].first = theta[i];
        angles[i].second = phi[i];
    }
    // calculate corresponding angle in complemental coordinates
    coordiTrans();
    // find all neighbors in complemental coordinates
    findNbs(oppoGrids);
    // find all weights of each neighbor for interpolation
    calWts(oppoGrids);
}

// calcualte theta', phi' in complemental coordinate
void GhostBC::coordiTrans() {
    for (int i = 0; i < angles.size(); i++) {
        double theta = angles[i].first;
        double phi = angles[i].second;
        
        double thetaPrim = acos(sin(theta)*sin(phi));
        double phiPrim;
        if (fabs(sin(theta)*cos(phi)/sin(thetaPrim)) >= 1.0) phiPrim = 0;
        else if (cos(theta)/sin(thetaPrim) > 0) phiPrim = acos(-sin(theta)*cos(phi)/sin(thetaPrim));
        else phiPrim = -acos(-sin(theta)*cos(phi)/sin(thetaPrim));
        
        anglesPrim[i].first = thetaPrim;
        anglesPrim[i].second = phiPrim;
    }
}

// from given complemental normal grids, find neighbors for interpolation
void GhostBC::findNbs(const NormalGrids &grids) {
    for (int i = 0; i < anglesPrim.size(); i++) {
        double thetaPrim = anglesPrim[i].first;
        double phiPrim = anglesPrim[i].second;
        
        int rowID = grids.thetaToIndex(thetaPrim);
        int colID = grids.phiToIndex(phiPrim);
        nbs[i][0] = pair<int, int>(rowID, colID);
        nbs[i][1] = pair<int, int>(rowID, colID+1);
        nbs[i][2] = pair<int, int>(rowID+1, colID);
        nbs[i][3] = pair<int, int>(rowID+1, colID+1);
    }
}

// from the neighbors' position, calculate interpolation weights
void GhostBC::calWts(const NormalGrids &grids) {
    for (int i = 0; i < weights.size(); i++) {
        double thetaPrim = anglesPrim[i].first;
        double phiPrim = anglesPrim[i].second;
        double dThetaPrim = grids.dTHETA();
        double dPhiPrim = grids.dPHI();
        double thetaPrim0 = grids.theta(nbs[i][0].first);
        double phiPrim0 = grids.phi(nbs[i][0].second);
        double s0 = dThetaPrim * dPhiPrim;
        
        weights[i][0] = (thetaPrim0+dThetaPrim-thetaPrim) * (phiPrim0+dPhiPrim-phiPrim) / s0;
        weights[i][1] = (thetaPrim0+dThetaPrim-thetaPrim) * (phiPrim-phiPrim0) / s0;
        weights[i][2] = (thetaPrim-thetaPrim0) * (phiPrim0+dPhiPrim-phiPrim) / s0;
        weights[i][3] = (thetaPrim-thetaPrim0) * (phiPrim-phiPrim0) / s0;
    }
}

// given f, interpolate f0
void GhostBC::interpolation(vector<double> &f0, const vector<vector<double> > &f) {
    assert(f0.size() == weights.size());
    for (int i = 0; i < f0.size(); i++) {
        f0[i] = 0;
        for (int j = 0; j < numNbs; j++) {
            f0[i] += weights[i][j] * f[nbs[i][j].first][nbs[i][j].second];
        }
    }
}









