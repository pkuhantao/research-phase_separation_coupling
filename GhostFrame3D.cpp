//
//  GhostFrame3D.cpp
//  
//
//  Created by Tao Han on 5/1/14.
//
//

#include "GhostFrame3D.h"
#include "GhostFrame2D.h"
#include "NormalGrids2D.h"
#include "NormalGrids3D.h"

#include <vector>

using namespace std;

// given 2D and 3D grids, construct 3D ghost frame
GhostFrame3D::GhostFrame3D(const NormalGrids3D &currGrids3D, const NormalGrids2D &currGrids2D, const NormalGrids2D &oppoGrids2D) {
    // construct 2D ghost frame, used for internal ghost boundary interpolation
    ghost2D = new GhostFrame2D(currGrids2D, oppoGrids2D);
    // create ghost boundaries
    const int Nr = currGrids3D.sizeInR();
    const int Ntheta = currGrids3D.sizeInTheta();
    const int Nphi = currGrids3D.sizeInPhi();
    
    f0_ob.resize(Ntheta, vector<double>(Nphi, 0)); // outer boundary
    f0_the_lo.resize(Nr, vector<double>(Nphi, 0)); // four internal ghost boundaries
    f0_the_hi.resize(Nr, vector<double>(Nphi, 0));
    f0_phi_lo.resize(Nr, vector<double>(Ntheta, 0));
    f0_phi_hi.resize(Nr, vector<double>(Ntheta, 0));
}

// destructor
GhostFrame3D::~GhostFrame3D() {
    delete ghost2D;
}

// given the bulk value in both its own grids and the opposite grids, interpolate the value in all ghost cells
void GhostFrame3D::interpAllGhost(const vector<vector<vector<double> > > &f0, const vector<vector<vector<double> > > &f1) {
    interpOuterBound(f0); // given the bulk value in its own grids, interpolate outer boundary
    interpIntern(f1); // given the bulk value in its opposite grids, interpolate four internal boudaries
}


// given the value in the bulk(its own grids), interpolate the value in the outer boundary
void GhostFrame3D::interpOuterBound(const vector<vector<vector<double> > > &f0) {
    const int n1 = f0.size();
    const int n2 = f0[0].size();
    const int n3 = f0[0][0].size();
    
    // here we just use the no flow boundary condition
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n3; j++) {
            f0_ob[i][j] = f0[n1-2][i][j]; // n1-2 is the second largest radius
        }
    }
}

// given the value in the bulk(the opposite grids), interpolate the value in the internal ghosts
// Note: here we interpolate the value from the in-plane cell, we don't consider the cells which are above or below the current plane
void GhostFrame3D::interpIntern(const vector<vector<vector<double> > > &f1) {
    const int n1 = f1.size(); // radius dimension
    // interpolate the internal boundaries shell by shell
    for (int k = 0; k < n1; k++) {
        ghost2D->interpolate(f0_the_lo[k], f0_the_hi[k], f0_phi_lo[k], f0_phi_hi[k], f1[k]);
    }
}














