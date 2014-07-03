//
//  GhostFrame3D.h
//  
//
//  Created by Tao Han on 4/30/14.
//
//

#ifndef ____GhostFrame3D__
#define ____GhostFrame3D__


#include "GhostFrame.h"
#include "NormalGrids.h"
#include "NormalGrids3D.h"

#include <iostream>
#include <vector>

using namespace std;

class GhostFrame3D {
public:
    GhostFrame3D(const NormalGrids3D &currGrids3D, const NormalGrids &currGrids2D, const NormalGrids &oppoGrids2D); // given 2D and 3D grids, construct 3D ghost frame
    ~GhostFrame3D();
    
    // given the bulk value in both its own grids and the opposite grids, interpolate the value in all ghost cells
    void interpAllGhost(const vector<vector<vector<double> > > &f0, const vector<vector<vector<double> > > &f1);
    
    vector<vector<double> > f0_ob; // outer boundary
    vector<vector<double> > f0_the_lo, f0_the_hi, f0_phi_lo, f0_phi_hi; // four internal ghost boundaries
    
private:
    // given the value in the bulk(its own grids), interpolate the value in the outer boundary
    void interpOuterBound(const vector<vector<vector<double> > > &f0);
    // given the value in the bulk(the opposite grids), interpolate the value in the internal ghosts
    // Note: here we interpolate the value from the in-plane cell, we don't consider the cells which are above or below the current plane
    void interpIntern(const vector<vector<vector<double> > > &f1);
    
    GhostFrame *ghost2D; // 2d internal ghost frame

};


#endif /* defined(____GhostFrame3D__) */
