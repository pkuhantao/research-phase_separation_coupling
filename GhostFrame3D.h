//
//  GhostFrame3D.h
//  
//
//  Created by Tao Han on 4/30/14.
//
//

// 3D ghost frame surrounding Yin/Yang patch, which includes four internal ghost boundaries,
// one outer ghost boundary(contacting with membrane) and one inner ghost boundary(contacting with core)

#ifndef ____GhostFrame3D__
#define ____GhostFrame3D__


#include "GhostFrame2D.h"
#include "NormalGrids2D.h"
#include "NormalGrids3D.h"

#include <vector>

using namespace std;

class GhostFrame3D {
public:
    GhostFrame3D(const NormalGrids3D &currGrids3D, const NormalGrids2D &currGrids2D, const NormalGrids2D &oppoGrids2D); // given 2D and 3D grids, construct 3D ghost frame
    ~GhostFrame3D();
    // given the bulk value in both its own grids and the opposite grids and at the center, interpolate the value in all ghost cells
    void interpAllGhost(const vector<vector<vector<double> > > &f0, const vector<vector<vector<double> > > &f1, const double f_ct);
    
    vector<vector<double> > f0_ob, f0_ib; // value at outer and inner boundaries
    vector<vector<double> > f0_the_lo, f0_the_hi, f0_phi_lo, f0_phi_hi; // four internal ghost boundaries
    
private:
    // given the value in the bulk(its own grids), interpolate the value in the outer boundary by 'no flux' boundary condition
    void interpOuterBound(const vector<vector<vector<double> > > &f0);
    // given the value at the center, interpolate the value in the inner boundary
    void interpInnerBoundFromCt(const double f_ct);
    // given the value in the bulk(the opposite grids), interpolate the value in the internal ghosts
    // Note: here we interpolate the value from the in-plane cell, we don't consider the cells which are above or below the current plane
    void interpIntern(const vector<vector<vector<double> > > &f1);
    
    GhostFrame2D *ghost2D; // 2d internal ghost frame

};


#endif /* defined(____GhostFrame3D__) */
