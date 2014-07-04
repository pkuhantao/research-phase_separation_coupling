//
//  GhostFrame2D.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#ifndef ____GhostFrame2D__
#define ____GhostFrame2D__

#include "GhostBC2D.h"
#include "NormalGrids2D.h"

#include <iostream>

using namespace std;

class GhostFrame2D {
public:
    GhostFrame2D(const NormalGrids2D &currGrids, const NormalGrids2D &oppoGrids);
    ~GhostFrame2D();
    void interpolate(vector<double> &f0_the_lo, vector<double> &f0_the_hi, vector<double> &f0_phi_lo, vector<double> &f0_phi_hi, const vector<vector<double> > &f); // from f in complemental grid, interpolate values in four boundaries of current grid
    
private:
    GhostBC2D *theta_lo, *theta_hi, *phi_lo, *phi_hi;
};


#endif /* defined(____GhostFrame2D__) */
