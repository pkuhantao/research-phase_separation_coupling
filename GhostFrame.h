//
//  GhostFrame.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#ifndef ____GhostFrame__
#define ____GhostFrame__

#include "GhostBC.h"
#include "NormalGrids.h"

#include <iostream>

using namespace std;

class GhostFrame {
public:
    GhostFrame(const NormalGrids &currGrids, const NormalGrids &oppoGrids);
    ~GhostFrame();
    void interpolate(vector<double> &f0_the_lo, vector<double> &f0_the_hi, vector<double> &f0_phi_lo, vector<double> &f0_phi_hi, const vector<vector<double> > &f); // from f in complemental grid, interpolate values in four boundaries of current grid
    
private:
    GhostBC *theta_lo, *theta_hi, *phi_lo, *phi_hi;
};


#endif /* defined(____GhostFrame__) */
