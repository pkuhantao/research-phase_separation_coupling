//
//  GhostFrame2D.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

// 2D ghost frame surrounding Yin/Yang patch which includes four ghost boundaries
// enabling the interpolation at the Yin/Yang patch boundaries

#ifndef ____GhostFrame2D__
#define ____GhostFrame2D__

#include "GhostBC2D.h"
#include "NormalGrids2D.h"

using namespace std;

class GhostFrame2D {
public:
    // given current and opposite Yin/Yang grids, build up interpolation scheme
    GhostFrame2D(const NormalGrids2D &currGrids, const NormalGrids2D &oppoGrids);
    ~GhostFrame2D();
    void interpolate(vector<double> &f0_the_lo, vector<double> &f0_the_hi, vector<double> &f0_phi_lo, vector<double> &f0_phi_hi, const vector<vector<double> > &f) const; // given value f in complemental Yin/Yang grids, interpolate values in four boundaries of current Yin/Yang grids, store into pass-in arguments
    void interpolate(const vector<vector<double> > &f); // given value f in complemental Yin/Yang grids, interpolate values in four boundaries of current Yin/Yang grids, store into its own member variables
    
    vector<double> f0_the_lo_, f0_the_hi_, f0_phi_lo_, f0_phi_hi_; // value in four ghost boundaries
    
private:
    GhostBC2D *theta_lo, *theta_hi, *phi_lo, *phi_hi; // four ghost boundaries
};


#endif /* defined(____GhostFrame2D__) */
