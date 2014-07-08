//
//  MembranePatch.h
//  
//
//  Created by Tao Han on 7/8/14.
//
//

// 2D Yin/Yang patch for membrane

#ifndef ____MembranePatch__
#define ____MembranePatch__

#include "Properties.h"
#include "NormalGrids2D.h"
#include "GhostFrame2D.h"

#include <vector>

using namespace std;

class MembranePatch {
public:
    MembranePatch(double dt, MembProp membProp, NormalGrids2D *cur2D, NormalGrids2D *opp2D);
    ~MembranePatch();
    
    void interpPsi(const vector<vector<double> > &oppPsi); // given psi in opposite grids, interpolate the ghost cells for psi in current grids
    void interpMu(const vector<vector<double> > &oppMu); // given mu in opposite grids, interpolate the ghost cells for mu in current grids
    
    void calcMu_dw(); // calculate mu from double well potential
    void calcMu_sd(); // calculate mu from simple diffusion model, where mu=psi
    void updatePsi(); // update psi from mu
    
    void initPsiGuass(const double ave, const double std); // initialize order parameter by Gaussian Distribution
    void initPsiConst(const double val); // initialize order parameter by constant
    
    vector<vector<double> > psi; // psi in the bulk
    vector<vector<double> > mu; // mu in the bulk

private:
    double dt; // time step
    MembProp props; // membrane properties
    
    const NormalGrids2D *grids2D; // 2D normal grids
    GhostFrame2D *ghost2D_psi, *ghost2D_mu; // 2D ghost frames for psi and mu
};



#endif /* defined(____MembranePatch__) */
