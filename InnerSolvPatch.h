//
//  InnerSolvPatch.h
//  
//
//  Created by Tao Han on 5/2/14.
//
//

// 3D complex inner solvent patch(Yin/Yang)

#ifndef ____InnerSolvPatch__
#define ____InnerSolvPatch__


#include "NormalGrids2D.h"
#include "NormalGrids3D.h"
#include "GhostFrame3D.h"
#include "Properties.h"

#include <vector>

using namespace std;

class InnerSolvPatch {
public:
    InnerSolvPatch(double dt, SolvProp solvProp, NormalGrids3D *cur3D, NormalGrids2D *cur2D, NormalGrids2D *opp2D);
    ~InnerSolvPatch();
    
    void interpPsi(const vector<vector<vector<double> > > &oppPsi, const double psiCt); // given psi in opposite grids and at the center, interpolate the ghost cells for psi in current grids
    void interpMu(const vector<vector<vector<double> > > &oppMu, const double muCt); // given mu in opposite grids and at the center, interpolate the ghost cells for mu in current grids
    
    void calcMu_dw(); // calculate mu from double well potential
    void calcMu_sd(); // calculate mu from simple diffusion model, where mu=psi
    void updatePsi(); // update psi from mu
    
    void initPsiGuass(const double ave, const double std); // given average and standard deviation, initialize psi by Guassian Random noise
    void initPsiSphere(double radius, double inVal, double otVal); // initialize psi=inVal in sphere with given radius and center(0, 0, 0), psi=otVal outside
    void initPsiSprSfGd(double x0, double y0, double z0, double radius, double inVal, double otVal); // initialize psi=inVal in sphere with given radius and center(x, y, z) in self grid, psi=otVal outside
    void initPsiSprOpGd(double xp0, double yp0, double zp0, double radius, double inVal, double otVal); // initialize psi=inVal in sphere with given radius and center(xp, yp, zp) in opposite grid, psi=otVal outside
    
    
    vector<vector<vector<double> > > psi; // psi in the bulk
    vector<vector<vector<double> > > mu; // mu in the bulk
    const NormalGrids3D *grids3D; // current 3D normal grids
    
private:
    SolvProp props; // inner solvent properties
    double dt; // time step
    
    GhostFrame3D *ghost3D_psi; // 3D ghost frame for psi
    GhostFrame3D *ghost3D_mu; // 3D ghost frame for mu
    
};










#endif /* defined(____InnerSolvPatch__) */
