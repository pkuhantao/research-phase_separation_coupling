//
//  Membrane.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

// 2D membrane, which includes Yin and Yang patches

#ifndef ____Membrane__
#define ____Membrane__

#include "Properties.h"
#include "NormalGrids2D.h"
#include "MembranePatch.h"
#include "InnerSolv.h"

using namespace std;

class Membrane {
public:
    Membrane(double dt, MembProp membProp); // constructor, given dt and membrane properties, create Yin/Yang patches and grids
    ~Membrane();
    
    void initPsiGuass(const double ave, const double std); // initialize order parameter by Gaussian Distribution
    void initPsiConst(const double val); // initialize order parameter by constant
    void initPsiDiskYinGd(double x0, double y0, double z0, double radius, double inVal, double otVal); // initialize psi=inVal within the sphere with given radius and center (x0, y0, z0) in Yin grid, otVal outside of sphere, effectively forms a disk on the membrane
    
    void calcMu_dw(); // calculate mu from double well potential
    void calcMu_dw_cp(const InnerSolv &inSolv); // calculate mu from double well potential and the local coupling with the inner solvent
    void calcMu_sd(); // calculate mu from simple diffusion model, where mu=psi
    void updatePsi(); // update psi from mu
    
    const MembranePatch* YinPart() const {return Yin;}; // Yin part
    const MembranePatch* YangPart() const {return Yang;}; // Yang part
    const NormalGrids2D* YinGrids() const {return Yin2D;}; // Yin grids
    const NormalGrids2D* YangGrids() const {return Yang2D;}; // Yang grids
    
private:
    void interpPsi(); // interpolate psi in ghost cells of Yin/Yang patches
    void interpMu();  // interpolate mu in ghost cells of Yin/Yang patches
    
    MembranePatch *Yin, *Yang; // Yin, Yang patches
    NormalGrids2D *Yin2D, *Yang2D; // 2D grids for Yin and Yang
};




#endif /* defined(____Membrane__) */
