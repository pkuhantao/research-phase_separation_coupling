//
//  InnerSolv.h
//  
//
//  Created by Tao Han on 5/14/14.
//
//

// Inner solvent

#ifndef ____InnerSolv__
#define ____InnerSolv__


#include "NormalGrids2D.h"
#include "NormalGrids3D.h"
#include "Properties.h"
#include "Center.h"
#include "InnerSolvPatch.h"
#include "Membrane.h"


class InnerSolv {
public:
    InnerSolv(double dt, SolvProp props); // constructor, given dt and inner solvent properties, create grids, Yin/Yang patches and center
    ~InnerSolv();
    
    void initPsiGauss(const double ave, const double std); // initialize psi by random variable from Gauss distribution
    void initPsiSphere(double radius, double inVal, double otVal); // initialize psi by inVal within the sphere with given radius, otVal outside of sphere
    void initPsiSprYinGd(double x0, double y0, double z0, double radius, double inVal, double otVal); // initialize psi=inVal within the sphere with given radius and center (x0, y0, z0) in Yin grid, otVal outside of sphere
    
    void calcMu_dw(); // calculate chemical potential by double-well potential
    void calcMu_dw_cp(const Membrane &memb); // calculate mu from double well potential and the local coupling with the membrane
    void calcMu_sd(); // calculate chemical potential by simple diffusion model, where mu=psi
    void updatePsi();
    
    const InnerSolvPatch* YinPart() const {return Yin;};
    const InnerSolvPatch* YangPart() const {return Yang;};
    const Center* CenterPart() const {return ct;};
    
private:
    void interpPsi(); // interpolate psi in ghost cells of Yin/Yang patches and the center
    void interpMu();  // interpolate mu in ghost cells of Yin/Yang patches and the center
    
    InnerSolvPatch *Yin, *Yang; // Yin, Yang patches
    Center *ct; // center
    NormalGrids2D *Yin2D, *Yang2D; // 2D grids for Yin and Yang
    NormalGrids3D *Yin3D, *Yang3D; // 3D grids for Yin and Yang
};




#endif /* defined(____InnerSolv__) */
