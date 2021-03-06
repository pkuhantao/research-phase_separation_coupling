//
//  Membrane.cpp
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#include "Membrane.h"
#include "Properties.h"
#include "NormalGrids2D.h"
#include "MembranePatch.h"
#include "InnerSolvPatch.h"

#include <cassert>

using namespace std;

// constructor, given dt and membrane properties, create Yin/Yang patches and grids
Membrane::Membrane(double dt, MembProp membProp) {
    const double radius = membProp.radius;
    const int nTheta = membProp.nTheta;
    const int nPhi = membProp.nPhi;
    // create 2D Yin/Yang grids
    Yin2D = new NormalGrids2D(nTheta, nPhi, radius);
    Yang2D = new NormalGrids2D(nTheta, nPhi, radius);
    // create 2D Yin/Yang patches
    Yin = new MembranePatch(dt, membProp, Yin2D, Yang2D);
    Yang = new MembranePatch(dt, membProp, Yang2D, Yin2D);
}

Membrane::~Membrane() {
    delete Yin2D;
    delete Yang2D;
    delete Yin;
    delete Yang;
}

// initialize order parameter by Gaussian Distribution
void Membrane::initPsiGuass(const double ave, const double std) {
    Yin->initPsiGuass(ave, std);
    Yang->initPsiGuass(ave, std);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

// initialize order parameter by constant
void Membrane::initPsiConst(const double val) {
    Yin->initPsiConst(val);
    Yang->initPsiConst(val);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

 // initialize psi=inVal within the sphere with given radius and center (x0, y0, z0) in Yin grid, psi=otVal outside of sphere, effectively forms a disk on the membrane
void Membrane::initPsiDiskYinGd(double x0, double y0, double z0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    Yin->initPsiDiskSfGd(x0, y0, z0, radius, inVal, otVal);
    Yang->initPsiDiskOpGd(x0, y0, z0, radius, inVal, otVal);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

// calculate mu from double well potential
void Membrane::calcMu_dw() {
    Yin->calcMu_dw();
    Yang->calcMu_dw();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpMu();
}

// calculate mu from double well potential and the local coupling with the inner solvent
void Membrane::calcMu_dw_cp(const InnerSolvPatch* inSolv_Yin, const InnerSolvPatch* inSolv_Yang) {
    // calculate mu from double well potential and the local coupling
    Yin->calcMu_dw_cp(*(inSolv_Yin->psi_at_R_pt()));
    Yang->calcMu_dw_cp(*(inSolv_Yang->psi_at_R_pt()));
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpMu();
}



// calculate mu from simple diffusion model, where mu=psi
void Membrane::calcMu_sd() {
    Yin->calcMu_sd();
    Yang->calcMu_sd();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpMu();
}

// update psi from mu
void Membrane::updatePsi() {
    Yin->updatePsi();
    Yang->updatePsi();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpPsi();
}

// interpolate psi in ghost cells of Yin/Yang patches
void Membrane::interpPsi() {
    Yin->interpPsi(Yang->psi);
    Yang->interpPsi(Yin->psi);
}

// interpolate mu in ghost cells of Yin/Yang patches
void Membrane::interpMu() {
    Yin->interpMu(Yang->mu);
    Yang->interpMu(Yin->mu);
}



















