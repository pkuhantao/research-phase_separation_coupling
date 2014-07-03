//
//  InnerSolv.cpp
//  
//
//  Created by Tao Han on 5/14/14.
//
//

#include "InnerSolv.h"
#include "NormalGrids.h"
#include "NormalGrids3D.h"
#include "Properties.h"
#include "Center.h"
#include "InnerSolvPatch.h"

#include <cassert>
#include <cmath>

InnerSolv::InnerSolv(double dt) {
    double M = 1;
    double a = 1;
    double b = 1;
    double w = 1;
    SolvProp props(M, w, a, b);
    
    double ro = 10;
    const int n1 = 5;
    const int n2 = 20;
    const int n3 = 24;
    Yin2D = new NormalGrids(n2, n3, ro);
    Yang2D = new NormalGrids(n2, n3, ro);
    Yin3D = new NormalGrids3D(ro, n1, n2, n3);
    Yang3D = new NormalGrids3D(ro, n1, n2, n3);
    Yin = new InnerSolvPatch(dt, props, Yin3D, Yin2D, Yang2D);
    Yang = new InnerSolvPatch(dt, props, Yang3D, Yang2D, Yin2D);
    ct = new Center(dt, props, *Yin, *Yang);
}

InnerSolv::~InnerSolv() {
    delete Yin2D;
    delete Yang2D;
    delete Yin3D;
    delete Yang3D;
    delete Yin;
    delete Yang;
    delete ct;
}

// initialize psi by random variable from Gauss distribution
void InnerSolv::initPsiGauss(const double ave, const double std) {
    Yin->initPsiGuass(ave, std);
    Yang->initPsiGuass(ave, std);
    ct->initPsi(ave);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

// initialize psi by inVal within the sphere with given radius, otVal outside of sphere
void InnerSolv::initPsiSphere(double radius, double inVal, double otVal) {
    Yin->initPsiSphere(radius, inVal, otVal);
    Yang->initPsiSphere(radius, inVal, otVal);
    if (radius > 0) ct->initPsi(inVal);
    else ct->initPsi(otVal);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

// initialize psi=inVal within the sphere with given radius and center (x0, y0, z0) in Yin grid, otVal outside of sphere
void InnerSolv::initPsiSprYinGd(double x0, double y0, double z0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    Yin->initPsiSprSfGd(x0, y0, z0, radius, inVal, otVal);
    Yang->initPsiSprOpGd(x0, y0, z0, radius, inVal, otVal);
    double dist = sqrt(x0*x0 + y0*y0 + z0*z0); // distance to the origin
    if (dist < radius) ct->initPsi(inVal);
    else ct->initPsi(otVal);
    // after initialization, we interpolate the value in ghost cells automatically
    interpPsi();
}

// calculate chemical potential by double-well potential
void InnerSolv::calcMu_dw() {
    Yin->calcMu_dw(ct->psiAtCenter());
    Yang->calcMu_dw(ct->psiAtCenter());
    ct->calcMu_dw();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpMu();
}

// calculate chemical potential by simple diffusion model, where mu=psi
void InnerSolv::calcMu_sd() {
    Yin->calcMu_sd();
    Yang->calcMu_sd();
    ct->calcMu_sd();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpMu();
}

void InnerSolv::updatePsi() {
    Yin->updatePsi(ct->muAtCenter());
    Yang->updatePsi(ct->muAtCenter());
    ct->updatePsi();
    // after calculate the value in the bulk, we interpolate the value in ghost cells automatically
    interpPsi();
}

void InnerSolv::interpPsi() {
    Yin->interpPsi(Yang->psi);
    Yang->interpPsi(Yin->psi);
}

void InnerSolv::interpMu() {
    Yin->interpMu(Yang->mu);
    Yang->interpMu(Yin->mu);
}




