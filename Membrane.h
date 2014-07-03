//
//  Membrane.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#ifndef ____Membrane__
#define ____Membrane__

#include "NormalGrids.h"
#include "GhostFrame.h"

#include <vector>
#include <iostream>

using namespace std;

class Membrane {
public:
    Membrane(double M, double dt, NormalGrids *grids, GhostFrame *frame);
    ~Membrane(){};
    void initGauRand(double ave, double std); // initialize order parameter by Gaussian Distribution
    void initDisk(int ic, int jc, int rad); // initialize order parameter by disk(+1) surrounded by -1
    void interpPsi(const vector<vector<double> > &oppoPsi); // interpolate psi in ghost cells by given psi in complemental grids
    void interpMu(const vector<vector<double> > &oppoMu); // interpolate mu in ghost cells by given mu in complemental grids
    void calcMu_simdiff(); // calculate mu from psi, simple diffusion model (mu=psi)
    void calcMu_doublewell(); // calculate mu from double well potential
    void updatePsi(); // update psi from mu
    
    vector<vector<double> > psi;
    vector<vector<double> > mu;
    
private:
    static const double w = 1.0;
    static const double a = 1.0;
    static const double b = 1.0;
    double M; // diffusivity
    double dt; // time step
    NormalGrids *normGrids; // normal grids
    GhostFrame *ghostFm; // ghost frame
    vector<double> psi_the_lo, psi_the_hi, psi_phi_lo, psi_phi_hi; // ghost cells for order parameter
    vector<double> mu_the_lo, mu_the_hi, mu_phi_lo, mu_phi_hi; // ghost cells for chemical potential
};

#endif /* defined(____Membrane__) */
