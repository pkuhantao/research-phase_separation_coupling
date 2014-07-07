//
//  Center.cpp
//  
//
//  Created by Tao Han on 5/2/14.
//
//

#include "Center.h"
#include "Properties.h"
#include "InnerSolvPatch.h"

#include <vector>


#define PI 3.14159265


using namespace std;

Center::Center(double dt, SolvProp solvProp, const InnerSolvPatch &Yin, const InnerSolvPatch &Yang) {
    props = solvProp;
    this->dt = dt;
    dr = Yin.grids3D->dRAD();
    // In the following, we find the positions of center's neighbors in Yin and Yang grids
    // There are three neighbors each in Yin and Yang grids,
    // nbs_Yin[i][0..2] stands for index of (r, theta, phi) in Yin grids for neighbor i
    nbs_Yin.resize(3, vector<int>(3, 0));
    nbs_Yang.resize(3, vector<int>(3, 0));
    // neighbors in Yin grids
    int rID = Yin.grids3D->radiusToIndex(dr);
    int thetaID = Yin.grids3D->thetaToIndex(PI/2.0);
    for (int index = 0; index < 3; index++) {
        double angle = (index-1)*PI/2.0;
        int phiID = Yin.grids3D->phiToIndex(angle);
        nbs_Yin[index][0] = rID;
        nbs_Yin[index][1] = thetaID;
        nbs_Yin[index][2] = phiID;
    }
    // neighbors in Yang grids
    rID = Yang.grids3D->radiusToIndex(dr);
    thetaID = Yang.grids3D->thetaToIndex(PI/2.0);
    for (int index = 0; index < 3; index++) {
        double angle = (index-1)*PI/2.0;
        int phiID = Yang.grids3D->phiToIndex(angle);
        nbs_Yang[index][0] = rID;
        nbs_Yang[index][1] = thetaID;
        nbs_Yang[index][2] = phiID;
    }
    // create vectors to store psi and mu in neighbor cells
    psi_Yin_nbs.resize(3, 0);
    psi_Yang_nbs.resize(3, 0);
    mu_Yin_nbs.resize(3, 0);
    mu_Yang_nbs.resize(3, 0);
}

// interpolate psi in neighboring cells from Yin and Yang grids
void Center::interpPsi(const vector<vector<vector<double> > > &psi_Yin, const vector<vector<vector<double> > > &psi_Yang) {
    // psi from Yin grids
    for (int i = 0; i < nbs_Yin.size(); i++) {
        int i1 = nbs_Yin[i][0];
        int i2 = nbs_Yin[i][1];
        int i3 = nbs_Yin[i][2];
        psi_Yin_nbs[i] = psi_Yin[i1][i2][i3];
    }
    // psi from Yang grids
    for (int i = 0; i < nbs_Yang.size(); i++) {
        int i1 = nbs_Yang[i][0];
        int i2 = nbs_Yang[i][1];
        int i3 = nbs_Yang[i][2];
        psi_Yang_nbs[i] = psi_Yang[i1][i2][i3];
    }
}

// interpolate mu in neighboring cells from Yin and Yang grids
void Center::interpMu(const vector<vector<vector<double> > > &mu_Yin, const vector<vector<vector<double> > > &mu_Yang) {
    // mu from Yin grids
    for (int i = 0; i < nbs_Yin.size(); i++) {
        int i1 = nbs_Yin[i][0];
        int i2 = nbs_Yin[i][1];
        int i3 = nbs_Yin[i][2];
        mu_Yin_nbs[i] = mu_Yin[i1][i2][i3];
    }
    // mu from Yang grids
    for (int i = 0; i < nbs_Yang.size(); i++) {
        int i1 = nbs_Yang[i][0];
        int i2 = nbs_Yang[i][1];
        int i3 = nbs_Yang[i][2];
        mu_Yang_nbs[i] = mu_Yang[i1][i2][i3];
    }
}

// calculate mu at the center from double well potential
void Center::calcMu_dw() {
    // calcualte laplacian of psi at first
    double laplace = 0;
    for (int i = 0; i < psi_Yin_nbs.size(); i++) {
        laplace += psi_Yin_nbs[i];
    }
    for (int i = 0; i < psi_Yang_nbs.size(); i++) {
        laplace += psi_Yang_nbs[i];
    }
    laplace = (laplace-6.0*psi)/(dr*dr);
    
    mu = -(props.w*props.w/2.0)*laplace - props.a*psi + props.b*psi*psi*psi;
}

// calculate mu at the center from simple diffusion model, where mu = psi
void Center::calcMu_sd() {
    mu = psi;
}

// update psi at the center from mu
void Center::updatePsi() {
    // calcualte laplacian of mu at first
    double laplace = 0;
    for (int i = 0; i < mu_Yin_nbs.size(); i++) {
        laplace += mu_Yin_nbs[i];
    }
    for (int i = 0; i < mu_Yang_nbs.size(); i++) {
        laplace += mu_Yang_nbs[i];
    }
    laplace = (laplace-6.0*mu)/(dr*dr);
    
    psi += dt*props.M*laplace;
}






