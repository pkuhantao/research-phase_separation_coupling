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
    psi_Yin = &Yin.psi;
    psi_Yang = &Yang.psi;
    mu_Yin = &Yin.mu;
    mu_Yang = &Yang.mu;
    
    this->dt = dt;
    dr = Yin.grids3D->dRAD();
    // in the following, we find the positions of center's neighbors in Yin and Yang grids
    // there are three neighbors each in Yin and Yang grids, nbs_Yin[i][0..2] stands for index of (r, theta, phi)
    // in Yin grids for neighbor i
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
}

// calculate mu at the center from double well potential
void Center::calcMu_dw() {
    // calcualte laplacian of psi at first
    double laplace = 0;
    for (int i = 0; i < nbs_Yin.size(); i++) {
        int i1 = nbs_Yin[i][0];
        int i2 = nbs_Yin[i][1];
        int i3 = nbs_Yin[i][2];
        laplace += (*psi_Yin)[i1][i2][i3];
    }
    for (int i = 0; i < nbs_Yang.size(); i++) {
        int i1 = nbs_Yang[i][0];
        int i2 = nbs_Yang[i][1];
        int i3 = nbs_Yang[i][2];
        laplace += (*psi_Yang)[i1][i2][i3];
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
    for (int i = 0; i < nbs_Yin.size(); i++) {
        int i1 = nbs_Yin[i][0];
        int i2 = nbs_Yin[i][1];
        int i3 = nbs_Yin[i][2];
        laplace += (*mu_Yin)[i1][i2][i3];
    }
    for (int i = 0; i < nbs_Yang.size(); i++) {
        int i1 = nbs_Yang[i][0];
        int i2 = nbs_Yang[i][1];
        int i3 = nbs_Yang[i][2];
        laplace += (*mu_Yang)[i1][i2][i3];
    }
    laplace = (laplace-6.0*mu)/(dr*dr);
    
    psi += dt*props.M*laplace;
}






