//
//  InnerSolvPatch.cpp
//  
//
//  Created by Tao Han on 5/2/14.
//
//

#include "InnerSolvPatch.h"
#include "NormalGrids2D.h"
#include "NormalGrids3D.h"
#include "GhostFrame3D.h"
#include "Properties.h"

#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;

InnerSolvPatch::InnerSolvPatch(double dt, SolvProp solvProp, NormalGrids3D *cur3D, NormalGrids2D *cur2D, NormalGrids2D *opp2D): grids3D(cur3D) {
    props = solvProp;
    this->dt = dt;
    // create matrix for psi and mu
    const int n1 = cur3D->sizeInR();
    const int n2 = cur3D->sizeInTheta();
    const int n3 = cur3D->sizeInPhi();
    psi.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    mu.resize(n1, vector<vector<double> >(n2, vector<double>(n3, 0)));
    // create ghost cells
    ghost3D_psi = new GhostFrame3D(*cur3D, *cur2D, *opp2D);
    ghost3D_mu = new GhostFrame3D(*cur3D, *cur2D, *opp2D);
}

InnerSolvPatch::~InnerSolvPatch() {
    delete ghost3D_psi;
    delete ghost3D_mu;
}

// given psi in opposite grids and at the center, interpolate the ghost cells for psi in current grids
void InnerSolvPatch::interpPsi(const vector<vector<vector<double> > > &oppPsi, const double psiCt) {
    ghost3D_psi->interpAllGhost(psi, oppPsi, psiCt);
}

// given mu in opposite grids and at the center, interpolate the ghost cells for mu in current grids
void InnerSolvPatch::interpMu(const vector<vector<vector<double> > > &oppMu, const double muCt) {
    ghost3D_mu->interpAllGhost(mu, oppMu, muCt);
}

// calculate mu from double well potential
void InnerSolvPatch::calcMu_dw() {
    // dimensions
    const int n1 = mu.size();
    const int n2 = mu[0].size();
    const int n3 = mu[0][0].size();
    // steps
    const double dR = grids3D->dRAD();
    const double dTheta = grids3D->dTHETA();
    const double dPhi = grids3D->dPHI();
    
    for (int k = 0; k < n1; k++) {
        double r = grids3D->rad(k); // r_k
        
        for (int i = 0; i < n2; i++) {
            double theta = grids3D->theta(i); // theta_i
            double sin_theta = sin(theta); // sin(theta_i)
            double tan_theta = tan(theta); // tan(theta_i)
            
            for (int j = 0; j < n3; j++) {
                // calculate laplace
                double r_lo = (k == 0) ? ghost3D_psi->f0_ib[i][j] : psi[k-1][i][j];
                double r_hi = (k == (n1-1)) ? ghost3D_psi->f0_ob[i][j] : psi[k+1][i][j];
                double the_lo = (i == 0) ? ghost3D_psi->f0_the_lo[k][j] : psi[k][i-1][j];
                double the_hi = (i == (n2-1)) ? ghost3D_psi->f0_the_hi[k][j] : psi[k][i+1][j];
                double phi_lo = (j == 0) ? ghost3D_psi->f0_phi_lo[k][i] : psi[k][i][j-1];
                double phi_hi = (j == (n3-1)) ? ghost3D_psi->f0_phi_hi[k][i] : psi[k][i][j+1];
                // in-plane derivative
                double laplace = (the_hi+the_lo-2.0*psi[k][i][j])/(dTheta*dTheta) + (phi_hi+phi_lo-2.0*psi[k][i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (the_hi-the_lo)/(2.0*tan_theta*dTheta);
                laplace /= (r*r);
                // out-of-plane derivative
                laplace += (r_hi+r_lo-2.0*psi[k][i][j])/(dR*dR) + (r_hi-r_lo)/(r*dR);
                
                mu[k][i][j] = -(props.w*props.w/2.0)*laplace - props.a*psi[k][i][j] + props.b*psi[k][i][j]*psi[k][i][j]*psi[k][i][j];
            }
        }
    }
}

// calculate mu in the bulk from simple diffusion model
void InnerSolvPatch::calcMu_sd() {
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                mu[k][i][j] = psi[k][i][j];
            }
        }
    }
}

// update psi in the bulk
void InnerSolvPatch::updatePsi() {
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    // steps
    const double dR = grids3D->dRAD();
    const double dTheta = grids3D->dTHETA();
    const double dPhi = grids3D->dPHI();
    
    for (int k = 0; k < n1; k++) {
        double r = grids3D->rad(k); // r_k
        
        for (int i = 0; i < n2; i++) {
            double theta = grids3D->theta(i); // theta_i
            double sin_theta = sin(theta); // sin(theta_i)
            double tan_theta = tan(theta); // tan(theta_i)
            
            for (int j = 0; j < n3; j++) {
                // calculate laplace
                double r_lo = (k == 0) ? ghost3D_mu->f0_ib[i][j] : mu[k-1][i][j];
                double r_hi = (k == (n1-1)) ? ghost3D_mu->f0_ob[i][j] : mu[k+1][i][j];
                double the_lo = (i == 0) ? ghost3D_mu->f0_the_lo[k][j] : mu[k][i-1][j];
                double the_hi = (i == (n2-1)) ? ghost3D_mu->f0_the_hi[k][j] : mu[k][i+1][j];
                double phi_lo = (j == 0) ? ghost3D_mu->f0_phi_lo[k][i] : mu[k][i][j-1];
                double phi_hi = (j == (n3-1)) ? ghost3D_mu->f0_phi_hi[k][i] : mu[k][i][j+1];
                // in-plane derivative
                double laplace = (the_hi+the_lo-2.0*mu[k][i][j])/(dTheta*dTheta) + (phi_hi+phi_lo-2.0*mu[k][i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (the_hi-the_lo)/(2.0*tan_theta*dTheta);
                laplace /= (r*r);
                // out-of-plane derivative
                laplace += (r_hi+r_lo-2.0*mu[k][i][j])/(dR*dR) + (r_hi-r_lo)/(r*dR);
                psi[k][i][j] += dt * props.M * laplace;
            }
        }
    }
}

// given average and standard deviation, initialize psi by Guassian Random noise
void InnerSolvPatch::initPsiGuass(const double ave, const double std) {
	double sum = 0.0;
	double x1, x2, w, y1, y2;
	int iset = 0;
    
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                if (iset == 0) {
                    do {
                        x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
                        x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
                        w = x1*x1 + x2*x2;
                    } while ( w >= 1.0 );
                    
                    w = sqrt(-2.0*log(w)/w);
                    y1 = x1 * w;
                    y2 = x2 * w;
                    psi[k][i][j] = y2;
                    iset = 1;
                }
                else {
                    psi[k][i][j] = y1;
                    iset = 0;
                }
                sum = sum + psi[k][i][j];
            }
        }
    }
	
	double average = sum / (n1*n2*n3);
	
	//rescale
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                psi[k][i][j] = ave + (psi[k][i][j]-average)*std;
            }
        }
    }
}

// initialize psi=inVal in sphere with given radius and center(0, 0, 0), psi=otVal outside
void InnerSolvPatch::initPsiSphere(double radius, double inVal, double otVal) {
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    
    const int nrCut = radius > 0 ? grids3D->radiusToIndex(radius) : 0;
    // inside of sphere
    for (int k = 0; k < nrCut; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                psi[k][i][j] = inVal;
            }
        }
    }
    // outside of sphere
    for (int k = nrCut; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                psi[k][i][j] = otVal;
            }
        }
    }
}

// initialize psi=inVal in sphere with given radius and center(x, y, z) in self grid, psi=otVal outside
void InnerSolvPatch::initPsiSprSfGd(double x0, double y0, double z0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                double x = grids3D->x(k, i, j);
                double y = grids3D->y(k, i, j);
                double z = grids3D->z(k, i, j);
                double dist = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
                if (dist < radius) psi[k][i][j] = inVal;
                else psi[k][i][j] = otVal;
            }
        }
    }
}

// initialize psi=inVal in sphere with given radius and center(xp, yp, zp) in opposite grid, psi=otVal outside
void InnerSolvPatch::initPsiSprOpGd(double xp0, double yp0, double zp0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    const int n3 = psi[0][0].size();
    
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                double xp = grids3D->xprime(k, i, j);
                double yp = grids3D->yprime(k, i, j);
                double zp = grids3D->zprime(k, i, j);
                double dist = sqrt((xp-xp0)*(xp-xp0) + (yp-yp0)*(yp-yp0) + (zp-zp0)*(zp-zp0));
                if (dist < radius) psi[k][i][j] = inVal;
                else psi[k][i][j] = otVal;
            }
        }
    }
}







