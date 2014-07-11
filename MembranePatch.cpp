//
//  MembranePatch.cpp
//  
//
//  Created by Tao Han on 7/8/14.
//
//

#include "MembranePatch.h"
#include "Properties.h"
#include "NormalGrids2D.h"
#include "GhostFrame2D.h"

#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cassert>

using namespace std;

// constructor
MembranePatch::MembranePatch(double dt, MembProp membProp, NormalGrids2D *cur2D, NormalGrids2D *opp2D): grids2D(cur2D) {
    this->dt = dt;
    props = membProp;
    // create the ghost frame for psi and mu
    ghost2D_psi = new GhostFrame2D(*cur2D, *opp2D);
    ghost2D_mu = new GhostFrame2D(*cur2D, *opp2D);
    // create psi and mu matrix
    const int nTheta = cur2D->sizeInTheta();
    const int nPhi = cur2D->sizeInPhi();
    psi.resize(nTheta, vector<double>(nPhi, 0));
    mu.resize(nTheta, vector<double>(nPhi, 0));
}

MembranePatch::~MembranePatch() {
    delete ghost2D_psi;
    delete ghost2D_mu;
}

// given psi in opposite grids, interpolate the ghost cells for psi in current grids
void MembranePatch::interpPsi(const vector<vector<double> > &oppPsi) {
    ghost2D_psi->interpolate(oppPsi);
}

// given mu in opposite grids, interpolate the ghost cells for mu in current grids
void MembranePatch::interpMu(const vector<vector<double> > &oppMu) {
    ghost2D_mu->interpolate(oppMu);
}

// calculate mu from double well potential
void MembranePatch::calcMu_dw() {
    const double dTheta = grids2D->dTHETA(); // steps
    const double dPhi = grids2D->dPHI();
    const double r = grids2D->rad(); // radius
    const int nrow = mu.size(); // dimension
    const int ncol = mu[0].size();
    
    for (int i = 0; i < nrow; i++) {
        double theta = grids2D->theta(i); // theta_i
        double sin_theta = sin(theta); // sin(theta_i)
        double tan_theta = tan(theta); // tan(theta_i)
        
        for (int j = 0; j < ncol; j++) {
            // calculate laplace
            double lt = (i == 0) ? ghost2D_psi->f0_the_lo_[j] : psi[i-1][j];
            double rt = (i == (nrow-1)) ? ghost2D_psi->f0_the_hi_[j] : psi[i+1][j];
            double dn = (j == 0) ? ghost2D_psi->f0_phi_lo_[i] : psi[i][j-1];
            double up = (j == (ncol-1)) ? ghost2D_psi->f0_phi_hi_[i] : psi[i][j+1];
            double laplace = (rt+lt-2.0*psi[i][j])/(dTheta*dTheta) + (up+dn-2.0*psi[i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (rt-lt)/(2.0*tan_theta*dTheta);
            laplace /= (r*r);
            
            mu[i][j] = -0.5*props.w*props.w*laplace - props.a*psi[i][j] + props.b*pow(psi[i][j], 3);
        }
    }
}

// calculate mu from simple diffusion model, where mu=psi
void MembranePatch::calcMu_sd() {
    for (int i = 0; i < mu.size(); i++) {
        for (int j = 0; j < mu[0].size(); j++) {
            mu[i][j] = psi[i][j];
        }
    }
}

// update psi from mu
void MembranePatch::updatePsi() {
    const double dTheta = grids2D->dTHETA(); // steps
    const double dPhi = grids2D->dPHI();
    const double r = grids2D->rad(); // radius
    const int nrow = psi.size(); // dimension
    const int ncol = psi[0].size();
    
    for (int i = 0; i < nrow; i++) {
        double theta = grids2D->theta(i); // theta_i
        double sin_theta = sin(theta); // sin(theta_i)
        double tan_theta = tan(theta); // tan(theta_i)
        
        for (int j = 0; j < ncol; j++) {
            // calculate laplace
            double lt = (i == 0) ? ghost2D_mu->f0_the_lo_[j] : mu[i-1][j];
            double rt = (i == (nrow-1)) ? ghost2D_mu->f0_the_hi_[j] : mu[i+1][j];
            double dn = (j == 0) ? ghost2D_mu->f0_phi_lo_[i] : mu[i][j-1];
            double up = (j == (ncol-1)) ? ghost2D_mu->f0_phi_hi_[i] : mu[i][j+1];
            double laplace = (rt+lt-2.0*mu[i][j])/(dTheta*dTheta) + (up+dn-2.0*mu[i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (rt-lt)/(2.0*tan_theta*dTheta);
            laplace /= (r*r);
            psi[i][j] += dt * props.M * laplace;
        }
    }
}

// initialize order parameter by Gaussian Distribution
void MembranePatch::initPsiGuass(const double ave, const double std) {
	double sum = 0.0;
	double x1, x2, w, y1, y2;
	int iset = 0;
    
    const int nrow = psi.size();
    const int ncol = psi[0].size();
	
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			if (iset == 0) {
				do {
					x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
					x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
					w = x1*x1 + x2*x2;
				} while ( w >= 1.0 );
				
				w = sqrt(-2.0*log(w)/w);
				y1 = x1 * w;
				y2 = x2 * w;
				psi[i][j] = y2;
				iset = 1;
			}
			else {
				psi[i][j] = y1;
				iset = 0;
			}
			sum = sum + psi[i][j];
		}
	}
	
	double average = sum / (nrow * ncol);
	
	//rescale
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			psi[i][j] = ave + (psi[i][j] - average) * std;
		}
	}
}

// initialize order parameter by constant
void MembranePatch::initPsiConst(const double val) {
    for (int i = 0; i < psi.size(); i++) {
        for (int j = 0; j < psi[0].size(); j++) {
            psi[i][j] = val;
        }
    }
}

// initialize psi=inVal within sphere with given radius and center(x, y, z) in self grid, psi=otVal outside, effectively form a disk on the membrane
void MembranePatch::initPsiDiskSfGd(double x0, double y0, double z0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            double x = grids2D->x(i, j);
            double y = grids2D->y(i, j);
            double z = grids2D->z(i, j);
            double dist = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
            if (dist < radius) psi[i][j] = inVal;
            else psi[i][j] = otVal;
        }
    }
}

// initialize psi=inVal within sphere with given radius and center(xp, yp, zp) in opposite grid, psi=otVal outside, effectively form a disk on the membrane
void MembranePatch::initPsiDiskOpGd(double xp0, double yp0, double zp0, double radius, double inVal, double otVal) {
    assert(radius >= 0);
    // dimensions
    const int n1 = psi.size();
    const int n2 = psi[0].size();
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            double xp = grids2D->xprime(i, j);
            double yp = grids2D->yprime(i, j);
            double zp = grids2D->zprime(i, j);
            double dist = sqrt((xp-xp0)*(xp-xp0) + (yp-yp0)*(yp-yp0) + (zp-zp0)*(zp-zp0));
            if (dist < radius) psi[i][j] = inVal;
            else psi[i][j] = otVal;
        }
    }
}








