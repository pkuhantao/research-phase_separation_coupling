//
//  Membrane.cpp
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#include "Membrane.h"
#include "NormalGrids2D.h"
#include "GhostFrame2D.h"

#include <vector>
#include <time.h>
#include <cstdlib>
#include <cmath>

using namespace std;

Membrane::Membrane(double M, double dt, NormalGrids2D *grids, GhostFrame2D *frame) {
    this->M = M;
    this->dt = dt;
    normGrids = grids;
    ghostFm = frame;
    
    const int nrow = normGrids->sizeInTheta();
    const int ncol = normGrids->sizeInPhi();
    // regular part
    psi.resize(nrow, vector<double>(ncol, 0));
    mu.resize(nrow, vector<double>(ncol, 0));
    // ghost part
    psi_the_lo.resize(ncol, 0);
    psi_the_hi.resize(ncol, 0);
    psi_phi_lo.resize(nrow, 0);
    psi_phi_hi.resize(nrow, 0);
    
    mu_the_lo.resize(ncol, 0);
    mu_the_hi.resize(ncol, 0);
    mu_phi_lo.resize(nrow, 0);
    mu_phi_hi.resize(nrow, 0);
}

// initialize order parameter by Gaussian Distribution
void Membrane::initGauRand(double ave, double std) {
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

// Note: this function is just used for a rough testing, it is not accurate at all, since
// the step in i and j may not be equal
// initialize order parameter by disk(+1) surrounded by -1
void Membrane::initDisk(int ic, int jc, int rad) {
    for (int i = 0; i < psi.size(); i++) {
		for (int j = 0; j < psi[0].size(); j++) {
			double dist = sqrt((double)(i-ic)*(i-ic) + (double)(j-jc)*(j-jc));
            if (dist <= rad) psi[i][j] = 1;
            else psi[i][j] = -1;
		}
	}
}


// interpolate psi in ghost cells by given psi in complemental grids
void Membrane::interpPsi(const vector<vector<double> > &oppoPsi) {
    ghostFm->interpolate(psi_the_lo, psi_the_hi, psi_phi_lo, psi_phi_hi, oppoPsi);
}

// interpolate mu in ghost cells by given mu in complemental
void Membrane::interpMu(const vector<vector<double> > &oppoMu) {
    ghostFm->interpolate(mu_the_lo, mu_the_hi, mu_phi_lo, mu_phi_hi, oppoMu);
}

// calculate mu from psi
// simple diffusion
void Membrane::calcMu_simdiff() {
    for (int i = 0; i < mu.size(); i++) {
        for (int j = 0; j < mu[0].size(); j++) {
            mu[i][j] = psi[i][j];
        }
    }
}

// calculate mu from double well potential
void Membrane::calcMu_doublewell() {
    const double dTheta = normGrids->dTHETA();
    const double dPhi = normGrids->dPHI();
    const double r = normGrids->rad();
    const double nrow = mu.size();
    const double ncol = mu[0].size();
    
    for (int i = 0; i < nrow; i++) {
        double theta = normGrids->theta(i); // theta_i
        double sin_theta = sin(theta); // sin(theta_i)
        double tan_theta = tan(theta); // tan(theta_i)
        
        for (int j = 0; j < ncol; j++) {
            // calculate laplace
            double lt = (i == 0) ? psi_the_lo[j] : psi[i-1][j];
            double rt = (i == (nrow-1)) ? psi_the_hi[j] : psi[i+1][j];
            double dn = (j == 0) ? psi_phi_lo[i] : psi[i][j-1];
            double up = (j == (ncol-1)) ? psi_phi_hi[i] : psi[i][j+1];
            double laplace = (rt+lt-2.0*psi[i][j])/(dTheta*dTheta) + (up+dn-2.0*psi[i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (rt-lt)/(2.0*tan_theta*dTheta);
            laplace /= (r*r);
            
            mu[i][j] = -0.5*w*w*laplace - a*psi[i][j] + b*pow(psi[i][j], 3);
        }
    }
}



// update psi from mu
void Membrane::updatePsi() {
    const double dTheta = normGrids->dTHETA();
    const double dPhi = normGrids->dPHI();
    const double r = normGrids->rad();
    const double nrow = psi.size();
    const double ncol = psi[0].size();
    
    for (int i = 0; i < nrow; i++) {
        double theta = normGrids->theta(i); // theta_i
        double sin_theta = sin(theta); // sin(theta_i)
        double tan_theta = tan(theta); // tan(theta_i)
        
        for (int j = 0; j < ncol; j++) {
            // calculate laplace
            double lt = (i == 0) ? mu_the_lo[j] : mu[i-1][j];
            double rt = (i == (nrow-1)) ? mu_the_hi[j] : mu[i+1][j];
            double dn = (j == 0) ? mu_phi_lo[i] : mu[i][j-1];
            double up = (j == (ncol-1)) ? mu_phi_hi[i] : mu[i][j+1];
            double laplace = (rt+lt-2.0*mu[i][j])/(dTheta*dTheta) + (up+dn-2.0*mu[i][j])/(sin_theta*sin_theta*dPhi*dPhi) + (rt-lt)/(2.0*tan_theta*dTheta);
            laplace /= (r*r);
            psi[i][j] += dt * M * laplace;
        }
    }
}
















