//
//  util_funcs.cpp
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#include "util_funcs.h"
#include "InnerSolv.h"
#include "InnerSolvPatch.h"
#include "Center.h"
#include "Membrane.h"
#include "MembranePatch.h"
#include "NormalGrids2D.h"
#include "NormalGrids3D.h"


#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>

#define PI 3.14159265

using namespace std;


// print out the whole inner solvent state in cartesian coordinates (x, y, z) of Yin grids
void printSolvState_car(const InnerSolv &inSolv, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/solv_YinYang_psi_" << stepnum << ".dat";
	
	FILE *fp_psi;
	
	//create files
	if ((fp_psi = fopen(oss.str().c_str(), "w")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	//output data
	fprintf(fp_psi, "VARIABLES = x, y, z, psi\n");
    // get all parts
    const InnerSolvPatch *Yin = inSolv.YinPart();
    const InnerSolvPatch *Yang = inSolv.YangPart();
    const Center *ct = inSolv.CenterPart();
    const NormalGrids3D *Yin3D = Yin->grids3D;
    const NormalGrids3D *Yang3D = Yang->grids3D;
    
    // dimensions (assume Yin and Yang have the same dimensions)
    const int n1 = Yin->psi.size();
    const int n2 = Yin->psi[0].size();
    const int n3 = Yin->psi[0][0].size();
    
    // print out the center part
    fprintf(fp_psi, "0 0 0 %lf\n", ct->psiAtCenter());
    // print out the Yin part
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                fprintf(fp_psi, "%lf %lf %lf %lf\n", Yin3D->x(k, i, j), Yin3D->y(k, i, j), Yin3D->z(k, i, j), Yin->psi[k][i][j]);
            }
        }
    }
    // print out the Yang part
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                fprintf(fp_psi, "%lf %lf %lf %lf\n", Yang3D->xprime(k, i, j), Yang3D->yprime(k, i, j), Yang3D->zprime(k, i, j), Yang->psi[k][i][j]);
            }
        }
    }
    
    //close the file
	fclose(fp_psi);
}



// print out the membrane Yin/Yang patch state in its own spherical coordinates
void printMembPatchState_sph(const vector<vector<double> > &psi, const NormalGrids2D &grids, bool isYinPatch, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/memb_" << (isYinPatch ? "Yin" : "Yang") << "_psi_" << stepnum << ".dat";
	
	FILE *fp_psi;
	
	//create files
	if ((fp_psi = fopen(oss.str().c_str(), "w")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	//output data
	fprintf(fp_psi, "VARIABLES = phi, r, theta, psi\n");
	fprintf(fp_psi, "ZONE I = %d, J = %d, K = %d\n", (int)psi[0].size(), 1, (int)psi.size());
    
    // print out the result
	double radius = grids.rad();
	for (int i = 0; i < psi.size(); i++) {
        double theta = grids.theta(i);
        for (int j = 0; j < psi[0].size(); j++) {
            double phi = grids.phi(j);
            fprintf(fp_psi, "%lf %lf %lf %lf\n", phi, radius, theta, psi[i][j]);
        }
    }
	//close the file
	fclose(fp_psi);
}



// print out the whole membrane state in cartesian coordinates (x, y, z) of Yin grids
void printMembState_car(const Membrane &memb, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/memb_YinYang_psi_" << stepnum << ".dat";
	
	FILE *fp_psi;
	
	//create files
	if ((fp_psi = fopen(oss.str().c_str(), "w")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	//output data
	fprintf(fp_psi, "VARIABLES = x, y, z, psi\n");
    
    // get all parts
    const MembranePatch *Yin = memb.YinPart();
    const MembranePatch *Yang = memb.YangPart();
    const NormalGrids2D *Yin2D = memb.YinGrids();
    const NormalGrids2D *Yang2D = memb.YangGrids();
    
    // dimensions (assume Yin and Yang have the same dimensions)
    const int n1 = Yin->psi.size();
    const int n2 = Yin->psi[0].size();
    
    
    // print out the Yin part
	for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            fprintf(fp_psi, "%lf %lf %lf %lf\n", Yin2D->x(i, j), Yin2D->y(i, j), Yin2D->z(i, j), Yin->psi[i][j]);
        }
    }
    // print out the Yang part
	for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            fprintf(fp_psi, "%lf %lf %lf %lf\n", Yang2D->xprime(i, j), Yang2D->yprime(i, j), Yang2D->zprime(i, j), Yang->psi[i][j]);
        }
    }
    
	//close the file
	fclose(fp_psi);
}

// print out the analysis for the whole 2D membrane
void printMembAnaly(const Membrane &memb, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/memb_analysis.dat";
	
	FILE *fp_analy;
	
	// open the file
	if ((fp_analy = fopen(oss.str().c_str(), "a")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	// add title if stepnum == 0
	if (stepnum == 0) fprintf(fp_analy, "VARIABLES = timestep, max_psi, min_psi, ave_psi\n");
    // output analysis
    fprintf(fp_analy, "%d %lf %lf %lf\n", stepnum, max_psi_memb(memb), min_psi_memb(memb), ave_psi_memb(memb));
    
    // close the file
	fclose(fp_analy);
}

// print out the analysis for the whole 3D inner solvent
void printInSolvAnaly(const InnerSolv &inSolv, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/inSolv_analysis.dat";
	
	FILE *fp_analy;
	
	// open the file
	if ((fp_analy = fopen(oss.str().c_str(), "a")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	// add title if stepnum == 0
	if (stepnum == 0) fprintf(fp_analy, "VARIABLES = timestep, max_psi, min_psi, ave_psi\n");
    // output analysis
    fprintf(fp_analy, "%d %lf %lf %lf\n", stepnum, max_psi_inso(inSolv), min_psi_inso(inSolv), ave_psi_inso(inSolv));
    
    // close the file
	fclose(fp_analy);
}


// max of given 2D matrix
double max_2D(const vector<vector<double> > &mat) {
    const int n1 = mat.size();
    const int n2 = mat[0].size();
    
    double max = mat[0][0];
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            if (mat[i][j] > max) max = mat[i][j];
        }
    }
    
    return max;
}

// min of given 2D matrix
double min_2D(const vector<vector<double> > &mat) {
    const int n1 = mat.size();
    const int n2 = mat[0].size();
    
    double min = mat[0][0];
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            if (mat[i][j] < min) min = mat[i][j];
        }
    }
    
    return min;
}

// max of given 3D matrix
double max_3D(const vector<vector<vector<double> > > &mat) {
    const int n1 = mat.size();
    const int n2 = mat[0].size();
    const int n3 = mat[0][0].size();
    
    double max = mat[0][0][0];
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                if (mat[k][i][j] > max) max = mat[k][i][j];
            }
        }
    }
    
    return max;
}

// min of given 3D matrix
double min_3D(const vector<vector<vector<double> > > &mat) {
    const int n1 = mat.size();
    const int n2 = mat[0].size();
    const int n3 = mat[0][0].size();
    
    double min = mat[0][0][0];
    for (int k = 0; k < n1; k++) {
        for (int i = 0; i < n2; i++) {
            for (int j = 0; j < n3; j++) {
                if (mat[k][i][j] < min) min = mat[k][i][j];
            }
        }
    }
    
    return min;
}

// max of psi on 2D membrane
double max_psi_memb(const Membrane &memb) {
    double max_Yin = max_2D(memb.YinPart()->psi);
    double max_Yang = max_2D(memb.YangPart()->psi);
    return max(max_Yin, max_Yang);
}

// min of psi on 2D membrane
double min_psi_memb(const Membrane &memb) {
    double min_Yin = min_2D(memb.YinPart()->psi);
    double min_Yang = min_2D(memb.YangPart()->psi);
    return min(min_Yin, min_Yang);
}

// average of psi on 2D membrane
double ave_psi_memb(const Membrane &memb) {
    const MembranePatch *Yin = memb.YinPart();
    const MembranePatch *Yang = memb.YangPart();
    const NormalGrids2D *YinGrids = memb.YinGrids();
    const NormalGrids2D *YangGrids = memb.YangGrids();
    
    double sum_psi = 0.0;  // sum of weighted psi, since each cell has different size
    double sum_area = 0.0; // sum of rescaled area
    const int n1 = Yin->psi.size(); // dimensions
    const int n2 = Yin->psi[0].size();
    // consider Yin patch
    for (int i = 0; i < n1; i++) {
        double wt = sin(YinGrids->theta(i)); // weight for each cell
        for (int j = 0; j < n2; j++) {
            sum_psi += Yin->psi[i][j]*wt;
            sum_area += wt;
        }
    }
    // consider Yang patch, assuming Yin and Yang have the same dimensions
    for (int i = 0; i < n1; i++) {
        double wt = sin(YangGrids->theta(i)); // weight for each cell
        for (int j = 0; j < n2; j++) {
            sum_psi += Yang->psi[i][j]*wt;
            sum_area += wt;
        }
    }
    
    return sum_psi/sum_area;
}


// max of psi in 3D inner solvent
double max_psi_inso(const InnerSolv &inSolv) {
    double max_Yin = max_3D(inSolv.YinPart()->psi);
    double max_Yang = max_3D(inSolv.YangPart()->psi);
    double max_Ct = inSolv.CenterPart()->psiAtCenter();
    return max(max(max_Yin, max_Yang), max_Ct);
}

// min of psi in 3D inner solvent
double min_psi_inso(const InnerSolv &inSolv) {
    double min_Yin = min_3D(inSolv.YinPart()->psi);
    double min_Yang = min_3D(inSolv.YangPart()->psi);
    double min_Ct = inSolv.CenterPart()->psiAtCenter();
    return min(min(min_Yin, min_Yang), min_Ct);
}

// average of psi in 3D inner solvent
// Note: there exists overlapping cells from Yin&Yang patches
// Note: the radius of the center part is ri-dr/2
double ave_psi_inso(const InnerSolv &inSolv) {
    const InnerSolvPatch *Yin = inSolv.YinPart();
    const InnerSolvPatch *Yang = inSolv.YangPart();
    const Center *ct = inSolv.CenterPart();
    const NormalGrids3D *YinGrids = Yin->grids3D;
    const NormalGrids3D *YangGrids = Yang->grids3D;
    
    // calculate the Yin part
    double sum_psi_Yin = 0.0;  // sum of weighted psi, since each cell has different size
    double sum_vol_Yin = 0.0; // sum of volume
    int n1 = Yin->psi.size(); // dimensions
    int n2 = Yin->psi[0].size();
    int n3 = Yin->psi[0][0].size();
    for (int k = 0; k < n1; k++) {
        double rad = YinGrids->rad(k); // radius
        for (int i = 0; i < n2; i++) {
            double wt = rad * rad * sin(YinGrids->theta(i)); // weight for each cell
            for (int j = 0; j < n3; j++) {
                sum_psi_Yin += Yin->psi[k][i][j]*wt;
                sum_vol_Yin += wt;
            }
        }
    }
    sum_psi_Yin *= YinGrids->dRAD() * YinGrids->dTHETA() * YinGrids->dPHI();
    sum_vol_Yin *= YinGrids->dRAD() * YinGrids->dTHETA() * YinGrids->dPHI();
    
    // calculate the Yang part
    double sum_psi_Yang = 0.0;  // sum of weighted psi, since each cell has different size
    double sum_vol_Yang = 0.0; // sum of volume
    n1 = Yang->psi.size(); // dimensions
    n2 = Yang->psi[0].size();
    n3 = Yang->psi[0][0].size();
    for (int k = 0; k < n1; k++) {
        double rad = YangGrids->rad(k); // radius
        for (int i = 0; i < n2; i++) {
            double wt = rad * rad * sin(YangGrids->theta(i)); // weight for each cell
            for (int j = 0; j < n3; j++) {
                sum_psi_Yang += Yang->psi[k][i][j]*wt;
                sum_vol_Yang += wt;
            }
        }
    }
    sum_psi_Yang *= YangGrids->dRAD() * YangGrids->dTHETA() * YangGrids->dPHI();
    sum_vol_Yang *= YangGrids->dRAD() * YangGrids->dTHETA() * YangGrids->dPHI();
    
    // calculate the center part
    double rad_ct = YinGrids->rad(0) - YinGrids->dRAD()/2; // rad = ri-dr/2
    double sum_vol_ct = (4.0*PI/3.0) * rad_ct * rad_ct * rad_ct;
    double sum_psi_ct = ct->psiAtCenter() * sum_vol_ct;
    
    // return the average
    return (sum_psi_Yin+sum_psi_Yang+sum_psi_ct) / (sum_vol_Yin+sum_vol_Yang+sum_vol_ct);
}














