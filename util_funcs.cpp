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
#include "NormalGrids2D.h"
#include "NormalGrids3D.h"


#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <stdlib.h>
#include <stdio.h>


using namespace std;


// print out the whole inner solvent state in (x, y, z) of Yin grids
void printMergeSolvState(const InnerSolv &inSolv, int stepnum, string foldername) {
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



//print out memb order parameter
void printMembState(const vector<vector<double> > &psi, const NormalGrids2D &grids, int YinYangID, int stepnum, string foldername) {
    ostringstream oss;
	oss << foldername << "/memb_" << (YinYangID==0 ? "Yin" : "Yang") << "_psi_" << stepnum << ".dat";
	
	FILE *fp_psi;
	
	//create files
	if ((fp_psi = fopen(oss.str().c_str(), "w")) == NULL) {
		printf("cannot open file\n");
		exit(1);
	}
	
	//output data
	fprintf(fp_psi, "VARIABLES = phi, r, theta, psi\n");
	fprintf(fp_psi, "ZONE I = %d, J = %d, K = %d\n", psi[0].size(), 1, psi.size());
    
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



// print out the whole membrane in (x, y, z) of Yin grids
void printMergeMembState(const vector<vector<double> > &psi_Yin, const NormalGrids2D &grids_Yin, const vector<vector<double> > &psi_Yang, const NormalGrids2D &grids_Yang, int stepnum, string foldername) {
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
    
    // print out the Yin part
	for (int i = 0; i < psi_Yin.size(); i++) {
        for (int j = 0; j < psi_Yin[0].size(); j++) {
            fprintf(fp_psi, "%lf %lf %lf %lf\n", grids_Yin.x(i, j), grids_Yin.y(i, j), grids_Yin.z(i, j), psi_Yin[i][j]);
        }
    }
    // print out the Yang part
	for (int i = 0; i < psi_Yang.size(); i++) {
        for (int j = 0; j < psi_Yang[0].size(); j++) {
            fprintf(fp_psi, "%lf %lf %lf %lf\n", grids_Yang.xprime(i, j), grids_Yang.yprime(i, j), grids_Yang.zprime(i, j), psi_Yang[i][j]);
        }
    }
    
	//close the file
	fclose(fp_psi);
}





