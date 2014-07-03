//
//  util_funcs.h
//  
//
//  Created by Tao Han on 4/13/14.
//
//

#ifndef ____util_funcs__
#define ____util_funcs__

#include "InnerSolv.h"
#include "NormalGrids.h"

#include <iostream>
#include <vector>

using namespace std;

// print out the whole inner solvent state in (x, y, z) of Yin grids
void printMergeSolvState(const InnerSolv &inSolv, int stepnum, string foldername);

//print out memb order parameter in its own spherical coordinates
void printMembState(const vector<vector<double> > &psi, const NormalGrids &grids, int YinYangID, int stepnum, string foldername);
// print out the whole membrane in (x, y, z) of Yin grids
void printMergeMembState(const vector<vector<double> > &psi_Yin, const NormalGrids &grids_Yin, const vector<vector<double> > &psi_Yang, const NormalGrids &grids_Yang, int stepnum, string foldername);



#endif /* defined(____util_funcs__) */
