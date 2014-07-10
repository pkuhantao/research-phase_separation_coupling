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
#include "Membrane.h"
#include "NormalGrids2D.h"

#include <vector>
#include <string>

using namespace std;

// print out the whole inner solvent state in cartesian coordinates (x, y, z) of Yin grids
void printSolvState_car(const InnerSolv &inSolv, int stepnum, string foldername);

// print out the membrane Yin/Yang patch state in its own spherical coordinates
void printMembPatchState_sph(const vector<vector<double> > &psi, const NormalGrids2D &grids, bool isYinPatch, int stepnum, string foldername);

// print out the whole membrane state in cartesian coordinates (x, y, z) of Yin grids
void printMembState_car(const Membrane &memb, int stepnum, string foldername);



#endif /* defined(____util_funcs__) */
