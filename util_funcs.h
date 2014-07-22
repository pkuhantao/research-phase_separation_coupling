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
void printInSolvState_car(const InnerSolv &inSolv, int stepnum, string foldername);

// print out the inner solvent Yin&Yang patches' state separately in their own spherical coordinates
// and the center part's state in cartesian coordinates
void printInSolvPatchesState_sfSph(const InnerSolv &inSolv, int stepnum, string foldername);

// print out the analysis for the whole 3D inner solvent
void printInSolvAnaly(const InnerSolv &inSolv, int stepnum, string foldername);

// print out the whole membrane state in cartesian coordinates (x, y, z) of Yin grids
void printMembState_car(const Membrane &memb, int stepnum, string foldername);

// print out the membrane Yin&Yang patches' state separately in their own spherical coordinates
void printMembPatchesState_sfSph(const Membrane &memb, int stepnum, string foldername);

// print out the analysis for the whole 2D membrane
void printMembAnaly(const Membrane &memb, int stepnum, string foldername);



// the followings are auxillary functions

// print out the inner solvent Yin/Yang patch state in its own spherical coordinates
void printInsolvPatchState_sph(const InnerSolvPatch &inSolvPatch, bool isYinPatch, int stepnum, string foldername);

// print out the inner solvent center part's state
void printInsolvCtState(const Center &ct, int stepnum, string foldername);

// print out the membrane Yin/Yang patch state in its own spherical coordinates
void printMembPatchState_sph(const MembranePatch &membPatch, bool isYinPatch, int stepnum, string foldername);



// max of given 2D matrix
double max_2D(const vector<vector<double> > &mat);

// min of given 2D matrix
double min_2D(const vector<vector<double> > &mat);

// max of given 3D matrix
double max_3D(const vector<vector<vector<double> > > &mat);

// min of given 3D matrix
double min_3D(const vector<vector<vector<double> > > &mat);

// max of psi on 2D membrane
double max_psi_memb(const Membrane &memb);

// min of psi on 2D membrane
double min_psi_memb(const Membrane &memb);

// average of psi on 2D membrane
// Note: there exists overlapping cells from Yin&Yang patches
double ave_psi_memb(const Membrane &memb);

// max of psi in 3D inner solvent
double max_psi_inso(const InnerSolv &inSolv);

// min of psi in 3D inner solvent
double min_psi_inso(const InnerSolv &inSolv);

// average of psi in 3D inner solvent
// Note: there exists overlapping cells from Yin&Yang patches
// Note: the radius of the center part is ri-dr/2
double ave_psi_inso(const InnerSolv &inSolv);










#endif /* defined(____util_funcs__) */
