//
//  Properties.h
//  
//
//  Created by Tao Han on 5/13/14.
//
//

#ifndef ____Properties__
#define ____Properties__


// inner solvent properties
struct SolvProp {
    // physical properties
    double M = 1; // mobility
    double w = 1; // coefficients in phase field energy
    double a = 1;
    double b = 1;
    double Lambda = 0; // coupling strength, by default, no coupling between inner solvent and membrane
    
    // geometric properties
    double ro = 2; // outer radius
    double ri = 1; // inner radius
    int nR = 1; // number of cells in r direction of Yin/Yang patch
    int nTheta = 1; // number of cells in theta direction of Yin/Yang patch
    int nPhi = 1; // number of cells in phi direction of Yin/Yang patch
};

// membrane properties
struct MembProp {
    // physical properties
    double M = 1; // mobility
    double w = 1; // coefficients in phase field energy
    double a = 1;
    double b = 1;
    double Lambda = 0; // coupling strength, by default, no coupling between inner solvent and membrane
    
    // geometric properties
    double thickness = 1; // thickness of the membrane
    double radius = 2; // curvature radius of the membrane shell
    int nTheta = 1; // number of cells in theta direction of Yin/Yang patch
    int nPhi = 1; // number of cells in phi direction of Yin/Yang patch
};






#endif /* defined(____Properties__) */
