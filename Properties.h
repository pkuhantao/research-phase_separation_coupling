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
    double M; // mobility
    double w; // coefficients in phase field energy
    double a;
    double b;
    double Lambda; // coupling strength
    double dPsi = 2.0; // dPsi = psi^+ - psi^-
    
    // geometric properties
    double ro; // outer radius
    double ri; // inner radius
    int nR; // number of cells in r direction of Yin/Yang patch
    int nTheta; // number of cells in theta direction of Yin/Yang patch
    int nPhi; // number of cells in phi direction of Yin/Yang patch
};

// membrane properties
struct MembProp {
    // physical properties
    double M; // mobility
    double w; // coefficients in phase field energy
    double a;
    double b;
    double Lambda; // coupling strength
    double dPsi = 2.0; // dPsi = psi^+ - psi^-
    
    // geometric properties
    double thickness; // thickness of the membrane
    double radius; // curvature radius of the membrane shell
    int nTheta; // number of cells in theta direction of Yin/Yang patch
    int nPhi; // number of cells in phi direction of Yin/Yang patch
};






#endif /* defined(____Properties__) */
