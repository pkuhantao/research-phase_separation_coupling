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
    double M; // mobility
    double w; // coefficients in phase field energy
    double a;
    double b;
    double Lambda; // coupling strength
    SolvProp() {M=1; w=1; a=1; b=1; Lambda=0;}; // by default, no coupling between inner solvent and membrane
    SolvProp(double M1, double w1, double a1, double b1, double Lambda1) : M(M1), w(w1), a(a1), b(b1), Lambda(Lambda1) {};
};

// membrane properties
struct MembProp {
    double M; // mobility
    double w; // coefficients in phase field energy
    double a;
    double b;
    double Lambda; // coupling strength
    double thickness; // thickness of the membrane
    MembProp() {M=1; w=1; a=1; b=1; Lambda=0; thickness=1;}; // by default, no coupling between inner solvent and membrane
    MembProp(double M1, double w1, double a1, double b1, double Lambda1, double thickness1) : M(M1), w(w1), a(a1), b(b1), Lambda(Lambda1), thickness(thickness1) {};
};



#endif /* defined(____Properties__) */
