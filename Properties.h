//
//  Properties.h
//  
//
//  Created by Tao Han on 5/13/14.
//
//

#ifndef ____Properties__
#define ____Properties__


struct SolvProp {
    double M;
    double w;
    double a;
    double b;
    SolvProp() {M=0; w=0; a=0; b=0;};
    SolvProp(double M1, double w1, double a1, double b1) : M(M1), w(w1), a(a1), b(b1) {};
};




#endif /* defined(____Properties__) */
