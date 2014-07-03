//
//  PhaseSeparation_solv.cpp
//  
//
//  Created by Tao Han on 5/21/14.
//
//

#include "InnerSolv.h"
#include "util_funcs.h"

#include <time.h>
#include <stdlib.h>

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Wrong input: program + dt + Nstep!" << endl;
        exit(1);
    }
    
    const double dt = atof(argv[1]);
    const int Nstep = atoi(argv[2]);
    const int NprintState = Nstep/20;
    const int NprintAnaly = Nstep/20;
    
    const string datafolder = "./PS_data";
    
    srand(time(NULL)); // initialize random seed
    
    InnerSolv inSolv(dt); // create the inner solvent
    // initialization
    double ave = 0.2;
    double std = 0.1;
    inSolv.initPsiGauss(ave, std);
    
    // evolve the dynamics
    for (int i = 0; i <= Nstep; i++) {
        inSolv.calcMu_dw();
        
        // store the state every NprintState steps
		if ((i % NprintState) == 0) {
            printMergeSolvState(inSolv, i, datafolder);
        }

        inSolv.updatePsi();
    }
    
    return 0;
}

