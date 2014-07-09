//
//  solv_SD.cpp
//  
//
//  Created by Tao Han on 5/14/14.
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
    
    const string datafolder = "./SD_data";
    
    srand(time(NULL)); // initialize random seed
    
    InnerSolv inSolv(dt); // create the inner solvent
    // initialization
    double x0 = 1;
    double y0 = 1;
    double z0 = 1;
    double radius = 1.5;
    inSolv.initPsiSprYinGd(x0, y0, z0, radius, -1, 1);
    
    // evolve the dynamics
    for (int i = 0; i <= Nstep; i++) {
        inSolv.calcMu_sd();
        
        // store the state every NprintState steps
		if ((i % NprintState) == 0) {
            printMergeSolvState(inSolv, i, datafolder);
        }

        inSolv.updatePsi();
    }
    
    return 0;
}

