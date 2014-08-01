//
//  cpl_PS.cpp
//  
//
//  Created by Tao Han on 8/1/14.
//
//

// coupled phase separation in both 2D membrane and 3D inner solvent
// coupled by local coupling effect

#include "Properties.h"
#include "InnerSolv.h"
#include "Membrane.h"
#include "util_funcs.h"


#include <iostream>
#include <cstdlib>
#include <tr1/unordered_map>
#include <string>
#include <ctime>

using namespace std;


int main(int argc, char* argv[]) {
    if (argc != 2) {
		cerr << "Wrong inputs! Should be (./$program $input_filename)" << endl;
		exit(1);
	}
    
    // parse the input file into hashtable at first
    tr1::unordered_map<string, string> paras = parseInputFile(argv[1]);
    // get the inner solvent and membrane properties
    SolvProp solv_Props = getInSolvProp(paras);
    MembProp memb_Props = getMembProp(paras);
    
    // parameters for parallelization
	const int nthreads = atoi(paras["nthreads"].c_str()); // number of openmp threads
    // parameters for time info
    const int nt = atoi(paras["nt"].c_str()); // number of time steps
    const int t_start = atoi(paras["t_start"].c_str()); // start time
    const double dt = atof(paras["dt"].c_str()); // time step
    // parameters for data storage and analysis
    const string datafolder = paras["datafolder"]; // folder name for storing data
	const int numConfig = atoi(paras["numConfig"].c_str()); // number of configurations
	const int numAnaly = atoi(paras["numAnaly"].c_str());  // nmber of analysis
	const int nPrintState = (nt - 1)/numConfig; // every nPrintState steps, store the state
	const int nPrintAnaly = (nt - 1)/numAnaly; // every nPrintAnaly steps, store the analysis
    // initialization info
    const double ave_psi_s = atof(paras["ave_psi_s"].c_str()); // average of psi in inner solvent
    const double std_psi_s = atof(paras["std_psi_s"].c_str()); // standard deviation of psi in inner solvent
    const double ave_psi_m = atof(paras["ave_psi_m"].c_str()); // average of psi on membrane
    const double std_psi_m = atof(paras["std_psi_m"].c_str()); // standard deviation of psi on membrane
    
    // initialize random seed
    srand(time(NULL));
    
    // create inner solvent
    InnerSolv inSol(dt, solv_Props);
    // create membrane
    Membrane memb(dt, memb_Props);
    
    // initialize psi in inner solvent by Guassian Random Distribution
    inSol.initPsiGauss(ave_psi_s, std_psi_s);
    // initialize psi on the membrane by Guassian Random Distribution
    memb.initPsiGuass(ave_psi_m, std_psi_m);
    
    // evolve the dynamics
    for (int i = 0; i < nt; i++) {
        // calculate chemical potential by double-well potential with local coupling
        inSol.calcMu_dw_cp(memb.YinPart(), memb.YangPart());
        memb.calcMu_dw_cp(inSol.YinPart(), inSol.YangPart());
        
        // store the state every nPrintState steps
		if ((i % nPrintState) == 0) {
			printInSolvState_sfSph(inSol, i+t_start, datafolder);
            printMembPatchesState_sfSph(memb, i+t_start, datafolder);
		}
        
        // store the analysis every nPrintAnaly steps
		if ((i % nPrintAnaly) == 0) {
			printInSolvAnaly(inSol, i+t_start, datafolder);
            printMembAnaly(memb, i+t_start, datafolder);
		}
        
        // update order parameters
        inSol.updatePsi();
        memb.updatePsi();
    }
    
    return 0;
}



