//
//  memb_SD.cpp
//  
//
//  Created by Tao Han on 7/11/14.
//
//

// simple diffusion in 2D membrane only

#include "Properties.h"
#include "Membrane.h"
#include "util_funcs.h"


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <tr1/unordered_map>
#include <string>
#include <sstream>
#include <utility>
#include <ctime>

using namespace std;


int main(int argc, char* argv[]) {
    if (argc != 2) {
		cerr << "Wrong inputs! Should be (./$program $input_filename)" << endl;
		exit(1);
	}
    
    ifstream inputfile(argv[1]); // input file
	tr1::unordered_map<string, string> paras; // parameters' map
    
    // open the file
	if (inputfile.is_open()) {
		// obtain all parameters and their value, and store into the map
		while (inputfile.good()) {
			string oneline;  // one pair of parameter and its value
			getline(inputfile, oneline); // get one line from the file
			if (oneline.empty()) continue;
			istringstream iss(oneline);
			string key, val;
			iss >> key >> val;
			paras.insert(make_pair<string, string>(key, val)); // add new pair into map
		}
		inputfile.close(); // close the parameters file
	}
	else {
		cerr << "Unable to open file" << endl;
		exit(1);
	}
    
    //parse parameters value
	const string datafolder = paras["datafolder"]; // folder name for storing data
	const int nthreads = atoi(paras["nthreads"].c_str()); // number of openmp threads
    // space and time info
    const double radius = atof(paras["radius"].c_str()); // radius
    const int nR = atoi(paras["nR"].c_str()); // space dimensions
    const int nTheta = atoi(paras["nTheta"].c_str());
    const int nPhi = atoi(paras["nPhi"].c_str());
    const int nt = atoi(paras["nt"].c_str()); // number of time steps
    const int t_start = atoi(paras["t_start"].c_str()); // start time
    const double dt = atof(paras["dt"].c_str()); // time step
    // parameters for data storage and analysis
	const int numConfig = atoi(paras["numConfig"].c_str()); // number of configurations
	const int numAnaly = atoi(paras["numAnaly"].c_str());  // nmber of analysis
	const int nPrintState = (nt - 1)/numConfig; // every nPrintState steps, store the state
	const int nPrintAnaly = (nt - 1)/numAnaly; // every nPrintAnaly steps, store the analysis
    
    // parameters for membrane
    MembProp memb_Props;
    memb_Props.M = atof(paras["M_m"].c_str()); // membrane's mobility
    memb_Props.w = atof(paras["w_m"].c_str()); // membrane's coefficients in phase field energy
    memb_Props.a = atof(paras["a_m"].c_str());
    memb_Props.b = atof(paras["b_m"].c_str());
    memb_Props.Lambda = atof(paras["Lambda"].c_str()); // coupling strength
    memb_Props.thickness = atof(paras["thickness"].c_str()); // membrane thickness
    memb_Props.radius = radius; // radius
    memb_Props.nTheta = nTheta;
    memb_Props.nPhi = nPhi;
    

    // create the membrane
    Membrane memb(dt, memb_Props);
    // initialize psi by a disk
    const double x0 = radius/2.0;
    const double y0 = 0;
    const double z0 = 0;
    const double nu_rad = radius*2.0/3.0;
    memb.initPsiDiskYinGd(x0, y0, z0, nu_rad, -1, 1);
    
    
    // evolve the dynamics
    for (int i = 0; i < nt; i++) {
        // calculate chemical potential by simple diffusion model, where mu=psi
        memb.calcMu_sd();
        
        // store the state every nPrintState steps
		if ((i % nPrintState) == 0) {
			printMembPatchesState_sfSph(memb, i+t_start, datafolder);
		}
        // store the analysis every nPrintAnaly steps
        if ((i % nPrintAnaly) == 0) {
			printMembAnaly(memb, i+t_start, datafolder);
		}
        
        // update order parameters
        memb.updatePsi();
    }
    
    return 0;
}



