//
//  GhostBC.h
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#ifndef ____GhostBC__
#define ____GhostBC__

#include <iostream>
#include <vector>
#include <utility>
#include "NormalGrids.h"

#define PI 3.14159265

using namespace std;
class GhostBC {
public:
    GhostBC(const vector<double> &theta, const vector<double> &phi, const NormalGrids &oppoGrids);
    ~GhostBC(){};
    void interpolation(vector<double> &f0, const vector<vector<double> > &f); // given f, interpolate f0

private:
    void coordiTrans(); // calcualte theta', phi' in complemental coordinate
    void findNbs(const NormalGrids &grids); // from given complemental normal grids, find neighbors for interpolation
    void calWts(const NormalGrids &grids); // from the neighbors' position, calculate interpolation weights
    
    vector<pair<double, double> > angles; // (theta, phi) position of ghost boundary
    vector<pair<double, double> > anglesPrim; // (theta', phi') position of ghost boundary in complemental coordinate
    vector<vector<pair<int, int> > > nbs; // index of neighbors for ghost cell in complemental coordinates
    vector<vector<double> > weights; // weights for each neighbor for interpolation
    static const int numNbs = 4;
};




#endif /* defined(____GhostBC__) */
