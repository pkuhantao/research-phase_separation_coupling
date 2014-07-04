//
//  GhostBC2D.h
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#ifndef ____GhostBC2D__
#define ____GhostBC2D__

#include "NormalGrids2D.h"

#include <iostream>
#include <vector>
#include <utility>


#define PI 3.14159265

using namespace std;
class GhostBC2D {
public:
    GhostBC2D(const vector<double> &theta, const vector<double> &phi, const NormalGrids2D &oppoGrids);
    ~GhostBC2D(){};
    void interpolation(vector<double> &f0, const vector<vector<double> > &f); // given f, interpolate f0

private:
    void coordiTrans(); // calcualte theta', phi' in complemental coordinate
    void findNbs(const NormalGrids2D &grids); // from given complemental normal grids, find neighbors for interpolation
    void calWts(const NormalGrids2D &grids); // from the neighbors' position, calculate interpolation weights
    
    vector<pair<double, double> > angles; // (theta, phi) position of ghost boundary
    vector<pair<double, double> > anglesPrim; // (theta', phi') position of ghost boundary in complemental coordinate
    vector<vector<pair<int, int> > > nbs; // index of neighbors for ghost cell in complemental coordinates
    vector<vector<double> > weights; // weights for each neighbor for interpolation
    static const int numNbs = 4;
};




#endif /* defined(____GhostBC2D__) */
