//
//  GhostBC2D.h
//  
//
//  Created by Tao Han on 4/12/14.
//
//

// 2D ghost Boundary (given one boundary, create its interpolation scheme from the opposite Yin/Yang grids)

#ifndef ____GhostBC2D__
#define ____GhostBC2D__

#include "NormalGrids2D.h"

#include <vector>
#include <utility>


using namespace std;

class GhostBC2D {
public:
    // given the boundary definition(r=const, theta, phi) in its own grids and the opposite Yin/Yang grids,
    // construct the interpolation scheme in this boundary
    GhostBC2D(const vector<double> &theta, const vector<double> &phi, const NormalGrids2D &oppoGrids);
    ~GhostBC2D(){};
    void interpolation(vector<double> &f0, const vector<vector<double> > &f) const; // given value f in opposite grids, interpolate f0 in this boundary in its own grids

private:
    void coordiTrans(); // calcualte theta', phi' in opposite/complemental Yin/Yang grids
    void findNbs(const NormalGrids2D &grids); // from given complemental normal grids, find neighbors for interpolation
    void calWts(const NormalGrids2D &grids); // from the neighbors' position, calculate interpolation weights
    
    vector<pair<double, double> > angles; // (theta, phi) position of ghost boundary
    vector<pair<double, double> > anglesPrim; // (theta', phi') position of ghost boundary in complemental coordinate
    vector<vector<pair<int, int> > > nbs; // index of neighbors for ghost cell in complemental coordinates
    vector<vector<double> > weights; // weights for each neighbor for interpolation
    static const int numNbs = 4;
};




#endif /* defined(____GhostBC2D__) */
