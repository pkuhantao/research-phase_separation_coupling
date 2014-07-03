//
//  NormalGrids3D.h
//  
//
//  Created by Tao Han on 4/30/14.
//
//

// normal grids for the inner solvent

#ifndef ____NormalGrids3D__
#define ____NormalGrids3D__

#include <iostream>
#include <vector>

#define PI 3.14159265


using namespace std;

class NormalGrids3D {
public:
    // pass in outer radius and number of cells in r, theta and phi direction
    NormalGrids3D(double ro, int NR, int NTheta, int NPhi);
    ~NormalGrids3D(){};
    int sizeInR() const {return rVal.size();}; // number of points on r axis
    int sizeInTheta() const {return thetaVal.size();}; // number of points on theta axis
    int sizeInPhi() const {return phiVal.size();}; // number of points on phi axis
    
    double dRAD() const {return dR;};
    double dTHETA() const {return dTheta;};
    double dPHI() const {return dPhi;};
    
    double rad(int k) const;
    double theta(int i) const;
    double phi(int i) const;
    
    double x(int k, int i, int j) const;
    double y(int k, int i, int j) const;
    double z(int k, int i, int j) const;
    double xprime(int k, int i, int j) const;
    double yprime(int k, int i, int j) const;
    double zprime(int k, int i, int j) const;
    
    int radiusToIndex(double r) const; // given radius, return the index for that radius(round to closest)
    int thetaToIndex(double angle) const; // given theta angle, return the index for that angle(round downwards)
    int phiToIndex(double angle) const; // given phi angle, return the index for that angle(round downwards)
    
private:
    static const double thetaMin = PI/4.0;
    static const double thetaMax = 3*PI/4.0;
    static const double phiMin = -3*PI/4.0;
    static const double phiMax = 3*PI/4.0;
    
    double ro_; // outer radius limit
    double dR;
    double dTheta;
    double dPhi;
    
    vector<double> rVal;
    vector<double> thetaVal;
    vector<double> phiVal;
    
    vector<vector<vector<double> > > x0, y0, z0; // (x, y, z) in its own coordinates
    vector<vector<vector<double> > > x1, y1, z1; // (x, y, z) in its opposite coordinates
};


#endif /* defined(____NormalGrids3D__) */
