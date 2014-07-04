//
//  NormalGrids2D.h
//  
//
//  Created by Tao Han on 4/12/14.
//
//

#ifndef ____NormalGrids2D__
#define ____NormalGrids2D__

#include <iostream>
#include <vector>

#define PI 3.14159265


using namespace std;

class NormalGrids2D {
public:
    // pass in number of cells in theta and phi direction
    NormalGrids2D(int NTheta, int NPhi, double r);
    ~NormalGrids2D(){};
    int sizeInTheta() const {return thetaVal.size();}; // number of points on theta axis
    int sizeInPhi() const {return phiVal.size();}; // number of points on phi axis
    double dTHETA() const {return dTheta;};
    double dPHI() const {return dPhi;};
    double rad() const {return radius;};
    double theta(int i) const;
    double phi(int i) const;
    double x(int i, int j) const;
    double y(int i, int j) const;
    double z(int i, int j) const;
    double xprime(int i, int j) const;
    double yprime(int i, int j) const;
    double zprime(int i, int j) const;
    
    int thetaToIndex(double angle) const; // given theta angle, return the index for that angle(round downwards)
    int phiToIndex(double angle) const; // given phi angle, return the index for that angle(round downwards)
    
private:
    static const double thetaMin = PI/4.0;
    static const double thetaMax = 3*PI/4.0;
    static const double phiMin = -3*PI/4.0;
    static const double phiMax = 3*PI/4.0;
    
    double dTheta;
    double dPhi;
    double radius;
    vector<double> thetaVal;
    vector<double> phiVal;
    vector<vector<double> > x0, y0, z0; // (x, y, z) in its own coordinates
    vector<vector<double> > x1, y1, z1; // (x, y, z) in its opposite coordinates
};


#endif /* defined(____NormalGrids2D__) */
