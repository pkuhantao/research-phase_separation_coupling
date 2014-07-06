//
//  NormalGrids3D.h
//  
//
//  Created by Tao Han on 4/30/14.
//
//

// normal grids for the 3D Yin/Yang patch

#ifndef ____NormalGrids3D__
#define ____NormalGrids3D__

#include <vector>

#define PI 3.14159265


using namespace std;

class NormalGrids3D {
public:
    // pass in outer and inner radius and number of cells in (r, theta, phi) directions
    NormalGrids3D(double ro, double ri, int NR, int NTheta, int NPhi);
    ~NormalGrids3D(){};
    int sizeInR() const {return rVal.size();}; // number of points on r axis != number of cells
    int sizeInTheta() const {return thetaVal.size();}; // number of points on theta axis
    int sizeInPhi() const {return phiVal.size();}; // number of points on phi axis
    
    double dRAD() const {return dR;};  // steps
    double dTHETA() const {return dTheta;};
    double dPHI() const {return dPhi;};
    
    double rad(int k) const;  // r
    double theta(int i) const;  // theta
    double phi(int j) const;  // phi
    
    double x(int k, int i, int j) const;  // x, y, z in itw own coordinates
    double y(int k, int i, int j) const;
    double z(int k, int i, int j) const;
    
    double xprime(int k, int i, int j) const;  // (x, y, z) in its opposite coordinates
    double yprime(int k, int i, int j) const;
    double zprime(int k, int i, int j) const;
    
    int radiusToIndex(double r) const; // given radius, return the index for that radius(round to closest)
    int thetaToIndex(double angle) const; // given theta angle, return the index for that angle(round downwards)
    int phiToIndex(double angle) const; // given phi angle, return the index for that angle(round downwards)
    
private:
    static const double thetaMin = PI/4.0;  // (theta, phi) range in Yin/Yang grids
    static const double thetaMax = 3*PI/4.0;
    static const double phiMin = -3*PI/4.0;
    static const double phiMax = 3*PI/4.0;
    
    double ro_; // outer radius
    double ri_; // inner radius
    
    double dR;  // steps
    double dTheta;
    double dPhi;
    
    vector<double> rVal;  // (r, theta, phi)
    vector<double> thetaVal;
    vector<double> phiVal;
    
    vector<vector<vector<double> > > x0, y0, z0; // (x, y, z) in its own coordinates
    
    vector<vector<vector<double> > > x1, y1, z1; // (x, y, z) in its opposite coordinates
};


#endif /* defined(____NormalGrids3D__) */
