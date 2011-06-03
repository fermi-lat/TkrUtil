/** @file Doca.h
@brief header file for Doca.cxx (q.v.)
@author Tracy Usher, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Doca.h,v 1.5 2006/03/06 02:32:43 lsrea Exp $
*/


#ifndef Doca_H
#define Doca_H

#include "geometry/Ray.h"

/** @class Doca

  @brief does some calculations about the DOCA between two Rays 
  
  A utility for calculating the distance of closest approach between two
  (assumed straight) tracks and the point on each line where this occurs.

  Also can do the distance between a line and a point
    
  The method is taken from "Distance between Lines and Segments with their
  Closest Point of Approach" found at 
  http://www.geometryalgorithms.com/Archive/algorithm_0106/algorithm_0106.htm
            
 @author Tracy Usher 03/05/02          
*/
namespace {
    const Point  p_nullRay(0., 0., 0.);
    const Vector v_nullRay(0., 0., 0.);
    const Ray    nullRay(p_nullRay, v_nullRay);
}
class Doca
{
public:
    Doca(const Ray& ray1, const Ray& ray2=nullRay);
    Doca(const Point& point1, const Vector& vector1,
        const Point& point2=p_nullRay, const Vector& vector2=v_nullRay);
    ~Doca() {}
    
    /// return doca between two Rays 
    double docaRay1Ray2()   {return doca;}
    /// return distance of doca from origin of 1st Ray
    double arcLenRay1()     {return s;}
    /// return distance of doca from origin of 2nd Ray
    double arcLenRay2()     {return t;}
    /// point of doca on 1st Ray
    Point  docaPointRay1();
    /// point of doca on 2nd Ray
    Point  docaPointRay2();
    /// doca of a point to the first line
    double docaOfPoint(const Point& p);
    
private:
    // stuff common to both constructors
    void ini();
    enum mode {LINELINE, LINEPOINT};
    /// origin of 1st Ray
    Point  P;
    /// direction of 1st Ray
    Vector u;
    /// origin of 2nd Ray
    Point  Q;
    /// direction of 2nd Ray
    Vector v;
    /// DOCA between the two Rays
    double doca;
    /// distance of DOCA point from origin 1
    double s;
    /// distance of DOCA point from origin 2
    double t;
    /// mode of doca call, either 2 lines or a line and a point
    mode m_mode;
    /// storage for point in LINEPOINT mode
    Point p1;
};

#endif
