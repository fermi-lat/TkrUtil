/** @file Doca.cxx
@brief calculates the DOCA between two straight lines, and related quantities
@author Tracy Usher, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Doca.cxx,v 1.4 2006/03/06 02:32:43 lsrea Exp $
*/

#include "Doca.h"

//
// ini() is where all the real work gets done...
//


Doca::Doca(const Ray& ray1, 
           const Ray& ray2) 
           : P(ray1.position()), u(ray1.direction()), 
           Q(ray2.position()), v(ray2.direction()), 
           p1(p_nullRay)
{
    ini();
}

Doca::Doca(const Point& point1, const Vector& vector1,
           const Point& point2, const Vector& vector2)
           : P(point1), u(vector1), Q(point2), v(vector2), 
           p1(p_nullRay)
{
    ini();
}

void Doca::ini() {

    u = u.unit();

    m_mode = LINELINE;
    if(v.mag()==0.0) {
        m_mode = LINEPOINT;
        return;
    }

    v = v.unit();

    //Determine vector from start point track 1 to start of track 2
    Vector w     = P - Q;

    //Projections of of tracks along vector between start points
    double d     = u.dot(w);
    double e     = v.dot(w);

    //Dot product between two tracks to check if parallel
    double b     = u.dot(v);
    double denom = 1. - b*b;

    //Lines are not parallel
    if (fabs(b) < 1.)
    {
        s    = (b*e - d  ) / denom;
        t    = (e   - b*d) / denom;
        w    = w + s * u - t * v;
        doca = w.magnitude();
    }
    //Lines are parallel
    else
    {
        s    = 0;
        t    = d / b;
        w    = w - t * v;
        doca = w.magnitude();
    }
}


Point Doca::docaPointRay1()
{
     return P + s * u;
}

Point Doca::docaPointRay2()
{
    if(m_mode==LINEPOINT) return p_nullRay;

    return Q + t * v;
}

double Doca::docaOfPoint(const Point& p)
{
    p1 = p;
    Vector dp =  p1 - P;
    s  = dp.dot(u);

    return (dp.cross(u)).mag();
}
