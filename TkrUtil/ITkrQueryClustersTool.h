/** @file ITkrQueryClustersTool.h
* @brief Abstract interface for methods to query TkrClusters

 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrQueryClustersTool.h,v 1.14 2005/01/30 07:12:28 lsrea Exp $
*/


#ifndef _H_ITkrQueryClustersTool
#define _H_ITkrQueryClustersTool

#include "GaudiKernel/IAlgTool.h"

#include "geometry/Point.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrQueryClustersTool("ITkrQueryClustersTool", 3 , 1); 

/** @class ITkrQueryClustersTool
* @brief Abstract interface for methods to query TkrClusters
*/

class   ITkrQueryClustersTool : virtual public IAlgTool {
public:
 
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrQueryClustersTool; }

    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within one tower in the other view, and a ref. 
    * to the id
    */
    /// real cluster
    virtual Point nearestHitOutside(int v, int layer, 
        double inDistance, const Point& Pini, int& id) const = 0;
    /// bad cluster
    virtual Point nearestBadHitOutside(int v, int layer, 
        double inDistance, const Point& Pini, int& id) const = 0;
    /** returns the nearest cluster found outside of "inDistance" of a point "Pini"
    * in the measured view, within one tower in the other view
    */
    /// real cluster
    virtual Event::TkrCluster* nearestClusterOutside(int v, int layer, 
                               double inDistance, const Point& Pini) const = 0;
    /// bad cluster
    virtual Event::TkrCluster* nearestBadClusterOutside(int v, int layer, 
                               double inDistance, const Point& Pini) const = 0;
    
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point, in requested bilayer
    virtual int numberOfHitsNear( int layer, double dX, double dY, 
        const Point& x0, const Vector dir = Vector(0., 0., 1.)) const = 0;

    /// Finds the number of unused clusters within 2*dX by 2*dY of a point 
    /// and within one tower, in requested layer and view
    virtual int numberOfUUHitsNear( int layer, double dX, double dY, 
        const Point& x0, const Vector dir = Vector(0., 0., 1.)) const = 0;
    virtual int numberOfHitsNear( int v, int layer, double inDistance, 
        const Point& x0, const Vector dir = Vector(0., 0., 1.)) const = 0;

    /// Access clusters by view and layer or by TkrId
    virtual const Event::TkrClusterVec  getClustersReverseLayer(int view, int layer) const = 0;
    virtual const Event::TkrClusterVec  getClusters(int view, int layer) const = 0;
    virtual const Event::TkrClusterVec& getClusters(const idents::TkrId& tkrId) const = 0;
    /// bad clusters
    virtual const Event::TkrClusterVec  getBadClusters(int view, int layer) const = 0;
    virtual const Event::TkrClusterVec& getBadClusters(const idents::TkrId& tkrId) const = 0;

    virtual double clusterWidth(Event::TkrCluster* cluster) const = 0;
};

#endif  // _H_ITkrQueryClustersTool
