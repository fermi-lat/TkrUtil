/** @file ITkrQueryClustersTool.h
* @brief Abstract interface for methods to query TkrClusters

 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrQueryClustersTool.h,v 1.7 2003/04/18 20:27:36 lsrea Exp $
*/


#ifndef _H_ITkrQueryClustersTool
#define _H_ITkrQueryClustersTool

#include "GaudiKernel/IAlgTool.h"

#include "geometry/Point.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrQueryClustersTool("ITkrQueryClustersTool", 2 , 0); 

/** @class ITkrQueryClustersTool
* @brief Abstract interface for methods to query TkrClusters

 Example of usage:

 @verbatim
  #include "GaudiKernel/IToolSvc.h"
  #include "GaudiKernel/AlgTool.h"
  #include "TkrUtil/ITkrQueryClustersTool.h"

...

  // private data member of Algorithm or Service
  IToolSvc m_pToolSvc;

  // in initialize
  m_pToolSvc = 0;
  sc = service("ToolSvc", m_pToolSvc, true);
  if (!sc.isSuccess ()){
      log << MSG::INFO << "Can't find ToolSvc, will quit now" << endreq;
      return StatusCode::FAILURE;
  }

  // in execute

  ITkrQueryClustersTool pQuery;
  StatusCode sc = m_pToolSvc->retrieveTool("TkrMeritTool", pQuery);
      if( sc.isFailure() ) {
          log << MSG::ERROR << "Unable to find a TkrQueryClustersTool" << endreq;
      }
  ...

  int nHits = pQuery->numberOfHitsNear(view, layer, inDistance, x0);
 @endverbatim
*/

class   ITkrQueryClustersTool : virtual public IAlgTool {
public:
 
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrQueryClustersTool; }

        /// returns the mean space point in for a given view and layer
    virtual Point meanHit(Event::TkrCluster::view v, int layer) const = 0;
    /** returns the mean space point for a given layer, view, within 
    * "inDistance" of a point Pini in the measurement view, and within 
    * one tower in the other view.
    */
    virtual Point meanHitInside (Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini) const = 0;
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within one tower in the other view, and a ref. 
    * to the id
    */
    virtual Point nearestHitOutside(Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini, int& id) const = 0;
    
    /// Finds the number of clusters with measured distances 
    /// inside a square of side 2*inDistance of a point, in requested bilayer
    virtual int numberOfHitsNear( int layer, double inDistance, const Point& x0) const = 0;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point, in requested bilayer
    virtual int numberOfHitsNear( int layer, double dX, double dY, const Point& x0) const = 0;
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within one tower, in requested layer and view
    virtual int numberOfHitsNear( Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& x0) const = 0;

    /// Finds the number of unused clusters within 2*dX by 2*dY of a point 
    /// and within one tower, in requested layer and view
    virtual int numberOfUUHitsNear( int layer, double dX, double dY, const Point& x0) const = 0;

};

#endif  // _H_ITkrQueryClustersTool
