
#ifndef _H_ITkrQueryClustersTool
#define _H_ITkrQueryCLustersTool

#include "GaudiKernel/IAlgTool.h"



// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrQueryClustersTool("ITkrQueryClustersTool", 1 , 0); 

/** @class IReconTool
* @brief Abstract interface for methods to query TkrClusters
*
* @author Leon Rochester
* $Header$
*
* Example of usage:
*
* #include "GaudiKernel/IToolSvc.h"
* #include "GaudiKernel/AlgTool.h"
* #include "TkrUtil/ITkrQueryClustersTool.h"
*
*...
*
* // private data member of Algorithm or Service
* IToolSvc* m_pToolSvc;
*
* // in initialize
* m_pToolSvc = 0;
* sc = service("ToolSvc", m_pToolSvc, true);
* if (!sc.isSuccess ()){
*     log << MSG::INFO << "Can't find ToolSvc, will quit now" << endreq;
*     return StatusCode::FAILURE;
* }
*
* // in execute
*
* ITkrQueryClustersTool* pQuery;
*       StatusCode sc = m_pToolSvc->retrieveTool("TkrMeritTool", pQuery);
*        if( sc.isFailure() ) {
*            log << MSG::ERROR << "Unable to find a TkrQueryClustersTool" << endreq;
*        }
*...
*
* int nHits = pQuery->numberOfHitsNear(view, layer, inDistance, x0);
*
*/

class   ITkrQueryClustersTool : virtual public IAlgTool {
public:
 
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrQueryClustersTool; }

        /// returns the mean space point in for a given view and layer
    virtual Point meanHit(Event::TkrCluster::view v, int layer) = 0;
    /** returns the mean space point for a given layer, view, within 
    * "inDistance" of a point Pini in the measurement view, and within 
    * one tower in the other view.
    */
    virtual Point meanHitInside (Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini) = 0;
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within one tower in the other view, and a ref. 
    * to the id
    */
    virtual Point nearestHitOutside(Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini, int& id) = 0;
    
    /// Finds the number of clusters with measured distances 
    /// inside a square of side 2*inDistance of a point
    virtual int numberOfHitsNear( int layer, double inDistance, const Point& x0) = 0;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point
    virtual int numberOfHitsNear( int layer, double dX, double dY, const Point& x0) = 0;
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within one tower.
    virtual int numberOfHitsNear( Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& x0) = 0;

};

#endif  // _H_ITkrQueryClustersTool
