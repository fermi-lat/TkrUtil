// $Header$

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

#include <vector>
#include "geometry/Point.h"

namespace 
{
    // fraction of towerPitch within which to accept a hit.  At 0.55, only points within
    // a couple of mm of the edge of an adjacent tower will be accepted, so this number
    // is probably too small.

    const double _towerFactor = 0.55;
}

class TkrQueryClustersTool : public AlgTool, virtual public ITkrQueryClustersTool 
{
public:
    
    TkrQueryClustersTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrQueryClustersTool() { }

    StatusCode initialize();
    
    /// returns the mean space point for a given view and layer
    Point meanHit(Event::TkrCluster::view v, int layer);
    /** returns the mean space point for a given layer, view, within 
    * "inDistance" of a point Pini in the measurement view, and within 
    * "one tower" in the other view.
    */
    Point meanHitInside (Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini) ;
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within "one tower" in the other view, and a ref. 
    * to the id
    */
    Point nearestHitOutside(Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini, int& id);
    
    /// Finds the number of clusters with measured distances 
    /// inside a square of side 2*inDistance of a point
    int numberOfHitsNear( int layer, double inDistance, const Point& x0);
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point
    int numberOfHitsNear( int layer, double dX, double dY, const Point& x0);
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within "one tower."
    int numberOfHitsNear( Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& x0);
    
private:

    /// Checks that a layer number is in the correct range, and sets some variables
    bool validLayer(int layer)
    {
        // pointer to clusters
        m_pClus = SmartDataPtr<Event::TkrClusterCol>(m_pEventSvc, 
        EventModel::TkrRecon::TkrClusterCol);
        // test distance (in unmeasured view)
        m_testDistance = _towerFactor*m_pGeom->towerPitch();
        // check for valid layer
        return (layer>=0 && layer < m_pGeom->numLayers());
    };

    // some pointers to services
    
    /// pointer to tracker geometry
    ITkrGeometrySvc*  m_pGeom;
    /// pointer to event data service
    IDataProviderSvc* m_pEventSvc;
    /// save test distance
    double m_testDistance;
    /// save pointer to clusters
    Event::TkrClusterCol* m_pClus;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrQueryClustersTool> s_factory;
const IToolFactory& TkrQueryClustersToolFactory = s_factory;

// Standard Constructor
TkrQueryClustersTool::TkrQueryClustersTool(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent)
                           : AlgTool( type, name, parent )
{    
    // Declare additional interface
    declareInterface<ITkrQueryClustersTool>(this); 
}

StatusCode TkrQueryClustersTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    if( serviceLocator() ) {
        sc = serviceLocator()->service( "TkrGeometrySvc", m_pGeom, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return sc;
        }

        sc = serviceLocator()->service( "EventDataSvc", m_pEventSvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find EventSvc" << endreq;
            return sc;
        }
    }
    log << MSG::INFO << "TkrQueryClustersTool successfully initialized" << endreq;
    return sc;
}

Point TkrQueryClustersTool::meanHit(Event::TkrCluster::view v, int layer)
{
    // Purpose and Method: Returns the mean position of all clusters in a 
    //       layer
    // Inputs:  view and layer number
    // Outputs:  mean position of all the clusters in the layer
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point Pini(0.,0.,0);
    
    if (!validLayer(layer)) return Pini;
 
    int nhits = m_pClus->nHits(v,layer);
    if (nhits == 0) return Pini;
    
    const std::vector<Event::TkrCluster*>& AuxList = m_pClus->getHits(v,layer);
    for (int ihit=0; ihit<nhits; ihit++){
        Pini += AuxList[ihit]->position();
    }
    Point Pini2(Pini.x()/nhits,Pini.y()/nhits,Pini.z()/nhits);
    return Pini2;
}

Point TkrQueryClustersTool::meanHitInside(Event::TkrCluster::view v, int layer, 
                                      double inDistance, 
                                      const Point& Pcenter)
{
    // Purpose and Method: Returns mean position of hits
    //    within a distance of a point in the measured dimension,
    //    and no more than one tower away
    // Inputs:  view and layer number, Distance and center
    // Outputs:  mean position of clusters satisfying criterion
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point P(0.,0.,0);
    if (!validLayer(layer)) return P;

    const std::vector<Event::TkrCluster*>& AuxList = m_pClus->getHits(v,layer);
    int nhits = AuxList.size();
    if (nhits == 0) return P;
    
    double nsum = 0.;
    double xsum = 0.;
    double ysum = 0.;
    double zsum = 0.;
    double hitDistance, twrDistance;
    
    for (int ihit=0; ihit<nhits; ihit++)
    {
        P = AuxList[ihit]->position();
        
        if        (v == Event::TkrCluster::X) {
            hitDistance = fabs(P.x() - Pcenter.x());
            twrDistance = fabs(P.y() - Pcenter.y());
        } else if (v == Event::TkrCluster::Y) {
            hitDistance = fabs(P.y() - Pcenter.y());
            twrDistance = fabs(P.x() - Pcenter.x());
        } else {
            hitDistance = (P-Pcenter).mag();
            twrDistance = 0.;
        }
        
        // Check that hit is close and within one tower
        if (hitDistance < inDistance && twrDistance < m_testDistance) 
        {
            nsum += 1.;
            xsum += P.x();
            ysum += P.y();
            zsum += P.z();
        }
    }
    
    if (nsum > 0.) P = Point(xsum/nsum, ysum/nsum, zsum/nsum);
   
    return P;
}

Point TkrQueryClustersTool::nearestHitOutside(Event::TkrCluster::view v, 
                                          int layer, double inDistance, 
                                          const Point& Pcenter, int& id)
{
    // Purpose and Method: returns the position of the closest cluster
    //    outside of a given distance from a point in the measured direction,
    //    and in the same or adjacent tower in the other direction.
    // Inputs:  view and layer, center and distance
    // Outputs:  Position of nearest cluster
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point Pnear(0.,0.,0.);
    id = -1;
    
    if (!validLayer(layer)) return Pnear;

    int nhits = m_pClus->nHits(v,layer);
    if (nhits == 0) return Pnear;
    
    const std::vector<Event::TkrCluster*>& AuxList = m_pClus->getHits(v,layer);
    
    double minDistance = inDistance;
    double maxDistance = 1e6;
    Point Pini(0.,0.,0.);
    for (int ihit = 0; ihit< nhits; ihit++) 
    {
        if (AuxList[ihit]->hitFlagged()) continue;
        
        Pini = AuxList[ihit]->position();
        
        // Kludge to prevent crashes when z layer incorrect
        //double zDistance   = fabs(Pini.z() - Pcenter.z());
        //if (zDistance > .3) continue;
        
        double hitDistance = fabs(Pini.x() - Pcenter.x());
        double twrDistance = fabs(Pini.y() - Pcenter.y());
        
        if      (v == Event::TkrCluster::Y) 
        {
            hitDistance = fabs(Pini.y() - Pcenter.y());
            twrDistance = fabs(Pini.x() - Pcenter.x());
        }
        else if (v != Event::TkrCluster::X) 
        {
            hitDistance = (Pini-Pcenter).mag();
            twrDistance = 0.;
        }
        
        if ( hitDistance >= minDistance && hitDistance < maxDistance 
                                        && twrDistance < m_testDistance) 
        {
            maxDistance = hitDistance;
            Pnear     = Pini;
            id        = AuxList[ihit]->id();
        }
    }
    return Pnear;
}

int TkrQueryClustersTool::numberOfHitsNear(int layer, double inDistance, 
                                       const Point& x0)
{
    // Purpose and Method: counts the number of hits in a bilayer 
    //       within a square of side 2*inDistance
    // Inputs:  layer number, distance, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    return numberOfHitsNear(layer, inDistance, inDistance, x0);
}

int TkrQueryClustersTool::numberOfHitsNear( int layer, double dX, double dY, 
                                       const Point& x0)
{
    // Purpose and Method: counts the number of hits in a bilayer 
    //      within a rectangle of sides 2*dX, 2*dY
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numHits = 0;
    
    if (!validLayer(layer)) return numHits;

    //Look for hits in the X view of desired layer
    const std::vector<Event::TkrCluster*>& xList = 
        m_pClus->getHits(Event::TkrCluster::X, layer);
    int nHitsInPlane = xList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - xList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - xList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < dX && fabs(hitDiffY) < m_testDistance) numHits++;
    }
    
    // Look for hits in the Y view of desired layer
    const std::vector<Event::TkrCluster*>& yList = 
        m_pClus->getHits(Event::TkrCluster::Y, layer);
    nHitsInPlane = yList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - yList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - yList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < m_testDistance && fabs(hitDiffY) < dY) numHits++;
    }
    
    return numHits;
}

int TkrQueryClustersTool::numberOfHitsNear( Event::TkrCluster::view v, int layer, 
                                       double inDistance, const Point& x0)
{
    // Purpose and Method: counts the number of hits within a distance 
    //     "inDistance" in the measurement direction, and within one tower 
    //     in the other direction
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numHits = 0;
    
    if (!validLayer(layer)) return numHits;

    // Look for hits in the desired view of the given layer
    const std::vector<Event::TkrCluster*> & auxList = m_pClus->getHits(v, layer);
    int nHitsInPlane = auxList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffV = v == Event::TkrCluster::X 
            ? x0.x() - auxList[nHitsInPlane]->position().x()
            : x0.y() - auxList[nHitsInPlane]->position().y();
        double hitDiffO = v == Event::TkrCluster::X 
            ? x0.y() - auxList[nHitsInPlane]->position().y()
            : x0.x() - auxList[nHitsInPlane]->position().x();
        
        if (fabs(hitDiffV) < inDistance && fabs(hitDiffO) < m_testDistance) 
            numHits++;
    }
    
    return numHits;
}

