// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrQueryClustersTool.cxx,v 1.4 2004/09/18 18:38:42 usher Exp $

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

#include <vector>
#include <map>
#include "geometry/Point.h"  

namespace 
{
    // fraction of towerPitch within which to accept a hit.  At 0.55, only points within
    // a couple of mm of the edge of an adjacent tower will be accepted, so this number
    // is probably too small.

    const double _towerFactor = 0.55;
}

typedef std::pair<int,int>         TkrViewLayerPair;
typedef std::vector<idents::TkrId> TkrIdVector;

struct CompareViewLayer
{
public:
    bool operator()(const TkrViewLayerPair left, const TkrViewLayerPair right) const
    {
        int leftPlane  = left.first  + 2 * left.second;
        int rightPlane = right.first + 2 * right.second;

        return leftPlane < rightPlane;
    }
};

typedef std::multimap<TkrViewLayerPair,idents::TkrId,CompareViewLayer> TkrViewLayerIdMap;

class TkrQueryClustersTool : public AlgTool, virtual public ITkrQueryClustersTool 
{
public:
    
    TkrQueryClustersTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrQueryClustersTool() { }

    StatusCode initialize();
    
    /// returns the mean space point for a given view and layer
    Point meanHit(int v, int layer) const;
    /** returns the mean space point for a given layer, view, within 
    * "inDistance" of a point Pini in the measurement view, and within 
    * "one tower" in the other view.
    */
    Point meanHitInside (int v, int layer, 
        double inDistance, const Point& Pini) const;
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within "one tower" in the other view, and a ref. 
    * to the id
    */
    Point nearestHitOutside(int v, int layer, 
        double inDistance, const Point& Pini, int& id) const;
    
    /// Finds the number of clusters with measured distances 
    /// inside a square of side 2*inDistance of a point
    int numberOfHitsNear( int layer, double inDistance, const Point& x0) const;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point
    int numberOfHitsNear( int layer, double dX, double dY, const Point& x0) const;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point which are not used
    int numberOfUUHitsNear( int layer, double dX, double dY, const Point& x0) const;
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within "one tower."
    int numberOfHitsNear( int v, int layer, double inDistance, const Point& x0) const;

    const Event::TkrClusterVec  getClustersReverseLayer(int view, int layer) const;
    const Event::TkrClusterVec  getClusters(int view, int layer) const;
    const Event::TkrClusterVec& getClusters(const idents::TkrId& tkrId) const;
    
private:

    /// Checks that a layer number is in the correct range, and sets some variables
    bool validLayer(int layer) const
    {
        m_idClusMap = SmartDataPtr<Event::TkrIdClusterMap>(m_pEventSvc, 
        EventModel::TkrRecon::TkrIdClusterMap);

        // check for valid layer
        return (layer>=0 && layer < m_pGeom->numLayers());
    };

    void initIdMap() const;

    // some pointers to services
    
    /// pointer to tracker geometry
    ITkrGeometrySvc*  m_pGeom;
    /// pointer to event data service
    IDataProviderSvc* m_pEventSvc;
    /// save test distance
    double m_testDistance;
    /// save pointer to clusters
    mutable Event::TkrIdClusterMap* m_idClusMap;

    /// THE table of life
    mutable TkrViewLayerIdMap m_ViewLayerIdMap;
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

    //m_pClus     = 0;
    m_idClusMap = 0;
    m_ViewLayerIdMap.clear();
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
        // test distance (in unmeasured view)
        m_testDistance = _towerFactor*m_pGeom->towerPitch();

        sc = serviceLocator()->service( "EventDataSvc", m_pEventSvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find EventSvc" << endreq;
            return sc;
        }
        // pointer to clusters
    }
    log << MSG::INFO << "TkrQueryClustersTool successfully initialized" << endreq;
    return sc;
}

void TkrQueryClustersTool::initIdMap() const
{
    // This will build the multi map converting view,layer pairs to TkrIds
    // Loop over views, view = 0 is x, view = 1 is y
    for(int view = 0; view < m_pGeom->numViews(); view++)
    {
        // Loop over layers, layer = 0 is at the bottom/back
        for(int layer = 0; layer < m_pGeom->numLayers(); layer++)
        {
            TkrViewLayerPair viewLayerPair(view,layer);
            int tray   = 0;
            int botTop = 0;

            // Convert to tray/bottom/top
            m_pGeom->layerToTray(layer, view, tray, botTop);

            // Two sets of loops over the towers
            for(int towerX = 0; towerX < m_pGeom->numXTowers(); towerX++)
            {
                for(int towerY = 0; towerY < m_pGeom->numYTowers(); towerY++)
                {
                    idents::TkrId tkrId(towerX, towerY, tray, botTop == 1, view);
                    m_ViewLayerIdMap.insert(std::pair<TkrViewLayerPair,idents::TkrId>(viewLayerPair,tkrId));
                }
            }
        }
    }
}


const Event::TkrClusterVec TkrQueryClustersTool::getClustersReverseLayer(int view, int reverseLayer) const
{
    int layer = m_pGeom->numLayers() - reverseLayer - 1;

    return getClusters(view, layer);
}

const Event::TkrClusterVec TkrQueryClustersTool::getClusters(int view, int layer) const
{
    Event::TkrClusterVec clusVec;

    if (!validLayer(layer)) return clusVec;

    if (m_ViewLayerIdMap.size() == 0) initIdMap();

    TkrViewLayerPair viewLayerPair(view,layer);

    std::pair<TkrViewLayerIdMap::const_iterator,TkrViewLayerIdMap::const_iterator> 
        clusIdRange = m_ViewLayerIdMap.equal_range(viewLayerPair);
    int numIds  = m_ViewLayerIdMap.count(viewLayerPair);

    for(TkrViewLayerIdMap::const_iterator clusIdIter = clusIdRange.first; clusIdIter != clusIdRange.second; clusIdIter++)
    {
        const idents::TkrId& newId = (*clusIdIter).second;

        const Event::TkrClusterVec& newClus = getClusters(newId);

        clusVec.insert(clusVec.end(),newClus.begin(),newClus.end());
    }

    return clusVec;
}
    
const Event::TkrClusterVec& TkrQueryClustersTool::getClusters(const idents::TkrId& tkrId) const
{
    return (*m_idClusMap)[tkrId];
}

Point TkrQueryClustersTool::meanHit(int view, int layer) const
{
    // Purpose and Method: Returns the mean position of all clusters in a 
    //       layer
    // Inputs:  view and layer number
    // Outputs:  mean position of all the clusters in the layer
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point Pini(0.,0.,0);
    
    if (!validLayer(layer)) return Pini;

    const Event::TkrClusterVec clusters = getClustersReverseLayer(view, layer);
    Point Pini3(0.,0.,0.);

    for(Event::TkrClusterVecConItr clusIter = clusters.begin(); clusIter != clusters.end(); clusIter++)
    {
        const Event::TkrCluster* cluster = *clusIter;
        Pini += cluster->position();
    }

    int nHits = clusters.size();
    Point Pini2(Pini.x()/nHits,Pini.y()/nHits,Pini.z()/nHits);

    return Pini2;
}

Point TkrQueryClustersTool::meanHitInside(int view, int layer, double inDistance, 
                                          const Point& Pcenter) const
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

    const Event::TkrClusterVec clusters = getClustersReverseLayer(view, layer);
    int nhits = clusters.size();
    if (nhits == 0) return P;
    
    double nsum = 0.;
    double xsum = 0.;
    double ysum = 0.;
    double zsum = 0.;
    double hitDistance, twrDistance;
    
    for(Event::TkrClusterVecConItr clusIter = clusters.begin(); clusIter != clusters.end(); clusIter++)
    {
        P = (*clusIter)->position();
        
        if        (view == idents::TkrId::eMeasureX) {
            hitDistance = fabs(P.x() - Pcenter.x());
            twrDistance = fabs(P.y() - Pcenter.y());
        } else if (view == idents::TkrId::eMeasureY) {
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

Point TkrQueryClustersTool::nearestHitOutside(int view, int layer, double inDistance, 
                                              const Point& Pcenter, int& id) const
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

    const Event::TkrClusterVec clusters = getClustersReverseLayer(view, layer);
    int nhits = clusters.size();
    if (nhits == 0) return Pnear;
    
    double minDistance = inDistance;
    double maxDistance = 1e6;
    Point Pini(0.,0.,0.);
    for(Event::TkrClusterVecConItr clusIter = clusters.begin(); clusIter != clusters.end(); clusIter++)
    {
        const Event::TkrCluster* cluster = (*clusIter);

        if (cluster->hitFlagged()) continue;
        
        Pini = cluster->position();
        
        // Kludge to prevent crashes when z layer incorrect
        //double zDistance   = fabs(Pini.z() - Pcenter.z());
        //if (zDistance > .3) continue;
        
        double hitDistance = fabs(Pini.x() - Pcenter.x());
        double twrDistance = fabs(Pini.y() - Pcenter.y());
        
        if      (view == idents::TkrId::eMeasureY) 
        {
            hitDistance = fabs(Pini.y() - Pcenter.y());
            twrDistance = fabs(Pini.x() - Pcenter.x());
        }
        else if (view != idents::TkrId::eMeasureX) 
        {
            hitDistance = (Pini-Pcenter).mag();
            twrDistance = 0.;
        }
        
        if ( hitDistance >= minDistance && hitDistance < maxDistance 
                                        && twrDistance < m_testDistance) 
        {
            maxDistance = hitDistance;
            Pnear     = Pini;
            id        = cluster->id();
        }
    }
    return Pnear;
}

int TkrQueryClustersTool::numberOfHitsNear(int layer, double inDistance, 
                                       const Point& x0) const
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
                                       const Point& x0) const
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
    int view = idents::TkrId::eMeasureX;
    const Event::TkrClusterVec xList = getClustersReverseLayer(view, layer);
    int nHitsInPlane = xList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - xList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - xList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < dX && fabs(hitDiffY) < m_testDistance) numHits++;
    }
    
    // Look for hits in the Y view of desired layer
    view = idents::TkrId::eMeasureY;
    const Event::TkrClusterVec yList = getClustersReverseLayer(view, layer);
    nHitsInPlane = yList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - yList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - yList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < m_testDistance && fabs(hitDiffY) < dY) numHits++;
    }
    
    return numHits;
}

int TkrQueryClustersTool::numberOfUUHitsNear( int layer, double dX, double dY, 
                                       const Point& x0) const
{
    // Purpose and Method: counts the number of un-used hits in a bilayer 
    //      within a rectangle of sides 2*dX, 2*dY
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numHits = 0;
    
    if (!validLayer(layer)) return numHits;

    //Look for hits in the X view of desired layer
    const Event::TkrClusterVec xList = getClustersReverseLayer(idents::TkrId::eMeasureX, layer);
    int nHitsInPlane = xList.size();
    
    while(nHitsInPlane--)
    {
        if(xList[nHitsInPlane]->hitFlagged()) continue; 

        double hitDiffX = x0.x() - xList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - xList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < dX && fabs(hitDiffY) < m_testDistance) numHits++;
    }
    
    // Look for hits in the Y view of desired layer
    const Event::TkrClusterVec yList = getClustersReverseLayer(idents::TkrId::eMeasureY, layer);
    nHitsInPlane = yList.size();
    
    while(nHitsInPlane--)
    {
        if(yList[nHitsInPlane]->hitFlagged()) continue;

        double hitDiffX = x0.x() - yList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - yList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < m_testDistance && fabs(hitDiffY) < dY) numHits++;
    }
    
    return numHits;
}

int TkrQueryClustersTool::numberOfHitsNear( int view, int layer, 
                                       double inDistance, const Point& x0) const
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
    const Event::TkrClusterVec clusters = getClustersReverseLayer(view, layer);
    
    for(Event::TkrClusterVecConItr clusIter = clusters.begin(); clusIter != clusters.end(); clusIter++)
    {
        const Event::TkrCluster* cluster = *clusIter;
        double hitDiffV = view == idents::TkrId::eMeasureX 
            ? x0.x() - cluster->position().x()
            : x0.y() - cluster->position().y();
        double hitDiffO = view == idents::TkrId::eMeasureY 
            ? x0.y() - cluster->position().y()
            : x0.x() - cluster->position().x();
        
        if (fabs(hitDiffV) < inDistance && fabs(hitDiffO) < m_testDistance) 
            numHits++;
    }
    
    return numHits;
}

