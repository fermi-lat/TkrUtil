// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrQueryClustersTool.cxx,v 1.11 2005/02/11 07:12:54 lsrea Exp $

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"

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
    
    enum clusterType { STANDARDCLUSTERS, BADCLUSTERS };
 
    TkrQueryClustersTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrQueryClustersTool() { }

    StatusCode initialize();
    
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within "one tower" in the other view, and a ref. 
    * to the id
    */
    /// real cluster
    Point nearestHitOutside(int v, int layer, 
        double inDistance, const Point& Pini, int& id) const;
    /// bad cluster
    Point nearestBadHitOutside(int v, int layer, 
        double inDistance, const Point& Pini, int& id) const;
    
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within "one tower" in the other view, and a ref. 
    * to the id
    */
    /// real cluster
    Event::TkrCluster* nearestClusterOutside(int v, int layer, 
                                     double inDistance, const Point& Pini) const;
    /// bad cluster
    Event::TkrCluster* nearestBadClusterOutside(int v, int layer, 
                                     double inDistance, const Point& Pini) const;
    
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point
    int numberOfHitsNear( int layer, double dX, double dY, const Point& x0, 
        const Vector dir) const;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point which are not used
    int numberOfUUHitsNear( int layer, double dX, double dY, const Point& x0, 
        const Vector dir) const;
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within "one tower."
    int numberOfHitsNear( int v, int layer, double inDistance, const Point& x0, 
        const Vector dir) const;

    const Event::TkrClusterVec  getClustersReverseLayer(int view, int layer) const;
    const Event::TkrClusterVec  getClusters(int view, int layer) const;
    const Event::TkrClusterVec& getClusters(const idents::TkrId& tkrId) const;
    const Event::TkrClusterVec  getBadClusters(int view, int layer) const;
    const Event::TkrClusterVec& getBadClusters(const idents::TkrId& tkrId) const;

    /// actual width of a cluster, including gaps
    double clusterWidth(Event::TkrCluster* cluster) const;

private:
    Event::TkrCluster* nearestClusterOutsideX(int v, int layer, 
        double inDistance, const Point& Pini, clusterType type) const;
    Point nearestHitOutsideX(int v, int layer, 
        double inDistance, const Point& Pini, int& id, clusterType type) const;

    int getNumNearHitsInPlane(int layer, int view, double dX,
        const Point& x0, const Vector& dir, bool checkFlag) const;

    int getNumNearHitsInLayer(int layer, double dX, double dY,
        const Point& x0, const Vector& dir, bool checkFlag) const;


    const Event::TkrClusterVec getClustersX(int view, int layer, clusterType type) const;
    const Event::TkrClusterVec& getClustersX(const idents::TkrId& tkrId, clusterType type) const;

    /// Checks that a layer number is in the correct range, and sets some variables
    bool validLayer(int layer, clusterType type=STANDARDCLUSTERS) const;
    void initIdMap() const;

    // some pointers to services
    
    /// pointer to tracker geometry
    ITkrGeometrySvc*  m_tkrGeom;
    /// pointer to event data service
    IDataProviderSvc* m_pEventSvc;
    /// save test distance
    double m_testDistance;
    /// current pointer
    mutable Event::TkrIdClusterMap* m_idClusMap;
    /// same for bad clusters
    mutable Event::TkrIdClusterMap* m_badIdClusMap;

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
        sc = serviceLocator()->service( "TkrGeometrySvc", m_tkrGeom, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return sc;
        }

        ITkrBadStripsSvc* pBadStrips = m_tkrGeom->getTkrBadStripsSvc();
        m_badIdClusMap = 0;
        if(pBadStrips) {
            m_badIdClusMap = pBadStrips->getBadIdClusterMap();
        }

        // test distance (in unmeasured view)
        m_testDistance = _towerFactor*m_tkrGeom->towerPitch();

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
    for(int view = 0; view < m_tkrGeom->numViews(); view++)
    {
        // Loop over layers, layer = 0 is at the bottom/back
        for(int layer = 0; layer < m_tkrGeom->numLayers(); layer++)
        {
            TkrViewLayerPair viewLayerPair(view,layer);
            int tray   = 0;
            int botTop = 0;

            // Convert to tray/bottom/top
            m_tkrGeom->layerToTray(layer, view, tray, botTop);

            // Two sets of loops over the towers
            for(int towerX = 0; towerX < m_tkrGeom->numXTowers(); towerX++)
            {
                for(int towerY = 0; towerY < m_tkrGeom->numYTowers(); towerY++)
                {
                    idents::TkrId tkrId(towerX, towerY, tray, botTop == 1, view);
                    m_ViewLayerIdMap.insert(std::pair<TkrViewLayerPair,idents::TkrId>(viewLayerPair,tkrId));
                }
            }
        }
    }
}

bool TkrQueryClustersTool::validLayer(int layer, clusterType type) const
{
    if(type==STANDARDCLUSTERS) {
        m_idClusMap = SmartDataPtr<Event::TkrIdClusterMap>(m_pEventSvc, 
            EventModel::TkrRecon::TkrIdClusterMap);
    } else {
        m_idClusMap = m_badIdClusMap;
    }

    // check for valid layer
    return (layer>=0 && layer < m_tkrGeom->numLayers());
};

const Event::TkrClusterVec TkrQueryClustersTool::getClustersReverseLayer(
    int view, int reverseLayer) const
{
    int layer = m_tkrGeom->reverseLayerNumber(reverseLayer);

    return getClusters(view, layer);
}

const Event::TkrClusterVec TkrQueryClustersTool::getClusters(
    int view, int layer) const
{
    return getClustersX(view, layer, STANDARDCLUSTERS);
}

const Event::TkrClusterVec TkrQueryClustersTool::getBadClusters(
    int view, int layer) const
{
    return getClustersX(view, layer, BADCLUSTERS);
}

const Event::TkrClusterVec TkrQueryClustersTool::getClustersX(
    int view, int layer, clusterType type) const
{
    Event::TkrClusterVec clusVec;

    if (!validLayer(layer, type)) return clusVec;

    if (m_ViewLayerIdMap.size() == 0) initIdMap();

    TkrViewLayerPair viewLayerPair(view,layer);

    std::pair<TkrViewLayerIdMap::const_iterator,TkrViewLayerIdMap::const_iterator> 
        clusIdRange = m_ViewLayerIdMap.equal_range(viewLayerPair);
    //int numIds  = m_ViewLayerIdMap.count(viewLayerPair);

    // Try to avoid too many calls to the data service
    Event::TkrIdClusterMap* idClusMap = 0;
    if(type==STANDARDCLUSTERS) 
    {
        idClusMap = SmartDataPtr<Event::TkrIdClusterMap>(m_pEventSvc, 
                                                         EventModel::TkrRecon::TkrIdClusterMap);
    }
    else idClusMap = m_badIdClusMap;

    TkrViewLayerIdMap::const_iterator clusIdIter = clusIdRange.first;
    for(; clusIdIter != clusIdRange.second; clusIdIter++)
    {
        const idents::TkrId& newId = (*clusIdIter).second;

        const Event::TkrClusterVec newClus = (*idClusMap)[newId];

        clusVec.insert(clusVec.end(),newClus.begin(),newClus.end());
    }

    return clusVec;
}
    
const Event::TkrClusterVec& TkrQueryClustersTool::getClusters(
    const idents::TkrId& tkrId) const
{
    return getClustersX(tkrId, STANDARDCLUSTERS);
}

const Event::TkrClusterVec& TkrQueryClustersTool::getBadClusters(
    const idents::TkrId& tkrId) const
{
    return getClustersX(tkrId, BADCLUSTERS);
}

const Event::TkrClusterVec& TkrQueryClustersTool::getClustersX(
    const idents::TkrId& tkrId, clusterType type) const
{
     if(type==STANDARDCLUSTERS) {
        m_idClusMap = SmartDataPtr<Event::TkrIdClusterMap>(m_pEventSvc, 
            EventModel::TkrRecon::TkrIdClusterMap);
    return (*m_idClusMap)[tkrId];
    } else {
        return (*m_badIdClusMap)[tkrId];
    }
}

Point TkrQueryClustersTool::nearestHitOutside(
    int view, int layer, double inDistance, 
    const Point& Pcenter, int& id) const
{
    return nearestHitOutsideX(view, layer, inDistance, 
        Pcenter, id, STANDARDCLUSTERS);
}

Point TkrQueryClustersTool::nearestBadHitOutside(
    int view, int layer, double inDistance, 
    const Point& Pcenter, int& id) const
{
    return nearestHitOutsideX(view, layer, inDistance, Pcenter, id, BADCLUSTERS);
}

Point TkrQueryClustersTool::nearestHitOutsideX(
    int view, int layer, double inDistance, 
    const Point& Pcenter, int& id, clusterType type) const
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
    Event::TkrCluster* clus = nearestClusterOutsideX(
        view, layer, inDistance, Pcenter, type);
    if (clus==0) return Pnear;
    
    Pnear = clus->position();
    id    = 0; //clus->id();
    return Pnear;
}

Event::TkrCluster* TkrQueryClustersTool::nearestClusterOutside(
    int view, int layer, double inDistance, const Point& Pcenter) const
{
    return nearestClusterOutsideX(view, layer, inDistance, Pcenter, STANDARDCLUSTERS);
}

Event::TkrCluster* TkrQueryClustersTool::nearestBadClusterOutside(
    int view, int layer, double inDistance, const Point& Pcenter) const
{ 
    return nearestClusterOutsideX(view, layer, inDistance, Pcenter, BADCLUSTERS);
}

Event::TkrCluster* TkrQueryClustersTool::nearestClusterOutsideX(
    int view, int layer, double inDistance, 
    const Point& Pcenter, clusterType type) const
{
    // Purpose and Method: returns the position of the closest cluster
    //    outside of a given distance from a point in the measured direction,
    //    and in the same or adjacent tower in the other direction.
    // Inputs:  view and layer, center and distance
    // Outputs:  Position of nearest cluster
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Event::TkrCluster* nearCluster = 0;
    
    if (!validLayer(layer, type)) return nearCluster;

    const Event::TkrClusterVec clusters = getClusters(view, layer);
    int nhits = clusters.size();
    if (nhits == 0) return nearCluster;
    
    double minDistance = inDistance;
    double maxDistance = 1e6;
    Point Pini(0.,0.,0.);
    Event::TkrClusterVecConItr clusIter = clusters.begin();
    for(; clusIter != clusters.end(); clusIter++)
    {
        const Event::TkrCluster* cluster = (*clusIter);
        if (cluster->hitFlagged()) continue;
        
        Pini = cluster->position();
        
        // Kludge to prevent crashes when z layer incorrect
        //double zDistance   = fabs(Pini.z() - Pcenter.z());
        //if (zDistance > .3) continue;
        
        Vector diff = Pini-Pcenter;
        double measDist;
        double orthDist;
        switch (view) {
            case idents::TkrId::eMeasureX:
                measDist = fabs(diff.x());
                orthDist = fabs(diff.y());
                break;
            case idents::TkrId::eMeasureY:
                measDist = fabs(diff.y());
                orthDist = fabs(diff.x());
                break;
            default:
                measDist = diff.mag();
                orthDist = 0.;
        }
                
        if ( measDist >= minDistance 
            && measDist < maxDistance && orthDist < m_testDistance) 
        {
            maxDistance = measDist;
            nearCluster = const_cast<Event::TkrCluster*>(cluster);
        }
    }
    return nearCluster;
}

int TkrQueryClustersTool::getNumNearHitsInPlane(int layer, int view, double dist,
                                            const Point& x0, const Vector& dir,
                                            bool checkFlag) const
{
    int numHits = 0;

    const Event::TkrClusterVec hitList = getClusters(view, layer);
    int nHitsInPlane = hitList.size();
    
    // move hit to z of requested plane
    double zPlane = m_tkrGeom->getLayerZ(layer, view);
    Vector delta = (zPlane - x0.z())/dir.z()*dir;
    Point xPlane = x0 + delta;

    while(nHitsInPlane--)
    {
        if(checkFlag && (hitList[nHitsInPlane]->hitFlagged())) continue; 

        Vector diff = xPlane - hitList[nHitsInPlane]->position();
        if (view==idents::TkrId::eMeasureX) {
            if (fabs(diff.x())<dist && fabs(diff.y())<m_testDistance) numHits++;
        } else {
            if (fabs(diff.x())<m_testDistance && fabs(diff.y())<dist) numHits++;
        }
    }
    return numHits;
}

int TkrQueryClustersTool::getNumNearHitsInLayer(int layer, double dX, double dY, 
                                       const Point& x0, const Vector& dir,
                                       bool checkFlag) const
{
    if (!validLayer(layer)) return 0;

    int numHits;

    numHits  = getNumNearHitsInPlane(layer, idents::TkrId::eMeasureX, dX, x0, dir, checkFlag);
    numHits += getNumNearHitsInPlane(layer, idents::TkrId::eMeasureY, dY, x0, dir, checkFlag);

    return numHits;
}

int TkrQueryClustersTool::numberOfHitsNear( int layer, double dX, double dY, 
                                       const Point& x0, const Vector dir) const
{
    // Purpose and Method: counts the number of hits in a bilayer 
    //      within a rectangle of sides 2*dX, 2*dY
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
        
    bool checkFlag = false;
    return getNumNearHitsInLayer(layer, dX, dY, x0, dir, checkFlag);
}

int TkrQueryClustersTool::numberOfUUHitsNear( int layer, double dX, double dY, 
                                       const Point& x0, const Vector dir) const
{
    // Purpose and Method: counts the number of un-used hits in a bilayer 
    //      within a rectangle of sides 2*dX, 2*dY
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None

    bool checkFlag = true;
    return getNumNearHitsInLayer(layer, dX, dY, x0, dir, checkFlag);
}

int TkrQueryClustersTool::numberOfHitsNear( int /*view*/, int layer, double inDistance, 
                                           const Point& x0, const Vector dir) const
{
    // Purpose and Method: counts the number of hits within a distance 
    //     "inDistance" in the measurement direction, and within one tower 
    //     in the other direction
    // Inputs:  layer number, inDistance, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    bool checkFlag = false;
    return getNumNearHitsInLayer(layer, inDistance, inDistance, x0, dir, checkFlag);
}

double TkrQueryClustersTool::clusterWidth(Event::TkrCluster* cluster) const
{
    double size = cluster->size();
    int stripsPerLadder = m_tkrGeom->ladderNStrips();
    int nGaps = cluster->lastStrip()/stripsPerLadder
        - cluster->firstStrip()/stripsPerLadder;
    double width = size*m_tkrGeom->siStripPitch() 
        + nGaps*(2*m_tkrGeom->siDeadDistance() + m_tkrGeom->ladderGap());
    return width;
}
