/**
* @class TkrReasonsTool
*
* @brief Tool for evaluating the status of missing clusters on tracks
*
* @author The Tracking Software Group
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrReasonsTool.cxx,v 1.1 2011/03/26 22:32:12 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "TkrUtil/ITkrReasonsTool.h"

// TkrRecon utilities
//#include "src/Track/TkrControl.h"
// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/ITkrHitTruncationTool.h"
#include "TkrUtil/TkrTrkParams.h"
#include "TkrUtil/TkrCovMatrix.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

using namespace Event;


namespace {
    int xPosIdx = TkrTrackParams::xPosIdx;
    int yPosIdx = TkrTrackParams::yPosIdx;
    int xSlpIdx = TkrTrackParams::xSlpIdx;
    int ySlpIdx = TkrTrackParams::ySlpIdx;

    double _siThickness;
    double _towerPitch;
    int    _numPlanes;
    int    _numXTowers;
    int    _numYTowers;
    double _trayWidth;
    double _sigma_alt;
    double _siResolution;
    double _ladderPitch, _waferPitch;
    double _ladderGap, _ladderInnerGap;
    double _deadGap;
    int    _nWafer;
    double _activeWaferSide;

    double _siStripPitch;
    int    _ladderNStrips;

    double _xTower;
    double _yTower;
    double _tower;
    int    _layer;
    int    _view;
    int    _tray;
    int    _face;
    idents::TkrId _tkrId;
    double _planeZ;
    int    _plane;

    Point  _end_pos;
}

class TkrReasonsTool : public AlgTool, virtual public ITkrReasonsTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrReasonsTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrReasonsTool() {}

    /// @brief 

    StatusCode initialize();

    virtual void setParams(const Point& end_pos, const int next_plane) const;
    virtual bool isFailed() const;
    virtual Vector getEdgeDistance() const;
    virtual Vector getGapDistance() const;
    virtual double getTruncDistance() const;
    virtual double getBadClusterDistance() const;
    virtual double getMinimumDistance(const Point& end_pos, const int next_plane) const;

private:

    IDataProviderSvc*    m_dataSvc;

        /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*     m_tkrGeom;

    /// Pointer to the failure service
    ITkrFailureModeSvc*  m_failSvc;

    /// Pointer to SplitsSvc
    ITkrSplitsSvc*       m_splitsSvc;

    /// pointer to GlastSvc
    IGlastDetSvc*        m_detSvc;

    /// pointer to truncation tool
    ITkrHitTruncationTool* m_truncTool;

        /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// Pointer to the TkrClusters
    TkrClusterCol*       m_clusters; 

};

static ToolFactory<TkrReasonsTool> s_factory;
const IToolFactory& TkrReasonsToolFactory = s_factory;

TkrReasonsTool::TkrReasonsTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrReasonsTool>(this);

    return;
}

StatusCode TkrReasonsTool::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    //Locate and store a pointer to the data service
    IService* iService = 0;

    if ((sc = serviceLocator()->getService("GlastDetSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
    m_detSvc = dynamic_cast<IGlastDetSvc*>(iService);

    //Locate and store a pointer to the geometry service
    iService = 0;
    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeom   = dynamic_cast<ITkrGeometrySvc*>(iService);
    m_splitsSvc = m_tkrGeom->getTkrSplitsSvc();
    m_failSvc   = m_tkrGeom->getTkrFailureModeSvc();
    if(!m_splitsSvc) {
        throw GaudiException("Service [TkrSplitsSvc] not found", name(), sc);
    }
    if(!m_failSvc) {
        throw GaudiException("Service [TkrFailureModeSvc] not found", name(), sc);
    }

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // Locate a pointer to the TrkHitTruncationTool
    if ((toolSvc()->retrieveTool("TkrHitTruncationTool", m_truncTool)).isFailure())
    {
        throw GaudiException("ToolSvc could not find TkrHitTruncationTool", name(), sc);
    }

    //Locate a pointer to the TrkQueryClusterTool
    if ((toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("ToolSvc could not find TkrQueryClusterTool", name(), sc);
    }

    _siThickness  = m_tkrGeom->siThickness();
    _towerPitch   = m_tkrGeom->towerPitch();
    _numPlanes    = m_tkrGeom->numPlanes();
    _numXTowers   = m_tkrGeom->numXTowers();
    _numYTowers   = m_tkrGeom->numYTowers();
    _trayWidth    = m_tkrGeom->trayWidth();
    _sigma_alt    = _trayWidth/sqrt(12.);
    _siResolution = m_tkrGeom->siResolution();

    _ladderPitch     = m_tkrGeom->ladderPitch();
    _waferPitch      = m_tkrGeom->waferPitch();
    _ladderGap       = m_tkrGeom->ladderGap();
    _ladderInnerGap  = m_tkrGeom->ladderInnerGap();
    _deadGap         = m_tkrGeom->siDeadDistance();
    _nWafer          = m_tkrGeom->nWaferAcross();
    _activeWaferSide = m_tkrGeom->siActiveWaferSide();
    _siStripPitch    = m_tkrGeom->siStripPitch();
    _ladderNStrips   = m_tkrGeom->ladderNStrips();

    return sc;
}

void TkrReasonsTool::setParams(const Point& end_pos, const int next_plane0) const 
{
    //test for next_plane
    int next_plane;
    if(next_plane0==-1) {
        next_plane = m_tkrGeom->getPlane(end_pos.z());
    } else {
        next_plane = next_plane0;
    }
    //if(planeGuess!=next_plane) std::cout << "nextplane/guess " 
    //    << next_plane << " " << planeGuess << std::endl;
    int iXTower, iYTower;
    _xTower = 
        m_tkrGeom->truncateCoord(end_pos.x(), _towerPitch, _numXTowers, iXTower);
     _yTower = 
        m_tkrGeom->truncateCoord(end_pos.y(), _towerPitch, _numYTowers, iYTower);

    _view  = m_tkrGeom->getView(next_plane);
    _layer = m_tkrGeom->getLayer(next_plane);
    m_tkrGeom->layerToTray(_layer, _view, _tray, _face);

    _tower = idents::TowerId(iXTower, iYTower).id();
    _tkrId = idents::TkrId(iXTower, iYTower, _tray, 
        (_face == idents::TkrId::eTKRSiTop), _view);
    _planeZ = m_tkrGeom->getPlaneZ(next_plane);
    _plane  = next_plane;
    _end_pos = end_pos;
}

bool TkrReasonsTool::isFailed() const
{
    // first check for failed plane
    bool failed = false;
    if(m_failSvc && !m_failSvc->empty() 
        && m_failSvc->isFailed(_tower, _layer, _view) ) failed = true;
    return failed;
}

Vector TkrReasonsTool::getEdgeDistance() const 
{
    //double xPitch, yPitch, xSiGap, ySiGap;

    double xActiveDistTower = 0;
    double yActiveDistTower = 0;
    double xPitch = _ladderPitch;
    double yPitch = _waferPitch;
    double xSiGap = _ladderGap;
    double ySiGap = _ladderInnerGap;
    if (_view==idents::TkrId::eMeasureX) {
        std::swap(xPitch, yPitch);
        std::swap(xSiGap, ySiGap);
    }

    // if these are negative, track misses active area of tower
    // probably no point in constraining hit in this plane
    double xGap = xSiGap + 2*_deadGap;
    xActiveDistTower = 0.5*(_nWafer*xPitch - xGap) - fabs(_xTower);
    //double xError = sqrt(next_params(xPosIdx,xPosIdx));

    double yGap = ySiGap + 2*_deadGap;
    yActiveDistTower = 0.5*(_nWafer*yPitch - yGap) - fabs(_yTower);
    //double yError = sqrt(next_params(yPosIdx,yPosIdx));

    Vector vec = Vector(xActiveDistTower, yActiveDistTower, 0);

    return vec;
}

Vector TkrReasonsTool::getGapDistance() const 
{
    double xPitch = _ladderPitch;
    double yPitch = _waferPitch;
    double xSiGap = _ladderGap;
    double ySiGap = _ladderInnerGap;
    if (_view==idents::TkrId::eMeasureX) {
        std::swap(xPitch, yPitch);
        std::swap(xSiGap, ySiGap);
    }


    int iXWafer;
    double xWafer = m_tkrGeom->truncateCoord(_xTower, xPitch, _nWafer, iXWafer);
    double xActiveDistWafer = 0.5*_activeWaferSide - fabs(xWafer);

    int iYWafer;
    double yWafer = m_tkrGeom->truncateCoord(_yTower, yPitch, _nWafer, iYWafer);
    double yActiveDistWafer = 0.5*_activeWaferSide - fabs(yWafer);

    Vector vec = Vector(xActiveDistWafer, yActiveDistWafer, 0);

    return vec;

}

double TkrReasonsTool::getTruncDistance() const 
{
    Vector towerPos = Vector(_xTower, _yTower, _planeZ);
    double dist = m_truncTool->getDistanceToTruncation(_tower, _plane, towerPos );

    return dist;
}

double TkrReasonsTool::getBadClusterDistance() const
{
    TkrCluster* badCluster = 
        m_clusTool->nearestBadClusterOutside(_view, _layer, 0.0, _end_pos);

    double distance = 9999.0;
    double width;

    if(badCluster) {
        // here is where we do something about the bad cluster
        width = m_clusTool->clusterWidth(badCluster);
        Point pos = badCluster->position();
        Vector diff = _end_pos - pos;
        distance = fabs(diff[_view]);
        // get the cluster width, including gaps
    }
    return distance;
}

double TkrReasonsTool::getMinimumDistance(const Point& end_pos, const int next_plane) const
{
    double distance = 9999.0;
    double minDist;
    Vector dist;

    setParams(end_pos, next_plane);

    if (isFailed()) {
        distance = 0.0;
        return distance;
    }
    dist = getEdgeDistance();
    minDist = std::min(dist.x(), dist.y());
    distance = std::min(distance, minDist);

    dist = getGapDistance();
    minDist = std::min(dist.x(), dist.y());
    distance = std::min(distance, minDist);

    distance = std::min(distance, getTruncDistance());
    distance = std::min(distance, getBadClusterDistance());

    return distance;
}
