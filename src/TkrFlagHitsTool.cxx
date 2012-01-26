// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrFlagHitsTool.cxx,v 1.3 2008/02/27 22:44:15 lsrea Exp $

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrFlagHitsTool.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/TkrTrkParams.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include <map>

using namespace Event;

class TkrFlagHitsTool : public AlgTool, virtual public ITkrFlagHitsTool 
{
public:

    TkrFlagHitsTool(const std::string& type, const std::string& name, const IInterface* parent);
    //virtual ~TkrFlagHitsTool();

    StatusCode initialize();


    int flagHits(idents::TkrId tkrId, 
        Event::TkrTrackParams inParams, double zIn, 
        double minError, double maxError, double nSigma, 
        Event::TkrTrackParams& outParams,  
        unsigned int& status_bits) const;


private:

/*    
    int flagHits(Point pos, int tower, int layer, int view, 
        double xTower, double yTower,
        double xError, double yError, double hitWidth, double sigmaCut,
        double *paramList, unsigned int& status_bits) const;
*/

    /// pointer to tracker geometry
    ITkrGeometrySvc*  m_tkrGeom;
    /// pointer to badStripsSvc
    //ITkrBadStripsSvc* m_pBadStrips;
    /// Pointer to the failure service
    ITkrFailureModeSvc*  m_failSvc;
    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;
    /// Pointer to SplitsSvc
    ITkrSplitsSvc*       m_splitsSvc;
    /// pointer to GlastSvc
    IGlastDetSvc*        m_detSvc;
    /// Pointer to the G4 propagator
    IPropagator*         m_propagatorTool;
    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

};

namespace {
    int xPosIdx = TkrTrackParams::xPosIdx;
    int yPosIdx = TkrTrackParams::yPosIdx;
    int xSlpIdx = TkrTrackParams::xSlpIdx;
    int ySlpIdx = TkrTrackParams::ySlpIdx;

    double sqrt12_inv = 1./sqrt(12.);
}

// Static factory for instantiation of algtool objects
//static ToolFactory<TkrFlagHitsTool> s_factory;
//const IToolFactory& TkrFlagHitsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrFlagHitsTool);

// Standard Constructor
TkrFlagHitsTool::TkrFlagHitsTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent)
                                 : AlgTool( type, name, parent )
{    
    // Declare additional interface
    declareInterface<ITkrFlagHitsTool>(this);

    return;
}


StatusCode TkrFlagHitsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    //Set the properties
    setProperties();

    if( serviceLocator() ) {
        sc = serviceLocator()->service( "TkrGeometrySvc", m_tkrGeom, true );
        if(sc.isFailure()) {
            throw GaudiException("Could not find TkrGeometrySvc", name(), sc);
            return sc;
        }

    }
    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    if ((sc = serviceLocator()->getService("GlastDetSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
    m_detSvc = dynamic_cast<IGlastDetSvc*>(iService);

    m_failSvc = m_tkrGeom->getTkrFailureModeSvc();
    m_splitsSvc = m_tkrGeom->getTkrSplitsSvc();

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    //Locate a pointer to the G4Propagator
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagatorTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("TkrQueryClustersTool", m_clusTool);
        if (sc.isFailure()) {
            throw GaudiException("Could not find TkrQueryClusterTool", 
                name(), sc);
        }
    }

    log << MSG::INFO << "TkrFlagHitsTool successfully initialized" << endreq;
    return sc;
}

int TkrFlagHitsTool::flagHits(idents::TkrId tkrId, 
                              Event::TkrTrackParams inParams, double zIn, 
                              double minError, double maxError, double sigmaCut, 
                              Event::TkrTrackParams& outParams,  
                              unsigned int& status_bits) const
{
    double towerPitch = m_tkrGeom->towerPitch();
    int xNum = m_tkrGeom->numXTowers();
    int yNum = m_tkrGeom->numYTowers();
    double siThickness = m_tkrGeom->siThickness();
    double siStripPitch = m_tkrGeom->siStripPitch();

    int tower = idents::TowerId(tkrId.getTowerX(), tkrId.getTowerY()).id();
    int layer = m_tkrGeom->trayToBiLayer(tkrId.getTray(), tkrId.getBotTop());
    int view = tkrId.getView();

    double xIn = inParams(xPosIdx);
    double yIn = inParams(yPosIdx);

    int iXTower, iYTower;
    double xTower = 
        m_tkrGeom->truncateCoord(xIn, towerPitch, xNum, iXTower);
    double yTower = 
        m_tkrGeom->truncateCoord(yIn, towerPitch, yNum, iYTower);

    double xError = std::min(maxError, sqrt(inParams(xPosIdx, xPosIdx)));
    double yError = std::min(maxError, sqrt(inParams(yPosIdx, yPosIdx)));
    xError = std::max(xError, minError);
    yError = std::max(yError, minError);

    // this is the projected width of the track in the measurement direction
    // in the silicon. Not clear what we actually want here, but for starters
    // we round up to the next strip width.

    double xHitWidth = inParams(xSlpIdx)*siThickness;
    double yHitWidth = inParams(ySlpIdx)*siThickness;
    double hitWidth;

    if (view==idents::TkrId::eMeasureX) {
        xHitWidth = ceil(xHitWidth/siStripPitch) * siStripPitch;
        hitWidth = xHitWidth;
    } else {
        yHitWidth = ceil(yHitWidth/siStripPitch) * siStripPitch;
        hitWidth = yHitWidth;
    }

    double halfHitWidth = 0.5*hitWidth;

    Point pos(xIn, yIn, zIn);

    int stage = 0;

    // The order of the tests is roughly in the order of increasing localization
    //   of the missing hit... from none for a dead plane, to potentially a
    //   single strip for a dead strips.

    // Do this here to satisfy the conditions for the dreaded goto

    SmartDataPtr<TkrTruncationInfo> truncInfo( 
        m_dataSvc, EventModel::TkrRecon::TkrTruncationInfo );
    Point posBad;

    int tray, face;
    m_tkrGeom->layerToTray(layer, view, tray, face);

    // +++++++++++++++++++++++++++
    // first check for failed plane
    // +++++++++++++++++++++++++++

    if(m_failSvc && !m_failSvc->empty() 
        && m_failSvc->isFailed(tower, layer, view) ) 
    {
        // nothing measured, just flag the hit
        status_bits |= TkrTrackHit::HITISDEADPLN;
        stage = 1;
        return stage;
    }

    // +++++++++++++++++++++++++++
    // Analyze truncated hits
    // +++++++++++++++++++++++++++

    //don't bother if not truncated

    if (truncInfo!=NULL&&truncInfo->isTruncated()) {
        TkrTruncationInfo::TkrTruncationMap*  truncMap = truncInfo->getTruncationMap();
        SortId id(tower, tray, face, view);
        TkrTruncationInfo::TkrTruncationMap::iterator iter = truncMap->find(id);
        if (iter!=truncMap->end() ) {
            TkrTruncatedPlane item = iter->second;
            //SortId sortId = iter->first;
            //std::cout << " FTT: T/T/F "<< sortId.getTower() << " " << sortId.getTray() << " " << sortId.getFace() << std::endl;
            if (item.isTruncated()) {
                // here's where the work begins!!
                // first try: compare extrapolated position to missing strip locations
                // check for RC truncation
                int status = item.getStatus();
                //const intVector&  stripNum = item.getStripNumber();
                int splitPoint  = m_splitsSvc->getSplitPoint(tower, layer, view);
                double stripPitch = m_tkrGeom->siStripPitch();
                int nStrips = m_tkrGeom->ladderNStrips()*m_tkrGeom->nWaferAcross();
                double splitPos = m_detSvc->stripLocalX(splitPoint);
                //double splitPos = (splitPoint - (nStrips-1)*0.5)/stripPitch;             
                bool lowSet  = ((status & TkrTruncatedPlane::END0SET)!=0);
                bool RCHighSet = ((status & TkrTruncatedPlane::RC1SET)!=0);

                double lowPos = splitPos, highPos = splitPos;

                // need a sigma as well as a position
                double pos_err = 
                    ((view==idents::TkrId::eMeasureX) ? xError : yError);
                // Error due to finite SSD thickness - matters at large angles
                // Add them in quadrature
                //    but better would be to add the width to the distance
                double rError=sqrt(pos_err*pos_err);
                double nError = rError*sigmaCut;

                double localX = (view == idents::TkrId::eMeasureX ? xTower : yTower);
                const floatVector stripLocalX = item.getLocalX();
                if(lowSet)    { lowPos  = stripLocalX[0];}
                if(RCHighSet) { highPos = stripLocalX[1];}
                bool truncBit = false;
                if(highPos>lowPos) {
                    if ((localX-lowPos + halfHitWidth)>-nError 
                        && (localX-highPos - halfHitWidth)<nError) truncBit = true;
                }
                // now do the same for the high end (CC1)
                if ((status & TkrTruncatedPlane::CC1SET)!=0) {
                    lowPos  = stripLocalX[2];
                    if ((localX-lowPos + halfHitWidth)>-nError) truncBit = true;
                }
                if (truncBit) {
                    status_bits |= TkrTrackHit::HITISTRUNCATED;

                    stage = 2;
                    return stage;
                }
            }
        }
    }

    // +++++++++++++++++++++++++++
    // Tower edges
    // +++++++++++++++++++++++++++

    // Get the signed distance to the active edge of the tower

    // if these are negative, track misses active area of tower
    // probably no point in constraining hit in this plane

    double xGap, yGap;
    double xActiveWidth, yActiveWidth;

    int nWafer = m_tkrGeom->nWaferAcross();
    double xPitch, yPitch, xSiGap, ySiGap;
    double deadGap = m_tkrGeom->siDeadDistance();

    xPitch = m_tkrGeom->ladderPitch();
    yPitch = m_tkrGeom->waferPitch();
    xSiGap = m_tkrGeom->ladderGap();
    ySiGap = m_tkrGeom->ladderInnerGap();
    if (view==idents::TkrId::eMeasureX) {
        std::swap(xPitch, yPitch);
        std::swap(xSiGap, ySiGap);
    }

    xGap = xSiGap + 2*deadGap;
    xActiveWidth = nWafer*xPitch - xGap;

    yGap = ySiGap + 2*deadGap;
    yActiveWidth = nWafer*yPitch - yGap;


    double xActiveDistTower = 0.5*xActiveWidth - fabs(xTower) + 0.5*xHitWidth;
    double yActiveDistTower = 0.5*yActiveWidth - fabs(yTower) + 0.5*yHitWidth;

    double sigmaXEdgeTower = xActiveDistTower/xError;
    bool nearXEdgeTower = (sigmaXEdgeTower < sigmaCut);

    double sigmaYEdgeTower = yActiveDistTower/yError;
    bool nearYEdgeTower = (sigmaYEdgeTower < sigmaCut);

    // bail here if near outer edge of tower
    if (nearXEdgeTower || nearYEdgeTower) {
        // let's say no measurement for this type of hit...
        // but for future reference:
        // double towerGap = towerPitch - nWafer*xPitch + xGap;

        status_bits |= TkrTrackHit::HITISTWR;

        stage = 3;
        return stage;
    }

    // We're in the active area of the tower, so:

    // +++++++++++++++++++++++++++
    // interwafer gap
    // +++++++++++++++++++++++++++

    // Starting here, we can assign some position information

    double activeWaferSide = m_tkrGeom->siActiveWaferSide();
    // look for internal gaps
    int iXWafer;
    double xWafer = m_tkrGeom->truncateCoord(xTower, xPitch, nWafer, iXWafer);
    double xActiveDistWafer = 0.5*activeWaferSide - fabs(xWafer) + 0.5*xHitWidth;
    bool nearXEdge = (xActiveDistWafer/xError < sigmaCut);

    int iYWafer;
    double yWafer = m_tkrGeom->truncateCoord(yTower, yPitch, nWafer, iYWafer);
    double yActiveDistWafer = 0.5*activeWaferSide - fabs(yWafer) + 0.5*yHitWidth;
    bool nearYEdge = (yActiveDistWafer/yError < sigmaCut);

    //std::cout << nearXEdge << " " << xActiveDistWafer << " " << xError << " " 
    //          << nearYEdge << " " << yActiveDistWafer << " " << yError << " "
    //          << sigmaCut << std::endl;


    if (nearXEdge || nearYEdge) {

        // for checking the handling of edges
        //if(nearXEdge) std::cout << "X: " << xActiveDistWafer << " " << xError << " " 
        //    << sigmaCut << " " << xError*sigmaCut << std::endl;
        //if(nearYEdge) std::cout << "Y: " << yActiveDistWafer << " " << yError << " " 
        //    << sigmaCut << " " << yError*sigmaCut << std::endl;

        if (nearXEdge&&nearYEdge) {
            //call it not measured, could be anything!
            status_bits |= TkrTrackHit::HITISGAP;

            stage = 4;
            return stage;
        }

        // not clear where to put the hit
        // for now:
        // move measured coord to the center of the gap
        // move other coord to the center of the tower

        int sign;
        double xPos, yPos;
        if(nearXEdge) {
            sign = (xWafer>0 ? 1  : -1);
            xPos = pos.x() + sign*(xActiveDistWafer + 0.5*xGap - 0.5*xHitWidth);
            yPos = pos.y() - yTower;
        } else {
            sign = (yWafer>0 ? 1  : -1);
            xPos = pos.x() - xTower;
            yPos = pos.y() + sign*(yActiveDistWafer + 0.5*yGap - 0.5*yHitWidth);
        }

        outParams(xPosIdx) = xPos;
        outParams(xSlpIdx) = 0.;
        outParams(yPosIdx) = yPos;
        outParams(ySlpIdx) = 0;

        if (nearXEdge) { yGap = towerPitch; }
        else           { xGap = towerPitch; }

        double sigmaX = xGap*sqrt12_inv;
        double sigmaY = yGap*sqrt12_inv;

        outParams(xPosIdx, xPosIdx) = sigmaX * sigmaX;
        outParams(yPosIdx, yPosIdx) = sigmaY * sigmaY;

        status_bits |= TkrTrackHit::HITISGAP;
        status_bits |= TkrTrackHit::HASMEASURED;
        if(nearXEdge) {status_bits |= TkrTrackHit::MEASURESX;}
        else          {status_bits |= TkrTrackHit::MEASURESY;}

        stage = 5;
        return stage;
    }

    // +++++++++++++++++++++++++++
    // BadClusters
    // +++++++++++++++++++++++++++

    // should we restrict ourselves to one tower?

    // To do:
    // The phase space for a gap is usually much bigger than for a dead strip
    // so need to compare distance to gap with distance to bad strip, and maybe
    // also look at the number of dead strips with 3 sigma of the hit position
    // before deciding if this is a gap or a dead strip.

    double minDist = 0.0;
    TkrCluster* badCluster = 0;

    for (;;) {
        badCluster = 
            m_clusTool->nearestBadClusterOutside(view, layer, minDist, pos);

        if(!badCluster) break; // no bad clusters found

        double distance, width;

        // here is where we do something about the bad cluster
        Point posBad = badCluster->position();
        Vector diff = pos - posBad;
        distance = fabs(diff[view]);

        // get the cluster width, including gaps
        width = m_clusTool->clusterWidth(badCluster);
        // move minimum distance past this cluster, for next try
        minDist = distance + 0.5*width;
     
        double deltaW = 0.5*(width - hitWidth);
        if(deltaW<0) continue; // cluster doesn't cover the hit

        distance += deltaW;

        double sig_bad = distance/(view==idents::TkrId::eMeasureX ? xError : yError);

        if(sig_bad >= sigmaCut) break; // cluster is too far away, we're done

        outParams(xPosIdx) = badCluster->position().x();
        outParams(xSlpIdx) = 0.;
        outParams(yPosIdx) = badCluster->position().y();
        outParams(ySlpIdx) = 0.;

        double sigma = 2.0*std::max(deltaW, siStripPitch)*sqrt12_inv;
        double sigma_alt = m_tkrGeom->trayWidth()*sqrt12_inv;

        double sigmaX = sigma;
        double sigmaY = sigma_alt;

        if(view==idents::TkrId::eMeasureX) {
            status_bits |= TkrTrackHit::MEASURESX;
        } else {
            std::swap(sigmaX, sigmaY);
            status_bits |= TkrTrackHit::MEASURESY;
        }
        status_bits |= TkrTrackHit::HASMEASURED;
        status_bits |= TkrTrackHit::HITISDEADST;

        outParams(xPosIdx, xPosIdx) = sigmaX * sigmaX;
        outParams(yPosIdx, yPosIdx) = sigmaY * sigmaY;

        stage = 6;
        return stage;
    } 

    // +++++++++++++++++++++++++++
    // nothing left to try, flag hit as unknown. 
    // Caller will handle this.
    // +++++++++++++++++++++++++++

    status_bits |= TkrTrackHit::HITISUNKNOWN;
    stage = -1;
    //std::cout << pos << " t/l/v " << tower << " " << layer << " " << view 
    //    << " width " << hitWidth << " xyErr " << xError << " " << yError << std::endl;

    return stage;
}
