/** @file TkrHitTruncationTool.cxx

* @class TkrHitTruncationTool
*
* @brief Generates hit truncation information for an event
*
* @author Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrHitTruncationTool.cxx,v 1.6 2011/01/19 00:44:22 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IDataProviderSvc.h"

#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"

#include "Event/TopLevel/EventModel.h"

#include "TkrHitTruncationTool.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

// constants defined at file scope
namespace {
    int _nStrips;
    int _nLadderStrips;
    double _activeGap;
    double _stripPitch;


    bool debug;

    int numRCTrunc;
    int numCCTrunc;
    int count;
}

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrHitTruncationTool::TkrHitTruncationTool(const std::string& type, 
                                           const std::string& name, 
                                           const IInterface* parent)
                                           : AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrHitTruncationTool>(this);

    return;
}

StatusCode TkrHitTruncationTool::initialize()
{
    // Purpose and Method: 
    // Outputs:
    // Dependencies: 
    // Restrictions and Caveats:  

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
        //Locate and store a pointer to the detector service
    if ( (sc = service("GlastDetSvc", m_detSvc, true)).isFailure() ) {
        throw GaudiException("GlastDetSvc not found", name(), sc);
    } 

    //Locate and store a pointer to the geometry service
    if ( (sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure() ) {
        throw GaudiException("TkrGeometrySvc not found", name(), sc);
    } 
    m_splitsSvc = m_tkrGeom->getTkrSplitsSvc();

    m_dataSvc = 0;
    sc = serviceLocator()->service( "EventDataSvc", m_dataSvc, true );
    if(sc.isFailure()){
        log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
        return sc;
    }

    //Locate and store a pointer to the incident service
    IIncidentSvc* incSvc = 0;
    if ((sc = serviceLocator()->service("IncidentSvc", incSvc, true)).isFailure())
    {
        throw GaudiException("Service [IncidentSvc] not found", name(), sc);
    }

    //set up listener for IncidentSvc
    incSvc->addListener(this, "BeginEvent", 100);

    _nLadderStrips = m_tkrGeom->ladderNStrips();
    _nStrips       = _nLadderStrips*m_tkrGeom->nWaferAcross();
    _activeGap     = 0.5*m_tkrGeom->ladderGap() + m_tkrGeom->siDeadDistance();
    _stripPitch    = m_tkrGeom->siStripPitch();

    log << MSG::DEBUG;
    debug = log.isActive();
    log << endreq;

    numRCTrunc = numCCTrunc = count = 0;

    return sc;
}

StatusCode TkrHitTruncationTool::analyzeDigis()
{
    // Purpose and Method: finds the truncated regions in the current event 
    //    by counting hits 
    // Inputs:  
    // Outputs: 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    using namespace idents;
    using namespace Event;

    MsgStream log(msgSvc(), name());

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    count++;

    // First, the collection of TkrDigis is retrieved from the TDS
    SmartDataPtr<TkrDigiCol> digiCol( m_dataSvc, EventModel::Digi::TkrDigiCol );

    if (digiCol == 0) { return sc; }

    TkrTruncationInfo* truncationInfo = SmartDataPtr<TkrTruncationInfo>(
        m_dataSvc,EventModel::TkrRecon::TkrTruncationInfo);
    if(truncationInfo==0) truncationInfo = new TkrTruncationInfo();
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTruncationInfo,
        truncationInfo);
    // Create the TkrTruncationInfo TDS object
    TkrTruncationInfo::TkrTruncationMap* truncMap = truncationInfo->getTruncationMap();

    TkrDigiCol::const_iterator pTkrDigi = digiCol->begin();
    // do the strip counts and generate the RC truncation information
    for (; pTkrDigi!= digiCol->end(); pTkrDigi++) {
        TkrDigi* pDigi = *pTkrDigi;

        TowerId towerId = pDigi->getTower();
        //int towerX = towerId.ix();
        //int towerY = towerId.iy();
        int tower  = towerId.id();

        int layer = pDigi->getBilayer();
        int view  = pDigi->getView();
        int tray, face;
        m_tkrGeom->layerToTray(layer, view, tray, face);
        int plane = m_tkrGeom->trayToPlane(tray, face);

        intVector stripCount(2,0);

        int maxStrips[2];
        maxStrips[0]   = m_splitsSvc->getMaxStrips(tower, layer, view, 0);
        maxStrips[1]   = m_splitsSvc->getMaxStrips(tower, layer, view, 1);
        int splitStrip = m_splitsSvc->getSplitPoint(tower, layer, view);

        int lastC0Strip  = pDigi->getLastController0Strip();

        intVector stripNumber(4,0);
        stripNumber[0] = -1;          // highest low-side strip  (for RC0 and CC0)
        stripNumber[1] = _nStrips;     // lowest high-side strip  (for RC1)
        stripNumber[2] = _nStrips;     // highest high-side strip (for CC1)
        stripNumber[3] = splitStrip;  // split-point strip

        //check for read controller truncation
        int nHits = pDigi->getNumHits();
        int status = 0;
        int end;
        int hit;
        bool firstHigh = true;
        for (hit=0; hit<nHits; ++hit) {
            int strip = pDigi->getHit(hit);
            end = (strip<=lastC0Strip ? 0 : 1);
            if(end==0) stripNumber[0] = strip;
            if(end==1) {
                stripNumber[2] = strip;
                if (firstHigh) {
                    stripNumber[1] = strip;
                    firstHigh = false;
                }
            }
            stripCount[end]++;
        }
        for(end=0; end<2; ++end) {
            // fill some stuff for the cable test
            if (stripCount[end]<maxStrips[end]) {
                //skip it
            } else if (stripCount[end]>maxStrips[end]) {
                TkrTruncatedPlane::addStatusBit(status, end, TkrTruncatedPlane::RCOVER);
            } else  {
                TkrTruncatedPlane::addStatusBit(status, end, TkrTruncatedPlane::RC);
            }
        }
   
        floatVector localX(4,0);

        // get the limits for the dead regions
        // tricky for -1 and nStrips (no limits) because neither is a legal strip number
        int strip = std::max(stripNumber[0], 0);
        localX[0] = m_detSvc->stripLocalX(strip) + 0.5*_stripPitch -
            (stripNumber[0]==-1 ? _stripPitch : 0);

        strip = std::min(stripNumber[1], _nStrips-1);
        localX[1] = m_detSvc->stripLocalX(strip) - 0.5*_stripPitch +
            (stripNumber[1]==_nStrips ? _stripPitch : 0);

        strip = std::min(stripNumber[2], _nStrips-1);
        localX[2] = m_detSvc->stripLocalX(strip) + 0.5*_stripPitch +
            (stripNumber[2]==_nStrips ? _stripPitch : 0);

        strip = std::max(stripNumber[3], -1);
        localX[3] = m_detSvc->stripLocalX(strip) + 0.5*_stripPitch -
            (stripNumber[3]==-1 ? _stripPitch : 0);

        //std::cout << localX[0] << " " << localX[1] << " "  << localX[2] << " "  
        //  << localX[3] << std::endl;
        float planeZ = m_tkrGeom->getPlaneZ(plane);
        // make and store the TkrTruncatedPlane
        TkrTruncatedPlane item(status, stripCount, stripNumber, localX, planeZ);
        truncationInfo->addItem(tower, tray, face, view, item);
    }

    // now go through the items and do the cable truncation calculation
    int cableHits[8] = {0,0,0,0,0,0,0,0};
    int cableBufferSize = m_splitsSvc->getCableBufferSize();
    int tower0 = -1;

    TkrTruncationInfo::TkrTruncationMap::iterator iter = truncMap->begin();
    for(; iter!=truncMap->end(); ++iter) {
        SortId id = iter->first;
        int tower = id.getTower();
        if (tower!=tower0) {
            int cable;
            for(cable=0; cable<8; ++cable) { cableHits[cable] = 0; }
            tower0 = tower;
        }
        int tray  = id.getTray();
        int face  = id.getFace();
        int layer, view;
        m_tkrGeom->trayToLayer(tray, face, layer, view);
        TkrTruncatedPlane& trunc = iter->second;
        int end;
        const intVector& numStrips = trunc.getStripCount();
        for (end=0; end<2; ++end) {
            int index = m_splitsSvc->getCableIndex(layer, view, end);
            cableHits[index] += numStrips[end];
            int numHits = cableHits[index];
            if (numHits<cableBufferSize) {
                // skip it
            } else if (numHits > cableBufferSize) {
                trunc.setStatusBit(end, TkrTruncatedPlane::CCOVER);
            } else {
                trunc.setStatusBit(end, TkrTruncatedPlane::CC);
            }
        }
    }
    // finish up:
    // set the truncation counts and deal with planes with no truncations
    //int sizeBefore = truncMap->size();
    int nRCTrunc = 0;
    int nCCTrunc = 0;
    
    std::vector<TkrTruncationInfo::TkrTruncationMap::iterator> iterVec;
    
    iter = truncMap->begin();
    
    for(;iter!=truncMap->end(); ++iter ) {
        const TkrTruncatedPlane trunc = iter->second;
        const int status = trunc.getStatus();
        if (status==0) {
            iterVec.push_back(iter);
            continue;
        }
        if( (status & TkrTruncatedPlane::RCSET)!=0 ) nRCTrunc++;
        if( (status & TkrTruncatedPlane::CCSET)!=0 ) nCCTrunc++;
    }
   
    // remove the TruncatedPlanes with status==0
    unsigned int i;
    for(i=0; i<iterVec.size(); ++i) {
        truncMap->erase(iterVec[i]);
    }
    
    truncationInfo->setCCTrunc(nCCTrunc);
    truncationInfo->setRCTrunc(nRCTrunc);
    // mark it done!
    truncationInfo->setDone();
    
    if (debug) {
        log << MSG::DEBUG;
        iter = truncMap->begin();
        tower0 = -1;
        for(; iter!=truncMap->end(); ++iter) {
            SortId id = iter->first;
            int tower = id.getTower();
            //if (tower!=tower0) std::cout << "Tower " << tower << std::endl;
            tower0 = tower;
            int tray  = id.getTray();
            int face  = id.getFace();
            int layer, view;
            m_tkrGeom->trayToLayer(tray, face, layer, view);
            TkrTruncatedPlane trunc = iter->second;
            const intVector& numStrips = trunc.getStripCount();
            const intVector& stripNumber = trunc.getStripNumber();
            const int status   = trunc.getStatus();
            const floatVector& localX = trunc.getLocalX();

            log << "Twr/Tray/face/layer/view " <<tower << "/" << tray << "/" << face
                << " " << layer << " " << view
                << ", status " << status 
                << ", #Strips " << numStrips[0] << "/" << numStrips[1]
                << ", # " << stripNumber[0]  << "/" << stripNumber[1] 
                << "/" << stripNumber[2]
                << ", X " << localX[0]  << "/" << localX[1]
                << "/" << localX[2]
                << endreq;
        }

        log << "Truncation count: " << truncationInfo->getNumRCTruncated() 
            << " " << truncationInfo->getNumCCTruncated() << endreq;
    }

    numRCTrunc += truncationInfo->getNumRCTruncated();
    numCCTrunc += truncationInfo->getNumCCTruncated();
    
    //std::cout << MSG::DEBUG << "event " << count << " truncs " 
    //    << truncationInfo->getNumRCTruncated() 
    //    << " " << truncationInfo->getNumCCTruncated() << std::endl;

    return sc;
}  

StatusCode TkrHitTruncationTool::finalize()
{
    std::cout << "number of truncation records " << numRCTrunc << " " 
        <<numCCTrunc << std::endl;
    std::cout << "number of calls " << count << std::endl;
    return StatusCode::SUCCESS;
}

void TkrHitTruncationTool::handle(const Incident & inc) 
{    
    MsgStream log(msgSvc(), name());

    if(inc.type()=="BeginEvent") {
        //std::cout << "handle called at BeginEvent" << std::endl;
        m_newEvent = true;
    }
}

double TkrHitTruncationTool::getDistanceToTruncation(
    int tower, int tray, int face, int view, double localX)
{

    // gets the distance to the edge of the nearest truncated region 
    //    in the requested plane.
    // the search is confined to the current tower

    using namespace Event;

    if(m_newEvent) {
        TkrTruncationInfo* truncationInfo = 
            SmartDataPtr<TkrTruncationInfo>(m_dataSvc, 
            EventModel::TkrRecon::TkrTruncationInfo);
        m_truncMap = 0;
        if (truncationInfo->isTruncated()) {
             m_truncMap = truncationInfo->getTruncationMap();
        }
        m_newEvent = false;
    }
 
    // no truncation returns a large distance
    // which means that it is far away from any truncated region

    double distance = 9999.0;

    if(m_truncMap) {
        SortId id(tower, tray, face, view);
        TkrTruncationInfo::TkrTruncationMap::iterator iter = m_truncMap->find(id);
        if (iter!=m_truncMap->end() ) {
            TkrTruncatedPlane item = iter->second;
            //SortId sortId = iter->first;
            //std::cout << " FTT: T/T/F "<< sortId.getTower() << " " 
            //   << sortId.getTray() << " " << sortId.getFace() << std::endl;
            if (item.isTruncated()) {
                // here's where the work begins!!
                // first try: compare extrapolated position to missing strip locations
                // check for RC truncation
                int status = item.getStatus();
                int layer = m_tkrGeom->trayToBiLayer(tray, face);
                int splitPoint  = m_splitsSvc->getSplitPoint(tower, layer, view);
                double splitPos = m_detSvc->stripLocalX(splitPoint);

                bool lowSet  = ((status & TkrTruncatedPlane::END0SET)!=0);
                bool RCHighSet = ((status & TkrTruncatedPlane::RC1SET)!=0);

                double lowPos = splitPos;
                double highPos = splitPos;

                const floatVector stripLocalX = item.getLocalX();
                if(lowSet)    { lowPos  = stripLocalX[0];}
                if(RCHighSet) { highPos = stripLocalX[1];}
                
                // in the other tests, negative means in the insensitive region
                // so here, negative should mean in the truncated range
                // that is, between the last possible low hit 
                //   and the first possible high hit
                // If either of the two distances is near zero or above, it's
                //   in the *good* area, and should have fired, so we want the 
                //   largest of the distances.
                // Then if any of the tests yields a ~ negative number it means
                //   that the hit was in one of the dead regions.
                //
                // Sorry, this is very confusing...

                double dist1 = lowPos - localX ;
                double dist2 = localX - highPos;
                distance = std::max(dist1, dist2);
               
                // now do the same for the high end (CC1)
                if ((status & TkrTruncatedPlane::CC1SET)!=0) {
                    lowPos  = stripLocalX[2];
                    //highPos = 0.5*_nStrips*_stripPitch;

                    // again, if this is negative, it's in the truncated region
                    // There is no "highPos" because all hits are below the top
                    //   by definition.

                    dist1 = lowPos - localX;
                    //dist2 = localX - highPos;

                    // The hit can't be in both dead regions at the same time
                    // So if either is in a dead region, that's the one we want.
                    // which is to say, the most negative of the two.
                    distance = std::min(distance, dist1);

                    //double distance1;
                    //if (dist1>=0&&dist2>=0) { distance1 = std::min(dist1, dist2); }
                    //else if (dist1<0)     { distance1 = dist1; }
                    //else                  { distance1 = dist2; }

                    //if(fabs(distance1)<fabs(distance)) distance = distance1;
                }
            }
        }
    }

    return distance;
}

double TkrHitTruncationTool::getDistanceToTruncation(idents::TkrId id, Vector towerPos)
{
    int tray = id.getTray();
    int face = id.getBotTop();
    int view = id.getView();
    int iX = id.getTowerX();
    int iY = id.getTowerY();
    int tower = idents::TowerId(iX, iY).id();
    double localX = (view==0 ? towerPos.x() : towerPos.y());
    double result = getDistanceToTruncation(tower, tray, face, view, localX);
    return result;
}

double TkrHitTruncationTool::getDistanceToTruncation(
    int tower, int plane, Vector towerPos)
{
    int tray = m_tkrGeom->planeToTray(plane);
    int face = m_tkrGeom->planeToBotTop(plane);
    int view = m_tkrGeom->getView(plane);
    double localX = (view==0 ? towerPos.x() : towerPos.y());
    double result = getDistanceToTruncation(tower, tray, face, view, localX);
    return result;
}
