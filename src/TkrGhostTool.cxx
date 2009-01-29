/**
* @class TkrGhostTool
*
* @brief Implements a Gaudi Tool for setting the candidate track energies before 
*        the track fit
*
* @author The Tracking Software Group
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrGhostTool.cxx,v 1.3 2009/01/22 01:39:14 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/TkrDigi.h"

#include "LdfEvent/Gem.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGhostTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include <iomanip>
#include <map>

class TkrGhostTool : public AlgTool, virtual public ITkrGhostTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrGhostTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrGhostTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    StatusCode getTkrVector(unsigned short& tkrVector);
    StatusCode calculateTkrVector(
        Event::TkrClusterCol* pCol, unsigned short& towerBits);
    StatusCode calculateTkrVector(
        Event::TkrDigiCol* pCol, unsigned short& towerBits);
    StatusCode flagSingles();
    StatusCode flagEarlyHits(Event::TkrClusterCol* col=0);
    StatusCode flagEarlyTracks();
    StatusCode flagEarlyVertices();

    /// @brief Tool for identifying and flagging ghost clusters

private:

    Event::DigiEvent*      m_digiEvent;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
    IGlastDetSvc*          m_pDetSvc;
    ITkrGeometrySvc*       m_tkrGeom;

    int m_numTowers;
    int m_numLayers;
    std::vector<unsigned int> towerBits;
};

static ToolFactory<TkrGhostTool> s_factory;
const IToolFactory& TkrGhostToolFactory = s_factory;

TkrGhostTool::TkrGhostTool(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrGhostTool>(this);

    return;
}

StatusCode TkrGhostTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    sc = service("GlastDetSvc", m_pDetSvc);
    if(sc.isFailure()) {
        log << MSG::ERROR << "GlastDetSvc not found!" << endreq;
        return sc;
    }
    sc = service("TkrGeometrySvc", m_tkrGeom);
    if(sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc not found!" << endreq;
        return sc;
    }

    int numX, numY;
    m_pDetSvc->getNumericConstByName("xNum", &numX);
    m_pDetSvc->getNumericConstByName("yNum", &numY);  
    m_numTowers = numX*numY;

    m_numLayers = m_tkrGeom->numLayers();

    towerBits.resize(m_numTowers,0);

    return sc;
}

StatusCode TkrGhostTool::getTkrVector(unsigned short& tkrVector)
{
    MsgStream log(msgSvc(), name());

    tkrVector = 0;

    //get the tkrVector
    // Retrieve the Event Summary data for this event
    SmartDataPtr<LdfEvent::Gem> gem(m_dataSvc, "/Event/Gem");

    if (!gem) {
        log << MSG::DEBUG << "No GEM found on TDS" << endreq;
        return StatusCode::FAILURE;
    }

    tkrVector= gem->tkrVector();
    return StatusCode::SUCCESS;
}

StatusCode TkrGhostTool::calculateTkrVector(
    Event::TkrClusterCol* clusterCol,
    unsigned short& trigBits)
{
    //set up the tower bits vector
    towerVec clusterTrigger(m_numTowers);
    int i;
    for(i=0;i<m_numTowers;++i) {
        clusterTrigger[i] = new TkrTowerBits();
    }

    //fill the bits vector
    int clusSize = clusterCol->size();
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        clus->clearStatusBits(Event::TkrCluster::maskZAPGHOSTS);
        int tower = clus->tower();
        clusterTrigger[tower]->setBit(clus);
    }

    //generate the tower trigger word
    trigBits = 0;
    for(i=0;i<m_numTowers;++i) {
        towerBits[i] = clusterTrigger[i]->getTriggeredBits();
        if(towerBits[i]>0) {trigBits |= (1<<i);}
    }
    return StatusCode::SUCCESS;
}

StatusCode TkrGhostTool::calculateTkrVector(
    Event::TkrDigiCol* digiCol,
    unsigned short& trigBits)
{
    //set up the tower bits vector
    towerVec digiTrigger(m_numTowers);

    int i;
    for(i=0;i<m_numTowers;++i) {
        digiTrigger[i] = new TkrTowerBits();
    }

    trigBits = 0;

    //fill the bits vector
    int digiSize = digiCol->size();
    for (i=0;i<digiSize;++i) {
        Event::TkrDigi* digi = (*digiCol)[i];
        if (digi->getNumHits()==0) continue;
        int tower = (digi->getTower()).id();
        digiTrigger[tower]->setBit(digi);
    }

    //generate the tower trigger word
    trigBits = 0;
    for(i=0;i<m_numTowers;++i) {
        towerBits[i] = digiTrigger[i]->getTriggeredBits();
        if(towerBits[i]>0) {trigBits |= (1<<i);}
    }
    return StatusCode::SUCCESS;
}


StatusCode TkrGhostTool::flagSingles() 
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    std::map<int, int> hC[3];

    //get the clusters
    SmartDataPtr<Event::TkrClusterCol> 
        clusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);

    int clusSize = clusterCol->size();
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());

        int tower = twrid.id();
        int layer = clus->getLayer();
        int view  = tkrid.getView();
        int end   = clus->getEnd();
        int index = 1000*tower + 2*layer + view ;

        if(hC[end].find(index)==hC[end].end()) hC[end][index] = 0;
        hC[end][index]++;
    }

    // Now mark lone hits

    for(i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        int tower = clus->tower();
        int layer = clus->getLayer();
        int end   = clus->getEnd();
        idents::TkrId id = clus->getTkrId();
        int view = id.getView();
        int index = 1000*tower + 2*layer + view;

        int totCount;
        int endCount[2] = {0,0};
        if(hC[0].find(index)!=hC[0].end()) {
            endCount[0] = hC[0][index];
        }
        if(hC[1].find(index)!=hC[1].end()) {
            endCount[1] = hC[1][index];
        }
        totCount = endCount[0] + endCount[1];

        if(end==2) { // this might never happen... 
            if(totCount==0) {
                clus->setStatusBits(Event::TkrCluster::maskALONE);
                clus->setStatusBits(Event::TkrCluster::maskALONEEND);
            }
        } else {
            if(totCount==1) clus->setStatusBits(Event::TkrCluster::maskALONE);
            if(endCount[end]==1) clus->setStatusBits(Event::TkrCluster::maskALONEEND);
        }
    }
    return sc;
}

StatusCode TkrGhostTool::flagEarlyHits(Event::TkrClusterCol* clusterCol)
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    unsigned short tkrVector;
    sc = getTkrVector(tkrVector);
    if(sc.isFailure()) return sc;

    //get the clusters
    if(clusterCol==0) {
        SmartDataPtr<Event::TkrClusterCol> 
            pClusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);
        clusterCol = pClusterCol;
    }

    unsigned short trigBits;
    sc = calculateTkrVector(clusterCol, trigBits);
    if(sc.isFailure()) return sc;

    if((tkrVector&trigBits)!=tkrVector) {
        log << MSG::DEBUG;
        if (log.isActive()) {
            log <<"tkrVector and calculated tower trigger disagree:" << endreq 
                << "tkrVector = " << std::hex << tkrVector 
            << ", calculation = " << trigBits << std::dec << endreq;
        }
    }

    // Flag the hits that would have made the trigger
    // This is an assumption... for a three-in-a-row, any one
    //   missing hit would disable the trigger
    // Also, flag ToT==255

    int clusSize = clusterCol->size();
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        // flat the ToT225s
        if(clus->getRawToT()==255) clus->setStatusBits(Event::TkrCluster::mask255);
        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
        int tower = twrid.id();
        // we only look for ghosts in towers with software trigger but no hardware trigger
        if((tkrVector&(1<<tower))!=0) continue;
        if((trigBits&(1<<tower))==0) continue;

        int layer = clus->getLayer();
        if(towerBits[tower]&(1<<layer)) {
            clus->setStatusBits(Event::TkrCluster::maskGHOST);

            log << MSG::DEBUG << "Ghost bit set for cluster " 
                << i << ", t/l "  << tower << ", " << layer  << ", isGhost "  
                << clus->isSet(Event::TkrCluster::maskGHOST) <<endreq;
            log << "tower bits " << std::hex << towerBits[tower] 
            << " "  << (1<<layer) << std::dec << endreq;
        }
    }
    return sc;
}

StatusCode TkrGhostTool::flagEarlyTracks()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Flag the hits on tracks containing ghosts or 255's
    //get the tracks
    SmartDataPtr<Event::TkrTrackCol> 
        trackCol(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);

    int trackCount = 0;
    Event::TkrTrackColConPtr tcolIter = trackCol->begin();
    for(tcolIter; tcolIter!=trackCol->end();++tcolIter,++trackCount) {
        Event::TkrTrack* track = *tcolIter;
        track->clearStatusBits(Event::TkrTrack::GHOST);
        Event::TkrTrackHitVecItr pHit = track->begin();

        int ghostCount = 0;
        int _255Count = 0;
        while(pHit != track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;
            bool is255      = pClus->isSet(Event::TkrCluster::mask255);
            bool isGhost    = pClus->isSet(Event::TkrCluster::maskGHOST);
            bool isAlone    = pClus->isSet(Event::TkrCluster::maskALONE);
            bool isAloneEnd = pClus->isSet(Event::TkrCluster::maskALONEEND);
            if(is255&&isAloneEnd) _255Count++;
            if(isGhost)           ghostCount++;
        }

        if(_255Count==0&&ghostCount==0) continue;

        // this is a ghost track!
        track->setStatusBit(Event::TkrTrack::GHOST);

        // flag the remaining clusters on this track
        pHit = track->begin();
        while(pHit!=track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;
            if(_255Count>0||ghostCount>0) {
                bool is255   = pClus->isSet(Event::TkrCluster::mask255);
                bool isGhost = pClus->isSet(Event::TkrCluster::maskGHOST);
                bool isAlone = pClus->isSet(Event::TkrCluster::maskALONE);
                bool isAloneEnd = pClus->isSet(Event::TkrCluster::maskALONEEND);
                if(!is255&&!isGhost) {
                    pClus->setStatusBits(Event::TkrCluster::maskSAMETRACK);
                }
            } 
        }
    }
    return sc;
}

StatusCode TkrGhostTool::flagEarlyVertices()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Flag the hits on tracks containing ghosts or 255's
    //get the tracks
    SmartDataPtr<Event::TkrVertexCol> 
        vertexCol(m_dataSvc, EventModel::TkrRecon::TkrVertexCol);

    int vertexCount = 0;
    Event::TkrVertexConPtr tcolIter = vertexCol->begin();
    for(tcolIter; tcolIter!=vertexCol->end();++tcolIter,++vertexCount) {
        Event::TkrVertex* vertex = *tcolIter;
        vertex->clearStatusBits(Event::TkrVertex::GHOST);
        SmartRefVector<Event::TkrTrack>::const_iterator trackIter;
        SmartRefVector<Event::TkrTrack>::const_iterator trackBegin = vertex->getTrackIterBegin();
        SmartRefVector<Event::TkrTrack>::const_iterator trackEnd = vertex->getTrackIterEnd();

        for(trackIter=trackBegin;trackIter!=trackEnd;++trackIter) {
            const Event::TkrTrack* track = *(trackIter);
            if((track->getStatusBits())&Event::TkrTrack::GHOST) {
                vertex->setStatusBit(Event::TkrVertex::GHOST);
                continue;
            }
        }
    }
    return sc;
}
