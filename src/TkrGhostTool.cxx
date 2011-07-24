/**
* @class TkrGhostTool
*
* @brief Tool for setting the ghost attributes of hits and tracks
*
* @author The Tracking Software Group
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrGhostTool.cxx,v 1.12 2011/06/09 04:09:02 heather Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/TkrDigi.h"

#include "LdfEvent/Gem.h"
#include "LdfEvent/DiagnosticData.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGhostTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrDiagnosticTool.h"
#include "TkrUtil/ITkrMapTool.h"
#include "TkrUtil/ITkrTrackVecTool.h"

#include "Doca.h"

#include <iomanip>
#include <map>

class TkrGhostTool : public AlgTool, virtual public ITkrGhostTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrGhostTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrGhostTool() {}

    /// @brief 

    StatusCode initialize();

    StatusCode getTkrVector(unsigned short& tkrVector);
    StatusCode getTkrGemDataCondArrivalTime(unsigned short& tkrTime);
    StatusCode calculateTkrVector(
        Event::TkrClusterCol* pCol, unsigned short& towerBits);
    StatusCode calculateTkrVector(
        Event::TkrDigiCol* pCol, unsigned short& towerBits);
    StatusCode flagSingles();
    StatusCode flagEarlyHits(Event::TkrClusterCol* col=0);
    StatusCode flagEarlyTracks();
    StatusCode flagEarlyVertices();
    StatusCode flagEarlyCalClusters();

private:
    int elecToGeo(int gtcc, int gtrc);
    int flagHitsFromDiag( Event::TkrClusterCol* clusterCol);
    void flagToT255Hits(Event::TkrClusterCol* clusterCol);

    Event::DigiEvent*      m_digiEvent;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
    IGlastDetSvc*          m_pDetSvc;
    ITkrGeometrySvc*       m_tkrGeom;
    ITkrDiagnosticTool*    m_pDiagTool;
    ITkrTrackVecTool*       m_trackVecTool;

    int m_numTowers;
    int m_numLayers;
    std::vector<unsigned int> m_towerBits;

    std::map< int, int>  m_tkrMap;
    std::map< int, bool> m_triggerMap;
    std::vector<unsigned int> m_diagBits;

    towerVec m_clusterTrigger;
    towerVec m_digiTrigger;

    bool m_useDiagInfo;
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
    
    declareProperty("UseDiagnosticInfo", m_useDiagInfo=true);

    return;
}

StatusCode TkrGhostTool::initialize()
{
    // Purpose and Method: 
    // Inputs:  
    // Outputs: 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    setProperties();

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

    m_pDiagTool = 0;
    if (toolSvc()->retrieveTool("TkrDiagnosticTool", m_pDiagTool).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrDiagnosticTool" << endreq;
        return StatusCode::FAILURE;
    }

    m_trackVecTool = 0;
    if (toolSvc()->retrieveTool("TkrTrackVecTool", m_trackVecTool).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrTrackVecTool" << endreq;
        return StatusCode::FAILURE;
    }

    int numX, numY;
    m_pDetSvc->getNumericConstByName("xNum", &numX);
    m_pDetSvc->getNumericConstByName("yNum", &numY);  
    m_numTowers = numX*numY;

    m_numLayers = m_tkrGeom->numLayers();

    m_towerBits.resize(m_numTowers,0);
    m_diagBits.resize(m_numTowers,0);

        //set up the tower bits vector
    m_clusterTrigger.resize(m_numTowers);
    m_digiTrigger.resize(m_numTowers);

    int i;
    for(i=0;i<m_numTowers;++i) {
        m_clusterTrigger[i] = new TkrTowerBits();
        m_digiTrigger[i]    = new TkrTowerBits();
    }


    return sc;
}

StatusCode TkrGhostTool::getTkrVector(unsigned short& tkrVector)
{
    MsgStream log(msgSvc(), name());

    tkrVector = 0;

    //get the tkrVector
    SmartDataPtr<LdfEvent::Gem> gem(m_dataSvc, "/Event/Gem");

    if (!gem) {
        log << MSG::DEBUG << "No GEM found on TDS" << endreq;
        return StatusCode::FAILURE;
    }

    tkrVector= gem->tkrVector();
    return StatusCode::SUCCESS;
}

StatusCode TkrGhostTool::getTkrGemDataCondArrivalTime(unsigned short& tkrTime)
{
    MsgStream log(msgSvc(), name());

    tkrTime = 0;

    //get the GemDataCondArrivalTime
    SmartDataPtr<LdfEvent::Gem> 
        gem(m_dataSvc, "/Event/Gem");

    if (!gem) {
        log << MSG::DEBUG << "No Gem found on TDS" << endreq;
        return StatusCode::FAILURE;
    }

    tkrTime = gem->condArrTime().tkr();
    return StatusCode::SUCCESS;
}

StatusCode TkrGhostTool::calculateTkrVector( 
    Event::TkrClusterCol* clusterCol,
    unsigned short& trigBits)
{
    //set up the tower bits vector
    int i;
    for(i=0;i<m_numTowers;++i) {
        m_clusterTrigger[i]->clear();
    }

    //fill the bits vector
    int clusSize = clusterCol->size();
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        int tower = clus->tower();
        m_clusterTrigger[tower]->setBit(clus);
    }

    //generate the tower trigger word
    trigBits = 0;
    for(i=0;i<m_numTowers;++i) {
        m_towerBits[i] = m_clusterTrigger[i]->getTriggeredBits();
        if(m_towerBits[i]>0) {trigBits |= (1<<i);}
    }
    return StatusCode::SUCCESS;
}

StatusCode TkrGhostTool::calculateTkrVector(
    Event::TkrDigiCol* digiCol,
    unsigned short& trigBits)
{
    //set up the tower bits vector
    int i;
    for(i=0;i<m_numTowers;++i) {
        m_digiTrigger[i]->clear();
    }

    trigBits = 0;

    //fill the bits vector
    int digiSize = digiCol->size();
    for (i=0;i<digiSize;++i) {
        Event::TkrDigi* digi = (*digiCol)[i];
        if (digi->getNumHits()==0) continue;
        int tower = (digi->getTower()).id();
        m_digiTrigger[tower]->setBit(digi);
    }

    //generate the tower trigger word
    trigBits = 0;
    for(i=0;i<m_numTowers;++i) {
        m_towerBits[i] = m_digiTrigger[i]->getTriggeredBits();
        if(m_towerBits[i]>0) {trigBits |= (1<<i);}
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

    log << MSG::DEBUG << "****************New Event****************" << endreq;

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
        int index = towerMult*tower + layerMult*layer + viewMult*view ;

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
        int index = towerMult*tower + layerMult*layer + viewMult*view ;

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
    StatusCode sc;

    //get the clusters
    if(clusterCol==0) {
        SmartDataPtr<Event::TkrClusterCol> 
            pClusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);
        clusterCol = pClusterCol;
    }

    // clear the bits to start
    int clusSize = clusterCol->size();
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        clus->clearStatusBits(Event::TkrCluster::maskZAPGHOSTS);
    }

    // do the 255's
    flagToT255Hits(clusterCol);

    // try to use the diagnostic info
    bool isDiagData = false;
    int nDiagBits = 0;
    unsigned short diagBits = 0;

    if(m_useDiagInfo) {
        sc = m_pDiagTool->getTkrDiagnosticData();
        if(sc.isSuccess()) {
            isDiagData = true;
            m_pDiagTool->calculateTkrVector(diagBits);
            nDiagBits = flagHitsFromDiag(clusterCol);
        }
    }


    // get the tkrVector
    unsigned short tkrVector;
    sc = getTkrVector(tkrVector);
    if(sc.isFailure()) return sc;

    unsigned short trigBits;
    sc = calculateTkrVector(clusterCol, trigBits);
    if(sc.isFailure()) return sc;

    if((tkrVector&trigBits)!=tkrVector) {
        log << MSG::DEBUG;
        if (log.isActive()) {
            log <<"tkrVector and calculated tower trigger disagree:" << endreq 
                << "tkrVector = " << std::oct << tkrVector 
            << ", calculation = " << trigBits << std::dec; 
        }
        log << endreq;
    }

    // Flag the hits that would have made the trigger
    // This is an assumption... for a three-in-a-row, any one
    //   missing hit would disable the trigger
    // Also, flag ToT==255

    clusSize = clusterCol->size();
    bool isGhost = false;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
        int tower = twrid.id();
        // we only look at towers with software trigger but no hardware trigger
        if((tkrVector&(1<<tower))!=0) continue;
        if((trigBits&(1<<tower))==0) continue;

        int layer = clus->getLayer();
        if(m_towerBits[tower]&(1<<layer)) {
            clus->setStatusBits(Event::TkrCluster::maskGHOST);
            isGhost = true;
            log << MSG::VERBOSE;
            if(log.isActive()) {
                int plane = clus->getPlane();
                int end   = clus->getEnd();
                int status = clus->getStatusWord();
                log << "cluster " 
                    << i << ", t/pl/end "  << tower << ", " << plane  << " " << end 
                    << " " << std::hex << status <<endreq;
                log << "tower bits " << std::oct << m_towerBits[tower] 
                << " "  << (1<<layer) << std::dec << endreq;
                log << "TowerBits: t/l = " << tower << " " << layer;
            }
            log << endreq;
        }
    }
    log << MSG::DEBUG;
    if(log.isActive()) {
        unsigned short int tkrTime;
        sc = getTkrGemDataCondArrivalTime(tkrTime);
        log << "Ghost check: ghost = " << isGhost << ", nDiagBits = " << nDiagBits;
        log << std::oct << ", tkrVector = " << tkrVector 
            << ", diagBits = " << diagBits << std::dec;
        log << " arrival time " << tkrTime;
        if(isDiagData&&(diagBits!=tkrVector)) log << " ***";
    }
    log << endreq;


    log << MSG::VERBOSE;
    if(log.isActive()) {
        for (i=0;i<clusSize;++i) {
            Event::TkrCluster* clus = (*clusterCol)[i];
            idents::TkrId tkrid = clus->getTkrId();
            idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
            int tower = twrid.id();
            int plane = clus->getPlane();
            int layer = clus->getLayer();
            int end   = clus->getEnd();
            int status = clus->getStatusWord();
            log << "cluster " 
                << i << ", t/pl/end "  << tower << ", " << plane  << " " << end 
                << " " << std::hex << status <<endreq;
            log << "tower bits " << std::oct << m_towerBits[tower] 
            << " "  << (1<<layer) << std::dec << endreq;
        }
    }
    log << endreq;

    return sc;
}

void TkrGhostTool::flagToT255Hits(Event::TkrClusterCol* clusterCol)
{
    int clusSize = clusterCol->size();
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        if(clus->getRawToT()==255) clus->setStatusBits(Event::TkrCluster::mask255);
    }
}

int TkrGhostTool::flagHitsFromDiag(Event::TkrClusterCol* clusterCol)
{
    MsgStream log(msgSvc(), name());
 
    int clusSize = clusterCol->size();
    int nDiagBits = 0;
    int tower;

    int lastIndex = -1;
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];

        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
        tower = twrid.id();
        int plane = m_tkrGeom->getPlane(tkrid);
        int layer = clus->getLayer();
        int view  = tkrid.getView();
        int end   = clus->getEnd();
        int index = towerMult*tower + planeMult*plane + end;
        // If trigger bit is not set, hit is a ghost!
        if(!m_pDiagTool->isSetTrigger(tower, plane, end)) {
            clus->setStatusBits(Event::TkrCluster::maskDIAGNOSTIC);
            log << MSG::VERBOSE;
            if(log.isActive()) {
                log << "SetDiagGhosttBit for (t/p/e): " << tower << " " 
                    << plane << " " << layer << " " << view << " " << end  
                    << " isSet " << clus->isSet(Event::TkrCluster::maskDIAGNOSTIC);
            }
            log << endreq;
            if(lastIndex!=index) {
                nDiagBits++;
                lastIndex = index;
            }
        }
    }
    return nDiagBits;
}

StatusCode TkrGhostTool::flagEarlyTracks()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Flag the hits on tracks containing ghosts or 255's
    //get the tracks

    std::vector<Event::TkrTrack*> trackVec = m_trackVecTool->getTrackVec();

    //SmartDataPtr<Event::TkrTrackCol> 
    //    trackCol(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);

    int trackCount = 0;
    unsigned int itk;
    
    //Event::TkrTrackColConPtr tcolIter = trackCol->begin();
    //for(; tcolIter!=trackCol->end(); ++tcolIter,++trackCount) {
    for(itk=0; itk<trackVec.size(); ++itk, ++trackCount) {
        Event::TkrTrack* track = trackVec[itk];;
        track->clearStatusBits(Event::TkrTrack::GHOST);
        track->clearStatusBits(Event::TkrTrack::TRIGGHOST);
        track->clearStatusBits(Event::TkrTrack::DIAGNOSTIC);
        track->clearStatusBits(Event::TkrTrack::GHOST255);
        Event::TkrTrackHitVecItr pHit = track->begin();

        int ghostCount = 0;
        int _255Count = 0;
        int diagCount = 0;
        while(pHit != track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;
            bool is255       = pClus->isSet(Event::TkrCluster::mask255);
            bool isGhost     = pClus->isSet(Event::TkrCluster::maskGHOST);
            bool isDiagGhost = pClus->isSet(Event::TkrCluster::maskDIAGNOSTIC); 
            //bool isAlone     = pClus->isSet(Event::TkrCluster::maskALONE);
            bool isAloneEnd  = pClus->isSet(Event::TkrCluster::maskALONEEND);
            bool isIdentified255 = (is255&&isAloneEnd);
            if(isIdentified255) _255Count++;
            if(isGhost)       ghostCount++;
            if(isDiagGhost) {
                diagCount++;
            }
        }

        if(_255Count==0&&ghostCount==0&&diagCount==0) continue;

        // this is some kind of ghost track!
        track->setStatusBit(Event::TkrTrack::GHOST);

    	log << MSG::DEBUG;
        bool doDebug = (log.isActive());
        log << endreq;
        if(diagCount>0){
            track->setStatusBit(Event::TkrTrack::DIAGNOSTIC);
            if(doDebug) log << "Found " << diagCount << " diag ghosts in track # " << trackCount << endreq;
        }
                
        if(_255Count>0){
            track->setStatusBit(Event::TkrTrack::GHOST255);
            if(doDebug) log << "Found " << _255Count << " 255's in track # " 
                << trackCount << endreq;
        }
        if(ghostCount>0){
            track->setStatusBit(Event::TkrTrack::TRIGGHOST);
            if(doDebug) log  << "Found " << ghostCount 
                << " trigger ghosts in track # " << trackCount << endreq;
        }

        // flag the remaining clusters on this track
        // we know that there are already some ghost hits
        pHit = track->begin();
        while(pHit!=track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;

            bool is255       = pClus->isSet(Event::TkrCluster::mask255);
            bool isGhost     = pClus->isSet(Event::TkrCluster::maskGHOST);
            bool isDiagGhost = pClus->isSet(Event::TkrCluster::maskDIAGNOSTIC); 
            bool isAloneEnd  = pClus->isSet(Event::TkrCluster::maskALONEEND);
            bool isIdentified255 = (is255&&isAloneEnd);
            if(!isIdentified255&&!isGhost&&!isDiagGhost) {
                pClus->setStatusBits(Event::TkrCluster::maskSAMETRACK);
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

    //int vertexCount = 0;
    Event::TkrVertexConPtr tcolIter = vertexCol->begin();
    for(; tcolIter!=vertexCol->end(); ++tcolIter) {
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

StatusCode TkrGhostTool::flagEarlyCalClusters()
{
    // Compute the minimum DOCA from all the ghost tracks
    // to all the CalClusters.  Eventually flag the ghost 
    // clusters (prescription yet to be developed).
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get the Cal Clusters
    Event::CalClusterCol* calClusterCol = 
        SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

    // No Cal Clusters, don't do anything
    if (!calClusterCol) return sc;

    // Get all the tracks from the tracker.
    std::vector<Event::TkrTrack*> trackVec = m_trackVecTool->getTrackVec();
    unsigned int trackCount, itk;
    double minDoca;

    // Loop over the clusters an find the min DOCA to a ghost track
    Event::CalClusterCol::const_iterator cluster ;
    for ( cluster = calClusterCol->begin() ;          
          cluster != calClusterCol->end() ;          
          cluster++) {          
        
        trackCount = 0;
        minDoca = 999999.;

        for(itk=0; itk<trackVec.size(); ++itk, ++trackCount) {
            Event::TkrTrack* track = trackVec[itk];;
            if( (track->getStatusBits())&Event::TkrTrack::GHOST ) {
                // Find the distance between Cal "centroid" and ghost track axis
                Doca trackDoca(track->getInitialPosition(), 
                               track->getInitialDirection());
                float tmpDoca = (float)trackDoca.docaOfPoint((*cluster)->getPosition());
                if ( tmpDoca < minDoca ) minDoca = tmpDoca;
            }
        }
        // If minDoca, unchanged set it with an unphysical value
        if (minDoca > 900000) minDoca = -1.;
        (*cluster)->getMomParamsRef().setMinGhostDoca(minDoca);
        log << MSG::DEBUG << "Found CalCluster minGhostDoca of: " 
            << minDoca << endreq;
    }
    return sc;
}
