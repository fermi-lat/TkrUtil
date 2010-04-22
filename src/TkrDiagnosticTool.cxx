/**
* @class TkrDiagnosticTool
*
* @brief Implements a Gaudi Tool for setting the candidate track energies before 
*        the track fit
*
* @author The Tracking Software Group
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrDiagnosticTool.cxx,v 1.2 2010/04/08 20:54:04 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/TkrDigi.h"

#include "LdfEvent/Gem.h"
#include "LdfEvent/DiagnosticData.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrDiagnosticTool.h"
#include "TkrUtil/ITkrMapTool.h"

#include <iomanip>
#include <map>

class TkrDiagnosticTool : public AlgTool, virtual public ITkrDiagnosticTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrDiagnosticTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrDiagnosticTool() {}

    /// @brief calculations associated with the TEM diagnostic info

    StatusCode initialize();

    StatusCode getTkrVector(unsigned short& tkrVector);
    StatusCode getTkrDiagnosticData();
    StatusCode calculateTkrVector(unsigned short& diagBits);
    void setTriggerInfo(int tower, int gtcc, int gtrc);
    bool isSetTrigger(int tower, int plane, int end);
    void clearTriggerInfo();


private:
    //int flagHitsFromDiag( Event::TkrClusterCol* clusterCol);
    void flagToT255Hits(Event::TkrClusterCol* clusterCol);
    //int elecToGeo(int gtcc, int gtrc);
    //int elecIndex(int gtcc, int gtrc) { return gtccMult*gtcc+gtrc; }

    Event::DigiEvent*      m_digiEvent;

    /// Pointer to the Gaudi data provider service
    IGlastDetSvc*          m_pDetSvc;
    ITkrGeometrySvc*       m_tkrGeom;
    ITkrMapTool*           m_mapTool;
    DataSvc*               m_dataSvc;
    int m_numTowers;
    int m_numLayers;
    std::vector<unsigned int> m_towerBits;

    towerVec m_diagTrigger;

    //std::map< int, int>  m_tkrMap;
    std::map< int, bool> m_triggerMap;
    std::vector<unsigned int> m_layerBits;
};

static ToolFactory<TkrDiagnosticTool> s_factory;
const IToolFactory& TkrDiagnosticToolFactory = s_factory;

TkrDiagnosticTool::TkrDiagnosticTool(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrDiagnosticTool>(this);

    return;
}

StatusCode TkrDiagnosticTool::initialize()
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

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("TkrMapTool", m_mapTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved TrkMapTool" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve TkrMapTool" << endreq;
        }

    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    } 


    int numX, numY;
    m_pDetSvc->getNumericConstByName("xNum", &numX);
    m_pDetSvc->getNumericConstByName("yNum", &numY);  
    m_numTowers = numX*numY;

    m_numLayers = m_tkrGeom->numLayers();

    m_towerBits.resize(m_numTowers,0);

    m_diagTrigger.resize(m_numTowers);
    int i;
    for(i=0; i<m_numTowers;++i) {
        m_diagTrigger[i] = new TkrTowerBits();
        m_diagTrigger[i]->setDebug(false);
    }

    return sc;
}

StatusCode TkrDiagnosticTool::getTkrVector(unsigned short& tkrVector)
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

StatusCode TkrDiagnosticTool::getTkrDiagnosticData()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    clearTriggerInfo();

    // Retrieve the Event data for this event
    SmartDataPtr<LdfEvent::DiagnosticData> diagTds(m_dataSvc, "/Event/Diagnostic");
    if (!diagTds) return StatusCode::FAILURE;

    int ind;
    int numTkrDiag = diagTds->getNumTkrDiagnostic();
    //int numCalDiag = diagTds->getNumCalDiagnostic();
    log << MSG::DEBUG;
    if(log.isActive()) {
        int nNonZero = 0;
        if (numTkrDiag>0) {
            for (ind = 0; ind < numTkrDiag; ind++) {
                LdfEvent::TkrDiagnosticData tkrDiagTds = diagTds->getTkrDiagnosticByIndex(ind);
                int dataword = tkrDiagTds.dataWord();
                if (dataword!=0) nNonZero++;
                std::cout << "TkrDiagData: " << ind << " " 
                    << tkrDiagTds.tower() << " " << tkrDiagTds.gtcc() << " " 
                    << dataword << std::endl;
            }
        }
        log << numTkrDiag 
            << " Tkr diagnostic records found, " 
            << nNonZero << " non-zero";
    }
    log << endreq;

    if(numTkrDiag==0) return StatusCode::FAILURE;

    for (ind = 0; ind < numTkrDiag; ind++) {
        LdfEvent::TkrDiagnosticData tkrDiagTds = diagTds->getTkrDiagnosticByIndex(ind);
        int dataword = tkrDiagTds.dataWord();
        if (dataword==0) continue;
        log << MSG::VERBOSE;
        if(log.isActive()) {
        log << ind << " " << "datum/gccc/gcrc " 
            << tkrDiagTds.tower() << " " 
            << tkrDiagTds.gtcc() << " " 
            << std::oct << tkrDiagTds.dataWord()
            << std::dec;
        }
        log << endreq;

        int tower = tkrDiagTds.tower();
        int cc = tkrDiagTds.gtcc();
        int rcBits = tkrDiagTds.dataWord();
        if(rcBits>0) {
            int rc;
            for(rc=0;rc<nRc;++rc) {
                if((rcBits&1)==1) setTriggerInfo(tower, cc, rc);
                rcBits >>= 1;
            }
        }
    }

    return sc;
}

//int TkrDiagnosticTool::elecToGeo(int gtcc, int gtrc)
//{
//    return m_tkrMap[elecIndex(gtcc, gtrc)];
//}

void TkrDiagnosticTool::setTriggerInfo(int tower, int gtcc, int gtrc)
{
    int geoIndex = towerMult*tower + planeMult*m_mapTool->elecToGeo(gtcc, gtrc) + endArray[gtcc];
    m_triggerMap[geoIndex] = true;
}

void TkrDiagnosticTool::clearTriggerInfo() {
    m_triggerMap.clear();
    //std::cout << "size of triggerMap after clear " 
    //    << m_triggerMap.size() << std::endl;
}


StatusCode TkrDiagnosticTool::calculateTkrVector(unsigned short &diagBits)
{
    
   MsgStream log(msgSvc(), name());

   //set up the tower bits vector
    diagBits = 0;
    int tower, plane, end;
    //std::cout << "num Towers " << m_numTowers << std::endl;
    for(tower=0;tower<m_numTowers;++tower) {
        m_diagTrigger[tower]->clear();
        for(plane=0;plane<m_tkrGeom->numPlanes();++plane) {
            for(end=0;end<2;++end) {
                if(isSetTrigger(tower, plane, end)) {
                    int layer = m_tkrGeom->getLayer(plane);
                    int view  = m_tkrGeom->getView(plane);
                    m_diagTrigger[tower]->setBit(layer, view);
                    log << MSG::VERBOSE;
                    if(log.isActive()) {
                        log << "diagTowerBits: " 
                            << tower << " " << plane << " " << end ;
                    }
                    log << endreq;
                }
            }
        }
        bool triggered = m_diagTrigger[tower]->isTriggered();
        int bits = m_diagTrigger[tower]->getBits();
        int trigBits = m_diagTrigger[tower]->getTriggeredBits();

        if(triggered) diagBits |= (1<<tower);
        if(triggered+bits+trigBits+diagBits==0) continue;
        log << MSG::VERBOSE;
        if(log.isActive()) {
            log << "tower " << tower << " trigrd/bits/trigbits " <<
                triggered << " " << bits << " " << trigBits << endreq;
            log << "diagVector " << diagBits;
        }
        log << endreq;
    }
    return StatusCode::SUCCESS;
}

//int TkrDiagnosticTool::flagHitsFromDiag(Event::TkrClusterCol* clusterCol)
//{
//    int clusSize = clusterCol->size();
//    int nDiagBits = 0;
//    int i;
//    int lastIndex = -1;
//    for (i=0;i<clusSize;++i) {
//        Event::TkrCluster* clus = (*clusterCol)[i];
//
//        idents::TkrId tkrid = clus->getTkrId();
//        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
//        int tower = twrid.id();
//        int plane = m_tkrGeom->getPlane(tkrid);
//        int layer = clus->getLayer();
//        int layerIndex = towerMult*tower + layerMult*layer;
//        int end   = clus->getEnd();
//        int index = towerMult*tower + planeMult*plane + end;
//        if(!isSetTrigger(tower, plane, end)) {
//            clus->setStatusBits(Event::TkrCluster::maskGHOST);
//            if(lastIndex!=index) {
//                nDiagBits++;
//                m_layerBits[layerIndex]; 
//                lastIndex = index;
//                std::cout << "Diagbit: t/p/e = " << tower << " " 
//                    << plane << " " << end << std::endl; 
//            }
//        }
//    }
//    return nDiagBits;
//}

bool TkrDiagnosticTool::isSetTrigger(int tower, int plane, int end)
{
    int index = towerMult*tower + planeMult*plane + end;
    if(m_triggerMap.find(index)!=m_triggerMap.end()) return true;
    return false;
}
