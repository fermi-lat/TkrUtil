/**
* @class TkrDiagnosticTool
*
* @brief Implements a Gaudi Tool for setting the candidate track energies before 
*        the track fit
*
* @author The Tracking Software Group
*
* $Header$
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
    int elecToGeo(int gtcc, int gtrc);
    int elecIndex(int gtcc, int gtrc) { return gtccMult*gtcc+gtrc; }

    Event::DigiEvent*      m_digiEvent;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
    IGlastDetSvc*          m_pDetSvc;
    ITkrGeometrySvc*       m_tkrGeom;

    int m_numTowers;
    int m_numLayers;
    std::vector<unsigned int> m_towerBits;

    towerVec m_diagTrigger;

    std::map< int, int>  m_tkrMap;
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

    //This maps the electronic space into the geometrical one
    //Two cables map on to each plane, one at each end
    //  Might want to refactor all this some day
    //
    // m_tkrMap[gtccMult*CC+RC] = PP; (gtccMult==10)
    //   CC = gtcc
    //   RC = gtrc
    //   PP = plane

    /*
    m_tkrMap[ 0] = m_tkrMap[10] =  0;    m_tkrMap[20] = m_tkrMap[30] =  1;
    m_tkrMap[60] = m_tkrMap[70] =  2;    m_tkrMap[40] = m_tkrMap[50] =  3;
    m_tkrMap[ 1] = m_tkrMap[11] =  4;    m_tkrMap[21] = m_tkrMap[31] =  5;
    m_tkrMap[61] = m_tkrMap[71] =  6;    m_tkrMap[41] = m_tkrMap[51] =  7;
    m_tkrMap[ 2] = m_tkrMap[12] =  8;    m_tkrMap[22] = m_tkrMap[32] =  9;
    m_tkrMap[62] = m_tkrMap[72] = 10;    m_tkrMap[42] = m_tkrMap[52] = 11;
    m_tkrMap[13] = m_tkrMap[ 3] = 12;    m_tkrMap[23] = m_tkrMap[33] = 13;
    m_tkrMap[63] = m_tkrMap[73] = 14;    m_tkrMap[43] = m_tkrMap[53] = 15;
    m_tkrMap[14] = m_tkrMap[ 4] = 16;    m_tkrMap[24] = m_tkrMap[34] = 17;
    m_tkrMap[74] = m_tkrMap[64] = 18;    m_tkrMap[44] = m_tkrMap[54] = 19;
    m_tkrMap[15] = m_tkrMap[ 5] = 20;    m_tkrMap[25] = m_tkrMap[35] = 21;
    m_tkrMap[75] = m_tkrMap[65] = 22;    m_tkrMap[45] = m_tkrMap[55] = 23;
    m_tkrMap[ 6] = m_tkrMap[16] = 24;    m_tkrMap[26] = m_tkrMap[36] = 25;
    m_tkrMap[76] = m_tkrMap[66] = 26;    m_tkrMap[46] = m_tkrMap[56] = 27;
    m_tkrMap[17] = m_tkrMap[ 7] = 28;    m_tkrMap[37] = m_tkrMap[27] = 29;
    m_tkrMap[77] = m_tkrMap[67] = 30;    m_tkrMap[57] = m_tkrMap[47] = 31;
    m_tkrMap[ 8] = m_tkrMap[18] = 32;    m_tkrMap[28] = m_tkrMap[38] = 33;
    m_tkrMap[68] = m_tkrMap[78] = 34;    m_tkrMap[58] = m_tkrMap[48] = 35;
    */

    // this is the same as the above
    for(i=0;i<9;++i) {
        m_tkrMap[i]    = m_tkrMap[10+i] = 4*i;
        m_tkrMap[20+i] = m_tkrMap[30+i] = 4*i+1;
        m_tkrMap[60+i] = m_tkrMap[70+i] = 4*i+2;
        m_tkrMap[40+i] = m_tkrMap[50+i] = 4*i+3;
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
    int numCalDiag = diagTds->getNumCalDiagnostic();
    log << MSG::DEBUG;
    if(log.isActive()) {
        int nNonZero = 0;
        if (numTkrDiag>0) {
            for (ind = 0; ind < numTkrDiag; ind++) {
                LdfEvent::TkrDiagnosticData tkrDiagTds = diagTds->getTkrDiagnosticByIndex(ind);
                int dataword = tkrDiagTds.dataWord();
                if (dataword!=0) nNonZero++;
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

int TkrDiagnosticTool::elecToGeo(int gtcc, int gtrc)
{
    return m_tkrMap[elecIndex(gtcc, gtrc)];
}

void TkrDiagnosticTool::setTriggerInfo(int tower, int gtcc, int gtrc)
{
    int geoIndex = towerMult*tower + planeMult*elecToGeo(gtcc, gtrc) + endArray[gtcc];
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
