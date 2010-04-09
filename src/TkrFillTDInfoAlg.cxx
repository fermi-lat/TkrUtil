
/** 
* @class TkrFillTDInfoAlg
*
* @brief Algorithm to construct TkrClusterCol/TkrCluster
*
* Adapted from SiCluster of Jose Hernando. 
*
* Handles bad strips
*
* @author Tracy Usher, Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrFillTDInfoAlg.cxx,v 1.1 2010/04/08 20:54:04 lsrea Exp $
*/

#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include <vector>
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Digi/TkrDigi.h"
#include "TkrUtil/ITkrMapTool.h"

#include "LdfEvent/Gem.h"
#include "LdfEvent/DiagnosticData.h"

class TkrFillTDInfoAlg : public Algorithm

{
public:
    TkrFillTDInfoAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFillTDInfoAlg() {}
    /// Looks for the geometry service (required) and the bad strips service 
    /// (optional)
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:

    /// pointer to geometry service
    ITkrGeometrySvc*         m_tkrGeom;

    IDataProviderSvc*        m_dataSvc;

    /// pointer to Tkr digis
    Event::TkrDigiCol*        m_TkrDigiCol;
    /// pointer to generated TkrClusterCol

    ITkrMapTool*              m_mapTool;

    LdfEvent::DiagnosticData* m_diagTds;

    int  m_testMode;
    int  m_nTowers;                     

};

static const AlgFactory<TkrFillTDInfoAlg>  Factory;
const IAlgFactory& TkrFillTDInfoAlgFactory = Factory;

TkrFillTDInfoAlg::TkrFillTDInfoAlg(const std::string& name, 
                                   ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    // expert parameter
    // -1 -> do nothing
    // 0 -> normal operation
    // 1 -> all ghosts
    // 2 -> mix of normal and ghosts
    declareProperty("testMode" , m_testMode=0);
}

using namespace Event;

StatusCode TkrFillTDInfoAlg::initialize()
{

    // Purpose and Method:  initializes TkrFillTDInfoAlg
    // Inputs:  None
    // Outputs: TkrGeometrySvc will be created if not already present
    // Dependencies:
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "TkrFillTDInfoAlg Initialization";
    if( (sc=setProperties()).isFailure()) log << " didn't work!";
    log << endreq;


    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        //throw GaudiException("Service [EventDataSvc] not found", name(), sc);
        return StatusCode::FAILURE;
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    //Look for the geometry service
    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." 
            << endreq;
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

    m_TkrDigiCol = 0;

    m_nTowers = m_tkrGeom->numXTowers()*m_tkrGeom->numYTowers();

    return StatusCode::SUCCESS;
}

StatusCode TkrFillTDInfoAlg::execute()
{
    // Purpose and Method: Fills the bits in TkrDiagnosticData
    // Inputs:
    // Outputs: 
    // TDS Input: TkrDigiCol
    // TDS Output: LdfEvent::DiagnosticData including TkrDiagnosticData
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;   
    MsgStream log(msgSvc(), name());

    // skip all this if mode==-1
    if(m_testMode==-1) return StatusCode::SUCCESS;

    // Recover a pointer to the raw digi objects, if none, nothing to do
    m_TkrDigiCol = SmartDataPtr<TkrDigiCol>(eventSvc(), EventModel::Digi::TkrDigiCol);
    if(!m_TkrDigiCol) return StatusCode::SUCCESS;

    // Retrieve the Diagnostic data for this event
    // If it's already there, nothing to do

    SmartDataPtr<LdfEvent::DiagnosticData> m_diagTds(m_dataSvc, "/Event/Diagnostic");
    if (m_diagTds) {
        if(m_diagTds->getNumTkrDiagnostic()>0) {return StatusCode::SUCCESS;}
    }
    else {
        // Create the DiagnosticData TDS object
        m_diagTds = new LdfEvent::DiagnosticData();
    }

    int size = m_nTowers*nCc;
    // this is the order that I found the diag objects in real data... 
    // probably doesn't matter but might as well keep them that way
    int order[nCc] =       { 6, 3, 7, 2, 5, 0, 4, 1 };
    // this goes the other way, so we can fill the array in the correct order
    int inverseOrder[nCc] = { 5, 7, 3, 1, 6, 4, 0, 2 };
    int tower, i;

    // prepare the vector of TkrDiagnosticData objects, in the correct order
    // we will be filling in the dataWords later
    for (tower= 0; tower<m_nTowers; ++tower) {
        for (i=0; i<nCc; ++i) {
            LdfEvent::TkrDiagnosticData tkr = 
                LdfEvent::TkrDiagnosticData(0, tower, order[i]);
            m_diagTds->addTkrDiagnostic(tkr); 
        }
    }

    // here we collect the bits in the dataWords, based on the digis and the mode
    // and add them to the existing words
    // skip this step for testMode==1

    if(m_testMode!=1) {
        Event::TkrDigiCol::const_iterator ppDigi = m_TkrDigiCol->begin();
        int count = 0;
        for (; ppDigi!= m_TkrDigiCol->end(); ppDigi++) {
            // each digi contains the digitized hits from one layer of one tower
            //  and from both ends

            Event::TkrDigi* pDigi = *ppDigi;

            int layer = pDigi->getBilayer();
            int view  = pDigi->getView();
            int tower = (pDigi->getTower()).id();
            int tray, botTop;
            m_tkrGeom->layerToTray(layer, view, tray, botTop);
            int plane = m_tkrGeom->trayToPlane(tray, botTop);
            int lastStrip = pDigi->getLastController0Strip();

            int numHits = pDigi->getNumHits();
            if(numHits==0) continue;

            int indEnd;
            bool iEnd[2] = {false, false};
            // no controller0 hits, so they must be controller 1
            iEnd[ (lastStrip<0 ? 1 : 0) ] = true;
            // remaining case is "the last strip" > lastStrip, controller 1
            if(iEnd[0] && numHits>1 && pDigi->getHit(numHits-1)>lastStrip) iEnd[1] = true;

            // if there's a hit, OR the bit into the correct dataWord
            for(indEnd=0;indEnd<2;++indEnd) {
                if(iEnd[indEnd]){
                    //set up the bit mask
                    int gtcc, gtrc;
                    m_mapTool->geoToElec(plane, indEnd, gtcc, gtrc); 
                    int dataWord = (1<<gtrc) ;
                    // get the index for this digi
                    int index = tower*nCc + inverseOrder[gtcc];
                    // mode 0 -> set all the bits
                    // mode 1 -> set ~half of the bits
                    if(m_testMode==0 || (m_testMode==2&&count%2==0)) {
                        dataWord 
                            |= m_diagTds->getTkrDiagnosticByIndex(index).dataWord();
                        m_diagTds->setTkrDataWordByIndex(index, dataWord);
                    }
                    count++;
                }
            }
        }
    }

    // Register the object in the TDS
    sc = eventSvc()->registerObject("/Event/Diagnostic", m_diagTds);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register " << "/Event/Diagnostic" << endreq;
        return sc;
    }

    // some checks and debug stuff here
    int ind;
    SmartDataPtr<LdfEvent::DiagnosticData> diagTds2(m_dataSvc, "/Event/Diagnostic");
    if(!diagTds2) {
        log << MSG::ERROR << "Failed to acquire TkrDiagnostic data" << endreq;
        return StatusCode::SUCCESS;
    }
    int numTkrDiag = diagTds2->getNumTkrDiagnostic();
    log << MSG::DEBUG;
    if(log.isActive()) {
        int nNonZero = 0;
        if (numTkrDiag>0) {
            for (ind = 0; ind < numTkrDiag; ind++) {
                LdfEvent::TkrDiagnosticData tkrDiagTDS 
                    = diagTds2->getTkrDiagnosticByIndex(ind);
                int dataword = tkrDiagTDS.dataWord();
                if (dataword!=0) nNonZero++;
                log << ind << " " 
                    << tkrDiagTDS.tower() << " " << tkrDiagTDS.gtcc() << " " 
                    << tkrDiagTDS.dataWord() << endreq;
            }
        }
        log << numTkrDiag 
            << " Tkr diagnostic records found, " 
            << nNonZero << " non-zero";
    }
    log << endreq;

    return StatusCode::SUCCESS;
}

StatusCode TkrFillTDInfoAlg::finalize()
{   
    return StatusCode::SUCCESS;
}
