
//$Header:$

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/BadStrips.h"

#include "TkrUtil/ITkrBadStripsSvcCalib.h"
#include "TkrUtil/ITkrFailureModeSvcCalib.h"

#include <string>


/**
@file TkrCalibAlg.cxx Simple algorithm organize the updating the bad strips data

  TkrCalibAlg is an algorithm which
  accesses hot and dead tracker strip calibration data in what is
  expected to be a standard manner. 
  
*/


/** 
@class TkrCalibAlg

  Algorithm that handles updating the bad strips calibrations
*/
class TkrCalibAlg : public Algorithm {
    
    
public:
    TkrCalibAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    
    StatusCode execute();
    
    StatusCode finalize();
    
private:
    
    /// pointer to data provider
    IDataProviderSvc*   m_pCalibDataSvc;
    /// pointer to bad strips service
    ITkrBadStripsSvcCalib*   m_pTkrBadStripsSvc;
    /// pointer to failure mode service
    ITkrFailureModeSvcCalib* m_pTkrFailureModeSvc;
    /// serial number of the hot strips calibration
    int                 m_serHot;
    /// serial number of the bad strips calibration
    int                 m_serDead;
    /// flavor of calibration required ("ideal" means: do nothing)
    std::string m_flavor;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<TkrCalibAlg> Factory;
const IAlgFactory& TkrCalibAlgFactory = Factory;


TkrCalibAlg::TkrCalibAlg(const std::string&  name, 
                         ISvcLocator*        pSvcLocator )
                         : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
                         m_pTkrBadStripsSvc(0), m_pTkrFailureModeSvc(0),
                         m_serHot(-1), m_serDead(-1)
{
    declareProperty("flavor", m_flavor="ideal");
}


StatusCode TkrCalibAlg::initialize() 
{
    StatusCode sc;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Initialize()" << endreq;
    
    // So far don't have any properties, but in case we do some day..
    setProperties();
    
    sc = service("CalibDataSvc", m_pCalibDataSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get IDataProviderSvc interface of CalibDataSvc" 
            << endreq;
        return sc;
    }

    unsigned long foo = 0;
    
    sc = service("TkrBadStripsSvc", m_pTkrBadStripsSvc, true);
    if ( !sc.isSuccess()) {
        log << MSG::ERROR 
            << "Could not get TkrBadStripsSvc" 
            << endreq;
        return sc;
    }
    
    sc = service("TkrFailureModeSvc", m_pTkrFailureModeSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrFailureModeSvc" 
            << endreq;
        return sc;
    }
        
    // Get properties from the JobOptionsSvc
    sc = setProperties();
    return StatusCode::SUCCESS;   
}


StatusCode TkrCalibAlg::execute( ) {
    
    MsgStream log(msgSvc(), name());

    if (m_flavor=="ideal") return StatusCode::SUCCESS;
    
    bool updateNow = false;

    // check the dead channels
    std::string fullDeadPath = "/Calib/TKR_DeadChan/vanilla";
    SmartDataPtr<CalibData::BadStrips> pDead(m_pCalibDataSvc, fullDeadPath);
    if (!pDead) {
        log << MSG::ERROR 
            << "Failed access to Dead strips via smart ptr" << endreq;
        return StatusCode::FAILURE;
    }

    m_pCalibDataSvc->updateObject((CalibData::BadStrips *)pDead);
    if (!pDead) {
        log << MSG::ERROR 
            << "Update of dead strips failed" << endreq;
        return StatusCode::FAILURE;
    }

    int newSerNo = pDead->getSerNo();
    if (newSerNo!=m_serDead) {
        log << MSG::INFO << "deadStrips serial number changed..." 
            << endreq;
        m_serDead = newSerNo;
        updateNow = true;
		log << MSG::INFO << "Retrieved with path " << fullDeadPath << endreq
			<< "Serial #" <<  pDead->getSerNo() << endreq; 
		log << MSG::INFO << "Vstart: " <<  (pDead->validSince()).hours()
			<< "  Vend: " << (pDead->validTill()).hours() << endreq;
		log << MSG::INFO << "Bad type: " << pDead->getBadType() 
			<< " has " << pDead->getBadTowerCount() << " bad towers " << endreq;				
    }

    // now the hot channels
    std::string fullHotPath = "/Calib/TKR_HotChan/vanilla";
    SmartDataPtr<CalibData::BadStrips> pHot(m_pCalibDataSvc, fullHotPath);
    if (!pHot) {
        log << MSG::ERROR 
            << "Failed access to Hot strips via smart ptr" << endreq;
        return StatusCode::FAILURE;
    }

    m_pCalibDataSvc->updateObject((CalibData::BadStrips *)pHot);
    if (!pHot) {
        log << MSG::ERROR 
            << "Update of Hot strips failed" << endreq;
        return StatusCode::FAILURE;
    }

    newSerNo = pHot->getSerNo();
    if (newSerNo!=m_serHot) {
        log << MSG::INFO << "hotStrips serial number changed..." 
            << endreq;
        m_serHot = newSerNo;
        updateNow = true;
		log << MSG::INFO << "Retrieved with path " << fullHotPath << endreq
			<< "Serial #" <<  pHot->getSerNo() << endreq; 
		log << MSG::INFO << "Vstart: " <<  (pHot->validSince()).hours()
			<< "  Vend: " << (pHot->validTill()).hours() << endreq;
		log << MSG::INFO << "Bad type: " << pHot->getBadType() 
			<< " has " << pHot->getBadTowerCount() << " bad towers " << endreq;
    }

    // just update everything if anything changes
    // no harm
    if (updateNow) {
        log << MSG::INFO <<" about to update constants" << endreq;
        m_pTkrBadStripsSvc->update(pDead, pHot);
        m_pTkrFailureModeSvc->update(pDead, pHot);
    }

    return StatusCode::SUCCESS;
}


StatusCode TkrCalibAlg::finalize( ) {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO 
        << "          Finalize TkrCalibAlg "
        << endreq;
    
    return StatusCode::SUCCESS;
}
