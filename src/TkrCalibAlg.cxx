
//$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrCalibAlg.cxx,v 1.17 2010/07/22 03:23:20 lsrea Exp $

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDetDataSvc.h"

#include "CalibData/CalibTime.h"
#include "CalibSvc/ICalibRootSvc.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Tkr/BadStrips.h"
#include "CalibData/Tkr/TkrSplitsCalib.h"
#include "CalibData/Tkr/TkrTot.h"
#include "CalibData/Tkr/TkrScale.h"
#include "CalibData/Tkr/TkrTowerAlignCalib.h"
#include "CalibData/Tkr/TkrInternalAlignCalib.h"

#include "TkrUtil/ITkrBadStripsSvcCalib.h"
#include "TkrUtil/ITkrFailureModeSvcCalib.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "TkrUtil/ITkrToTSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

#include <string>

/**
@file TkrCalibAlg.cxx Algorithm to organize the updating the calibration data

TkrCalibAlg is an algorithm which
accesses calibration data in what is expected to be a standard manner to inform
the services of the appearance of a new calibration. 
*/

/** 
@class TkrCalibAlg
Algorithm that handles updating the calibrations
*/
class TkrCalibAlg : public Algorithm {

public:
    TkrCalibAlg(const std::string& name, ISvcLocator* pSvcLocator);    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:

    void showCalibrationInfo(const std::string type, const std::string path, 
        const CalibData::CalibBase* ptr, const CalibData::BadStrips* bs_ptr=0) const; 
    StatusCode failedAccess(const std::string type) const;
    StatusCode failedUpdate(const std::string type) const;

    /// pointer to data provider
    IDataProviderSvc*   m_pCalibDataSvc;
    /// pointer to path provider
    ICalibPathSvc*      m_pCalibPathSvc;
    /// Handle to the IDetDataSvc interface of the CalibDataSvc
    IDetDataSvc*        m_detDataSvc;
    /// pointer to bad strips service
    ITkrBadStripsSvcCalib*   m_pTkrBadStripsSvc;
    /// pointer to failure mode service
    ITkrFailureModeSvcCalib* m_pTkrFailureModeSvc;
    /// pointer to splits service
    ITkrSplitsSvc*           m_pTkrSplitsSvc;
    /// pointer to splits service
    ITkrToTSvc*              m_pTkrToTSvc;
    /// pointer to alignment service
    ITkrAlignmentSvc*        m_pTkrAlignmentSvc;
    /// serial number of the hot strips calibration
    int m_serHot;
    /// serial number of the dead strips calibration
    int m_serDead;
    /// serial number of the splits calibration
    int m_serSplits;
    /// serial number of the charge injection calibration
    int m_serInjection;
    /// serial number of the ToT muon calibration
    int m_serMuons;
    /// serial number of the tower alignment calibration
    int m_serTowerAlign;
    /// serial number of the internal alignment calibration
    int m_serInternalAlign;

    /// flavor of calibration required ("ideal" means: do nothing)
    std::string m_flavor;
    // these for individual calib types, overrides the above.
    /// flavor of the dead strips calibration
    std::string m_deadStripsFlavor;
    /// flavor of the hot strips calibration
    std::string m_hotStripsFlavor;
    /// flavor of the splits calibration
    std::string m_splitsFlavor;
    /// flavor of the charge injection calibration
    std::string m_injectionFlavor;
    /// flavor of the ToT muon calibration
    std::string m_muonFlavor;
    /// flavor of the tower alignment calibration
    std::string m_towerAlignFlavor;
    /// flavor of the internal alignment calibration
    std::string m_internalAlignFlavor;

};

/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<TkrCalibAlg> Factory;
//const IAlgFactory& TkrCalibAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TkrCalibAlg);

TkrCalibAlg::TkrCalibAlg(const std::string&  name, 
                         ISvcLocator*        pSvcLocator )
                         : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
                         m_pTkrBadStripsSvc(0), m_pTkrFailureModeSvc(0),
                         m_serHot(-1), m_serDead(-1)
{
    declareProperty("calibFlavor",           m_flavor           = "ideal");
    declareProperty("deadStripsCalibFlavor", m_deadStripsFlavor = "notSet");
    declareProperty("hotStripsCalibFlavor",  m_hotStripsFlavor  = "notSet");
    declareProperty("splitsCalibFlavor",     m_splitsFlavor     = "notSet");
    declareProperty("chargeInjectionCalibFlavor",  
        m_injectionFlavor  = "notSet");
    declareProperty("muonCalibFlavor",       m_muonFlavor       = "notSet");

    declareProperty("internalAlignmentCalibFlavor", m_internalAlignFlavor = "notSet");
    declareProperty("towerAlignmentCalibFlavor",    m_towerAlignFlavor    = "notSet");
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

    sc = service("CalibDataSvc", m_pCalibPathSvc, true); 
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get ICalibPathSvc interface of CalibDataSvc"
            << endreq;
        return sc;
    }

    // Query the IDetDataSvc interface of the calib data service
    sc = m_pCalibDataSvc->queryInterface(IDetDataSvc::interfaceID(), 
        (void**) &m_detDataSvc);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not query IDetDataSvc interface of CalibDataSvc" 
            << endreq;
        return sc;
    } else {
        log << MSG::DEBUG 
            << "Retrieved IDetDataSvc interface of CalibDataSvc" 
            << endreq;
    }

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

    sc = service("TkrSplitsSvc", m_pTkrSplitsSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrSplitsSvc" 
            << endreq;
        return sc;
    }

    sc = service("TkrToTSvc", m_pTkrToTSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrToTSvc" 
            << endreq;
        return sc;
    }

    sc = service("TkrAlignmentSvc", m_pTkrAlignmentSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrAlignmentSvc" 
            << endreq;
        return sc;
    }

    // go through the individual flavors... 
    // set them equal to the overall flavor unless they've been set

    if (m_deadStripsFlavor=="notSet") m_deadStripsFlavor = m_flavor;
    if (m_hotStripsFlavor=="notSet")  m_hotStripsFlavor  = m_flavor;
    if (m_splitsFlavor=="notSet")     m_splitsFlavor     = m_flavor;
    if (m_injectionFlavor=="notSet")  m_injectionFlavor  = m_flavor;
    if (m_muonFlavor=="notSet")       m_muonFlavor       = m_flavor;
    if (m_internalAlignFlavor=="notSet")  m_internalAlignFlavor = m_flavor;
    if (m_towerAlignFlavor=="notSet")     m_towerAlignFlavor    = m_flavor;

    m_serHot       = -1;
    m_serDead      = -1;
    m_serSplits    = -1;
    m_serInjection = -1;
    m_serMuons     = -1;
    m_serInternalAlign = -1;
    m_serTowerAlign    = -1;

    return StatusCode::SUCCESS;   
}


StatusCode TkrCalibAlg::execute( ) {

    MsgStream log(msgSvc(), name());

    // set service for SIM or REC

    // look for "Rec" in the name
    if(name().find("Rec")!=std::string::npos) {
        m_pTkrBadStripsSvc->SetCalibType(ITkrBadStripsSvcCalib::REC);
        m_pTkrFailureModeSvc->SetCalibType(ITkrFailureModeSvcCalib::REC);
        m_pTkrAlignmentSvc->SetCalibType(ITkrAlignmentSvc::REC);
    } else {
        m_pTkrBadStripsSvc->SetCalibType(ITkrBadStripsSvcCalib::SIM);
        m_pTkrFailureModeSvc->SetCalibType(ITkrFailureModeSvcCalib::SIM);
        m_pTkrAlignmentSvc->SetCalibType(ITkrAlignmentSvc::SIM);
    }

    std::string fullPath;
    DataObject* pObject;
    std::string type;

    // check the dead channels
    type = "dead strips";
    CalibData::BadStrips* pDead = 0;

    if(m_deadStripsFlavor!="ideal" && m_deadStripsFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_DeadChan, m_deadStripsFlavor );
        pDead = SmartDataPtr<CalibData::BadStrips>(m_pCalibDataSvc, fullPath);
        if (!pDead) {return failedAccess(type);}

        m_pCalibDataSvc->updateObject((CalibData::BadStrips *)pDead);
        if (!pDead) {return failedAccess(type);}

        int newSerNo = pDead->getSerNo();
        if (newSerNo!=m_serDead) {
            m_serDead = newSerNo;
            showCalibrationInfo(type, fullPath, pDead, pDead); 
            m_pTkrBadStripsSvc->update(pDead, 0);
            m_pTkrFailureModeSvc->update(pDead, 0);
        }
    }

    // now the hot channels
    type = "hot strips";
    CalibData::BadStrips* pHot = 0;

    if(m_hotStripsFlavor!="ideal" && m_hotStripsFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_HotChan, m_hotStripsFlavor );
        pHot = SmartDataPtr<CalibData::BadStrips>(m_pCalibDataSvc, fullPath);
        if (!pHot) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject((CalibData::BadStrips *)pHot);
        if (!pHot) { return failedUpdate(type); }

        int newSerNo = pHot->getSerNo();
        if (newSerNo!=m_serHot) {
            m_serHot = newSerNo;
            showCalibrationInfo(type, fullPath, pHot, pHot); 
            m_pTkrBadStripsSvc->update(0, pHot);
            m_pTkrFailureModeSvc->update(0, pHot);
        }
    }

    // for the rest, we need to update on every call, because we may be alternating calibrations

    // next the splits
    type = "splits";
    CalibData::TkrSplitsCalib* pSplits = 0;

    if(m_splitsFlavor!="ideal" && m_splitsFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_Splits, m_splitsFlavor );
        m_pCalibDataSvc->retrieveObject(fullPath, pObject);
        pSplits = dynamic_cast<CalibData::TkrSplitsCalib*> (pObject);
        if (!pSplits) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject(pObject);
        pSplits = dynamic_cast<CalibData::TkrSplitsCalib*> (pObject);
        if (!pSplits) { return failedUpdate(type); }

        int newSerNo = pSplits->getSerNo();
        if (newSerNo!=m_serSplits) {
            m_serSplits = newSerNo;
            showCalibrationInfo(type, fullPath, pSplits); 
        }
    }
    m_pTkrSplitsSvc->update(pSplits);

    // now charge injection
    type = "charge injection";
    CalibData::TkrTotCol* pToT = 0;

    if(m_injectionFlavor!="ideal" && m_injectionFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_TOTSignal, m_injectionFlavor );
        m_pCalibDataSvc->retrieveObject(fullPath, pObject);
        pToT = dynamic_cast<CalibData::TkrTotCol*> (pObject);
        if (!pToT) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject(pObject);
        pToT = dynamic_cast<CalibData::TkrTotCol*> (pObject);
        if (!pToT) { return failedUpdate(type); }

        int newSerNo = pToT->getSerNo();        
        if (newSerNo!=m_serInjection) {
            m_serInjection = newSerNo;
            showCalibrationInfo(type, fullPath, pToT);
        }
    }
    m_pTkrToTSvc->update(pToT);

    // now muon calibration
    type = "muon scale";
    CalibData::TkrScaleCol* pScale = 0;

    if(m_muonFlavor!="ideal" && m_muonFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_ChargeScale, m_muonFlavor );
        m_pCalibDataSvc->retrieveObject(fullPath, pObject);
        pScale = dynamic_cast<CalibData::TkrScaleCol*> (pObject);
        if (!pScale) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject(pObject);
        pScale = dynamic_cast<CalibData::TkrScaleCol*> (pObject);
        if (!pScale) { return failedUpdate(type); }

        int newSerNo = pScale->getSerNo();
        if (newSerNo!=m_serMuons) {
            m_serMuons = newSerNo;
            showCalibrationInfo(type, fullPath, pScale);
        }
    }
    m_pTkrToTSvc->update(pScale);

    // now tower alignment calibration
    type = "tower alignment";
    bool updateAlign = false;
    CalibData::TkrTowerAlignCalib* pTowerAlign = 0;

    if(m_towerAlignFlavor!="ideal" && m_towerAlignFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_TowerAlign, m_towerAlignFlavor );
        m_pCalibDataSvc->retrieveObject(fullPath, pObject);
        pTowerAlign = dynamic_cast<CalibData::TkrTowerAlignCalib*> (pObject);
        if (!pTowerAlign) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject(pObject);
        pTowerAlign = dynamic_cast<CalibData::TkrTowerAlignCalib*> (pObject);
        if (!pTowerAlign) { return failedUpdate(type); }

        int newSerNo = pTowerAlign->getSerNo();
        if (newSerNo!=m_serTowerAlign) {
            m_serTowerAlign = newSerNo;
            showCalibrationInfo(type, fullPath, pTowerAlign);
            updateAlign = true;
        }
    }

    // now internal alignment calibration
    type = "internal alignment";
    CalibData::TkrInternalAlignCalib* pInternalAlign = 0;

    if(m_internalAlignFlavor!="ideal" && m_internalAlignFlavor!="") {

        fullPath = m_pCalibPathSvc->getCalibPath(
            ICalibPathSvc::Calib_TKR_InternalAlign, m_internalAlignFlavor );
        m_pCalibDataSvc->retrieveObject(fullPath, pObject);
        pInternalAlign = dynamic_cast<CalibData::TkrInternalAlignCalib*> (pObject);
        if (!pInternalAlign) { return failedAccess(type); }

        m_pCalibDataSvc->updateObject(pObject);
        pInternalAlign = dynamic_cast<CalibData::TkrInternalAlignCalib*> (pObject);
        if (!pInternalAlign) { return failedUpdate(type); }

        int newSerNo = pInternalAlign->getSerNo();
        if (newSerNo!=m_serInternalAlign) {
            m_serInternalAlign = newSerNo;
            showCalibrationInfo(type, fullPath, pInternalAlign);
            updateAlign = true;
        }
    }

    if(updateAlign) m_pTkrAlignmentSvc->update(pTowerAlign, pInternalAlign);

    return StatusCode::SUCCESS;
}

StatusCode TkrCalibAlg::finalize( ) {

    MsgStream log(msgSvc(), name());
    log << MSG::INFO 
        << "          Finalize TkrCalibAlg "
        << endreq;

    return StatusCode::SUCCESS;
}

void TkrCalibAlg::showCalibrationInfo(const std::string type, 
                                      const std::string path, 
                                      const CalibData::CalibBase* ptr,
                                      const CalibData::BadStrips* bs_ptr) const
{
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "New " << type << " serial number: " << ptr->getSerNo()<< endreq;  
    log << "path: " << path << endreq;
    log << "Vstart: " <<  (ptr->validSince()).hour(true)
        << "  Vend: " << (ptr->validTill()).hour(true) << endreq;
    if(bs_ptr!=0) {
        log << "Bad type: " << bs_ptr->getBadType() 
            << " has " << bs_ptr->getBadTowerCount() << " towers with " << type << endreq;				
    }
}

StatusCode TkrCalibAlg::failedAccess(const std::string type) const
{
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR 
        << "Failed to access " << type << endreq;
    return StatusCode::FAILURE;
}

StatusCode TkrCalibAlg::failedUpdate(const std::string type) const
{
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR 
        << "Failed to update " << type << endreq;
    return StatusCode::FAILURE;
}
