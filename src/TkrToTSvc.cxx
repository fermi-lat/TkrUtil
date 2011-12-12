/*
@file TkrToTSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrUtil/src/TkrToTSvc.cxx,v 1.19.618.1 2010/09/09 14:03:23 heather Exp $

*/

#include "src/TkrToTSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "idents/TowerId.h"

#include <algorithm>
#include <fstream>
#include "facilities/Util.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

// declare the service factories for the TkrToTSvc
//static SvcFactory<TkrToTSvc> a_factory;
//const ISvcFactory& TkrToTSvcFactory = a_factory; 
DECLARE_SERVICE_FACTORY(TkrToTSvc);

using namespace idents;

TkrToTSvc::TkrToTSvc(const std::string& name,ISvcLocator* svc) 
: Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // declare the properties


    // values adjusted 03-Jul-05, based on average for Tracker A and B
    declareProperty("defaultThreshold", m_defaultThreshold = 1.177);   
    declareProperty("defaultGain",      m_defaultGain      = 0.589); 
    declareProperty("defaultQuad",      m_defaultQuad     = 0.00490);
    declareProperty("defaultQuality",   m_defaultQuality   = 0.0);
    // use the MPV, not the mean, to agree with Hiro's calibration
    declareProperty("mevPerMip"       , m_mevPerMip        = 0.113);
    // wallet card says 4.667, but Hiro uses 5.0, so fix this to agree
    // otherwise the ToT's come out wrong!
    declareProperty("fCPerMip"        , m_fCPerMip         = 5.0);


    declareProperty("defaultMuonScale",      m_defaultMuonScale      = 1.0);
    declareProperty("countsPerMicrosecond",  m_countsPerMicrosecond  = 5.0);
    declareProperty("maxToT",                m_maxToT                = 250);
    declareProperty("useSingleTowerConsts" , m_useSingleTowerConsts  = false);
    declareProperty("baseTower",             m_baseTower             = 0);
    declareProperty("useDefaultIfMissing"  , m_useDefaultIfMissing   = false);
}

StatusCode  TkrToTSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
{
    if (IID_ITkrToTSvc == riid) {
        *ppvIF = dynamic_cast<ITkrToTSvc*> (this);
        return StatusCode::SUCCESS;
    } else {
        return Service::queryInterface (riid, ppvIF);
    }
}

const InterfaceID&  TkrToTSvc::type () const {
    return IID_ITkrToTSvc;
}

StatusCode TkrToTSvc::initialize () 
{
    // Purpose and Method: Initialize the lists of dead units
    //                     
    // Inputs: None        
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode  status = StatusCode::SUCCESS;

    // Open the message log
    MsgStream log( msgSvc(), name() );

    // Call super-class
    Service::initialize ();

    // set the calibration pointers to zero
    // if they become non-zero, that means that calib info has
    // been requested and found
    m_pToT   = 0;
    m_pScale = 0;

    m_tkrGeom = 0;
    if( service( "TkrGeometrySvc", m_tkrGeom, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't retrieve TkrGeometrySvc" << endreq;
        return StatusCode::FAILURE;
    }

    m_nTowers = m_tkrGeom->numXTowers()*m_tkrGeom->numYTowers();
    m_nLayers = m_tkrGeom->numLayers();
    m_nViews  = 2;
    m_nStrips = m_tkrGeom->ladderNStrips()*m_tkrGeom->nWaferAcross();

    m_callCount = 0;

    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }

    status = doInit();

    return status;
}

StatusCode TkrToTSvc::doInit()
{
    // Open the message log
    MsgStream log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

double TkrToTSvc::getGain(int tower, int layer, int view, int strip) const
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return -1.0; }
    if (m_pToT) {
        TkrId id = getId(tower, layer, view);
        const CalibData::TkrTotStrip* pInfo = m_pToT->getStripInfo(id, strip);
        if(pInfo!=0) {
            return  pInfo->getSlope();
        } else {
            if (m_useDefaultIfMissing) {return m_defaultGain;}
            else {return 0.0;}
        }
    } else {
        return m_defaultGain;
    }
}
double TkrToTSvc::getQuad(int tower, int layer, int view, int strip) const 
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return -1.0; }
    if (m_pToT) {
        TkrId id = getId(tower, layer, view);
        const CalibData::TkrTotStrip* pInfo = m_pToT->getStripInfo(id, strip);
        if(pInfo!=0) {
            return  pInfo->getQuad();
        } else {
            if (m_useDefaultIfMissing) {return m_defaultQuad;}
            else {return 0.0;}
        }
    } else {
        return m_defaultQuad;
    }
}
double TkrToTSvc::getThreshold(int tower, int layer, int view, int strip) const 
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return -1.0; }
    if (m_pToT) {
        TkrId id = getId(tower, layer, view);
        const CalibData::TkrTotStrip* pInfo = m_pToT->getStripInfo(id, strip);
        if(pInfo!=0) {
            return  pInfo->getIntercept();
        } else {
            if (m_useDefaultIfMissing) {return m_defaultThreshold;}
            else {return 0.0;}
        }
    } else {
        return m_defaultThreshold;
    }
}
double TkrToTSvc::getQuality(int tower, int layer, int view, int strip) const 
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return -1.0; }
    if (m_pToT) {
        TkrId id = getId(tower, layer, view);
        const CalibData::TkrTotStrip* pInfo = m_pToT->getStripInfo(id, strip);
        if(pInfo!=0) {
            return  pInfo->getChi2();
        } else {
            if (m_useDefaultIfMissing) {return m_defaultQuality;}
            else {return 0.0;}
        }
    } else {
        return m_defaultQuality;
    }
}

double TkrToTSvc::getMuonScale(int tower, int layer, int view, int strip) const 
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    double muonScale = -1.0;
    if (!valid(tower, layer, view, strip)) { return muonScale; }
    if (m_pScale) {
        TkrId id = getId(tower, layer, view);
        CalibData::TkrScaleObj pInfo = m_pScale->getStripInfo(id, strip);
        muonScale = static_cast<double>(pInfo.getScale());
        if(muonScale==0.0 && m_useDefaultIfMissing) muonScale = m_defaultMuonScale;
    } else {
        muonScale = m_defaultMuonScale;
    }
    return muonScale;
}

int TkrToTSvc::getRawToT(double eDep, int tower, int layer, int view, int strip) const
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return 0; }

    TkrId id = getId(tower, layer, view);
    double gain, quad, threshold, muonScale;
    getConsts(id, strip, threshold, gain, quad, muonScale);

    double charge    = eDep/m_mevPerMip*m_fCPerMip; // in fCs

    // consts are for ToT -> charge
    // Here is the inverse:

    if(gain==0 || muonScale==0) { return 0; }

    double term = (charge/muonScale-threshold)/gain;
    double test = quad/gain;
    double time;
    if (fabs(test)>1.e-6) {
        // just the quadratic formula
        time = 0.5*(-1.0 + sqrt(1.0 + 4.*test*term))/test;
    } else {
        // degenerate case
        time = term*(1.0 - test*term);
    }

    time *= m_countsPerMicrosecond;
    int iToT = static_cast<int> ( std::max( 0., time));
    return std::min(iToT, m_maxToT);
}

double TkrToTSvc::getCharge(double rawToT, int tower, int layer, int view, int strip) const
{
    if(m_useSingleTowerConsts) tower = m_baseTower;
    if (!valid(tower, layer, view, strip)) { return 0.0; }

    TkrId id = getId(tower, layer, view);
    double gain, quad, threshold, muonScale;
    getConsts(id, strip, threshold, gain, quad, muonScale);

    double time = rawToT/m_countsPerMicrosecond;
    double charge = muonScale*(threshold + time*(gain + time*quad));

    return charge;
}
idents::TkrId TkrToTSvc::getId(int tower, int layer, int view) const
{
    int tray, face;
    m_tkrGeom->layerToTray(layer, view, tray, face);
    TowerId towerId = TowerId(tower);
    return TkrId(towerId.ix(), towerId.iy(), tray, (face==TkrId::eTKRSiTop));
}


double TkrToTSvc::getMipsFromToT(double rawToT, 
                                 int tower, int layer, int view, int strip) const
{
    return getCharge(rawToT, tower, layer, view, strip)/getFCPerMip();
}

double TkrToTSvc::getMipsFromCharge(double charge) const
{
    return charge/getFCPerMip();
}

void TkrToTSvc::getConsts(idents::TkrId id, int strip, 
                          double& threshold, double& gain,
                          double& quad, double& muonScale) const 
{
    bool warned = false;
    if (m_pToT) {
        const CalibData::TkrTotStrip* pInfo = m_pToT->getStripInfo(id, strip);
        if (pInfo!=0) {
            gain = pInfo->getSlope();
            quad = pInfo->getQuad();
            threshold = pInfo->getIntercept();
        } else {
            if(m_useDefaultIfMissing) {
                gain = m_defaultGain;
                quad = m_defaultQuad;
                threshold = m_defaultThreshold;
            } else {
                gain = 0.0;
                quad = 0.0;
                threshold = 0.0;
            }
            doWarning("ToT consts", id, strip);
            m_callCount++;
            warned = true;
        }
    } else {
        gain      = m_defaultGain;
        quad      = m_defaultQuad;
        threshold = m_defaultThreshold;
    }

    if(m_pScale) {
        CalibData::TkrScaleObj pInfo1 = m_pScale->getStripInfo(id, strip);
        muonScale = static_cast<double>(pInfo1.getScale());
        if(muonScale==0.0) {
            doWarning("muon scale", id, strip);
            if(m_useDefaultIfMissing) {
                muonScale = m_defaultMuonScale;
            }
            if (!warned) m_callCount++;
        }
    } else {
        muonScale = m_defaultMuonScale;
    }

    return;
}

void TkrToTSvc::doWarning(const std::string type, const idents::TkrId id, int strip) const
{
    const int maxWarn = 20;
    if (m_callCount>=maxWarn) return;
    MsgStream log(msgSvc(), name());
    log << MSG::WARNING << "Info missing for " << type 
        << ", TkrId: TowerX/Y " << id.getTowerX() << "/" << id.getTowerY()
        << " Tray/botTop " << id.getTray() << "/" << id.getBotTop() 
        << " strip " << strip;
    log << (m_useDefaultIfMissing ? ", default used" : ", set to zero") << endreq;
    if (m_callCount==maxWarn-1) log << " Warning suppressed for the rest of this job" << endreq;
    return;
}
StatusCode TkrToTSvc::finalize() {

    MsgStream log(msgSvc(), name());
    if(m_callCount>0) {
        log << MSG::WARNING << m_callCount 
            << "Call(s) to getConsts() with missing info;";
        if (m_useDefaultIfMissing){ 
            log << " default values used";
        } else {
            log << " consts set to zero";
        }
        log << endreq;
    }

    return StatusCode::SUCCESS;
}
