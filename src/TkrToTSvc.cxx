/*
@file TkrToTSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrToTSvc.cxx,v 1.7 2004/12/16 23:28:30 usher Exp $

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
static SvcFactory<TkrToTSvc> a_factory;
const ISvcFactory& TkrToTSvcFactory = a_factory; 

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

    declareProperty("ToTFile",          m_ToTFile          = ""  );
    // The default gains and threshold are chosen to reproduce the 
    // previous behavior of the ToT
    declareProperty("defaultGain",      m_defaultGain      = 2.50267833 );
    declareProperty("defaultGain2",     m_defaultGain2     = 0.0);
    declareProperty("defaultThreshold", m_defaultThreshold = -2.92);
    declareProperty("defaultQuality",   m_defaultQuality   = 0.0);
    declareProperty("defaultMuonFactor", m_defaultMuonFactor = 1.0);
    declareProperty("mode"            , m_mode             = "ideal");
    declareProperty("countsPerMicrosecond", m_countsPerMicrosecond = 5.0);
    declareProperty("mevPerMip"       , m_mevPerMip        = 0.155);
    declareProperty("fCPerMip"        , m_fCPerMip         = 4.667);
    declareProperty("maxToT"          , m_maxToT           = 250);
}

StatusCode  TkrToTSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrToTSvc == riid) {
        *ppvIF = dynamic_cast<ITkrToTSvc*> (this);
        return StatusCode::SUCCESS;
    } else {
        return Service::queryInterface (riid, ppvIF);
    }
}

const IID&  TkrToTSvc::type () const {
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
    
    m_tkrGeom = 0;
    if( service( "TkrGeometrySvc", m_tkrGeom, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't retrieve TkrGeometrySvc" << endreq;
        return StatusCode::FAILURE;
    }
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
        
    if (m_ToTFile!="") {
        int ret =  facilities::Util::expandEnvVar(&m_ToTFile);
        if (ret>=0) {
            log << MSG::INFO << "Input file for splits: " 
                << m_ToTFile << endreq;
        } else {
            log << MSG::ERROR << "Input filename " << m_ToTFile << " not resolved" << endreq;
            return StatusCode::FAILURE;
        }
    }

    status = doInit();
    
    return status;
}

StatusCode TkrToTSvc::doInit()
{
    // Open the message log
    MsgStream log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;

    // can be removed when geometry is iterfaced here

    const int nChipsPerLadder = m_tkrGeom->chipsPerLadder();
    const int nChips  = nChipsPerLadder*m_tkrGeom->nWaferAcross();
    const int nStrips = m_tkrGeom->ladderNStrips()/nChipsPerLadder;
    const int numLayers = m_tkrGeom->numLayers();

    int tower, layer, view, strip, chip;

    if(m_mode.substr(0,5)=="ideal") {
    // all gains and thresholds set to the same value, to reproduce the standard MC ToT
        for(tower=0;tower<NTOWERS;++tower) {
            for (layer=0;layer<numLayers;++layer) {
                for (view=0;view<NVIEWS;++view) {
                    for(chip=0;chip<nChips;++chip) {
                        for (strip=0;strip<nStrips;++strip) {
                            int theStrip = chip*nStrips + strip;
                            m_ToTGain[tower][layer][view][theStrip] = 
                                m_defaultGain;
                            m_ToTGain2[tower][layer][view][theStrip] = 
                                m_defaultGain2;
                            m_ToTThreshold[tower][layer][view][theStrip] = 
                                m_defaultThreshold;
                            m_ToTQuality  [tower][layer][view][theStrip] = 
                                m_defaultQuality;
                            m_ToTMuonFactor[tower][layer][view][theStrip] = 
                                m_defaultMuonFactor;
                        }
                    }
                }
            }
    }
    } else if (m_mode.substr(0,2)=="EM") {
        // thresholds and gains set randomly to reproduce EM1 values
        // chips and strips are randomized separately
        int mySeed = 123456789;
        HepRandom::setTheSeed(mySeed);
        for(tower=0;tower<NTOWERS;++tower) {
            for (layer=0;layer<numLayers;++layer) {
                for (view=0;view<NVIEWS;++view) {
                    for(chip=0;chip<nChips;++chip) {
                        // generate chip thresholds and gains
                        double chipGain = RandGauss::shoot(1.789, 0.3088);
                        double chipThresh = -0.8546 - 0.2142*chipGain + RandGauss::shoot(0.013, 0.167);
                        for (strip=0;strip<nStrips;++strip) {
                            int theStrip = chip*nStrips + strip;
                            double test = RandFlat::shoot(432.);
                            double devGain;
                            if (test>320.) {
                                devGain = RandGauss::shoot(0.188, 0.399);
                            } else {
                                devGain = RandGauss::shoot(-0.107, 0.262);
                            }
                            double devThresh = -0.541*devGain + RandGauss::shoot(0, 0.262);

                            m_ToTGain[tower][layer][view][theStrip] = 
                                chipGain + devGain;
                            m_ToTGain2[tower][layer][view][theStrip] = 
                                m_defaultGain2;
                            m_ToTThreshold[tower][layer][view][theStrip] = 
                                chipThresh + devThresh;
                            m_ToTQuality[tower][layer][view][theStrip] = 
                                m_defaultQuality;
                            m_ToTMuonFactor[tower][layer][view][theStrip] = 
                                m_defaultMuonFactor;
                        }
                    }
                }
            }
        }
    } else {
        log << MSG::ERROR << "Called with mode: """ << m_mode 
            << """, should be ""default"" or ""EM"" "
            << std::endl;
        sc = StatusCode::FAILURE;
    }

    return sc;
}

double TkrToTSvc::getCharge(double ToT, int tower, int layer, int view, int strip) const
{
    double gain      = m_ToTGain[tower][layer][view][strip];
    double gain2     = m_ToTGain2[tower][layer][view][strip];
    double threshold = m_ToTThreshold[tower][layer][view][strip];
    double charge;
    // constants are for: ToT = threshold + charge*(gain + charge*gain2))
    // here is the inverse:
    double term = (threshold-ToT/m_countsPerMicrosecond)/gain;
    double test = gain2/gain;
    if (fabs(test)>1.e-6) {
        // just the quadratic formula
        charge = 0.5*(-1.0 + sqrt(1.0 - 4.*test*term))/test;
    } else {
        // degenerate case
        charge = -term*(1.0 + test*term);
    }
    return charge;
}

double TkrToTSvc::getMipsFromToT(double ToT, int tower, int layer, int view, int strip) const
{
    double muonFactor = m_ToTMuonFactor[tower][layer][view][strip];
    return muonFactor/getFCPerMip()*getCharge(ToT, tower, layer, view, strip);
}

double TkrToTSvc::getMipsFromCharge(double charge, int tower, int layer, int view, int strip) const
{
    double muonFactor = m_ToTMuonFactor[tower][layer][view][strip];
    return muonFactor/getFCPerMip()*charge;
}

StatusCode TkrToTSvc::finalize() {
    
    MsgStream log(msgSvc(), name());
    return StatusCode::SUCCESS;
}
