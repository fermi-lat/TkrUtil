/*
@file TkrToTSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrToTSvc.cxx,v 1.2 2004/03/12 05:49:22 lsrea Exp $

*/

#include "src/TkrToTSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "idents/TowerId.h"

#include <algorithm>
#include <fstream>
#include "facilities/Util.h"
#include "xml/IFile.h"

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
    // The default gain and threshold are chosen to reproduce the 
    // previous behavior of the ToT
    declareProperty("defaultGain",      m_defaultGain      = 2.50267833 );
    declareProperty("defaultThreshold", m_defaultThreshold = -2.92);
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
    
    /*
    m_geoSvc = 0;
    if( service( "TkrGeometrySvc", m_geoSvc, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't retrieve TkrGeometrySvc" << endreq;
        return StatusCode::FAILURE;
    }
    */
    
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

    // test of getting TkrGeometrySvc from inside TkrToTSvc... It works!
    //int stripsPerLadder  = m_geoSvc->ladderNStrips();

    // can be removed when geometry is iterfaced here
    const int NSTRIPS = 64;
    const int NCHIPS  = 24;

    int tower, layer, view, strip, chip;
    for(tower=0;tower<NTOWERS;++tower) {
        for (layer=0;layer<NLAYERS;++layer) {
            for (view=0;view<NVIEWS;++view) {
                for(chip=0;chip<NCHIPS;++chip) {
                    for (strip=0;strip<NSTRIPS;++strip) {
                        int theStrip = chip*NSTRIPS + strip;
                        m_ToTGain[tower][layer][view][theStrip] = 
                            m_defaultGain;
                        m_ToTThreshold[tower][layer][view][theStrip] = 
                            m_defaultThreshold;
                    }
                }
            }
        }
    }
    return sc;
}

StatusCode TkrToTSvc::finalize() {
    
    MsgStream log(msgSvc(), name());  
    return StatusCode::SUCCESS;
}
