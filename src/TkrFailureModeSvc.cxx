// Implementation file for TkrFailureModeSvc which handles the failure mode testing
// for the Tkr.
// 
//
// $Header$
//
// Author: L. Rochester (after Richard Dubois)


#include "TkrUtil/TkrFailureModeSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include <algorithm>


// declare the service factories for the TkrFailureModeSvc
static SvcFactory<TkrFailureModeSvc> a_factory;
const ISvcFactory& TkrFailureModeSvcFactory = a_factory; 

TkrFailureModeSvc::TkrFailureModeSvc(const std::string& name,ISvcLocator* svc) : Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // declare the properties

    declareProperty("towerList", m_towerListProperty);
    declareProperty("layerList", m_layerListProperty);
}

StatusCode  TkrFailureModeSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrFailureModeSvc == riid) {
        *ppvIF = dynamic_cast<ITkrFailureModeSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

const IID&  TkrFailureModeSvc::type () const {
    return IID_ITkrFailureModeSvc;
}

StatusCode TkrFailureModeSvc::initialize () 
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
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
        
    m_failureModes = 0;
    processTowerList();
    processLayerList();
    
    return StatusCode::SUCCESS;
}

StatusCode TkrFailureModeSvc::finalize () {return StatusCode::SUCCESS;}

void TkrFailureModeSvc::processLayerList() {
    // Purpose and Method: process the jobOptions input lists of (tower,layer) pairs
    //                     
    
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theLayers = m_layerListProperty.value( );
    if (theLayers.empty()) return;

    m_failureModes = m_failureModes || 1 << LAYER;

    log << MSG::DEBUG << "Layers to kill " << endreq;

    std::vector<std::string>::const_iterator it;
    std::vector<std::string>::const_iterator itend = theLayers.end( );
    for (it = theLayers.begin(); it != itend; it++) {
        int len      = (*it).size();
        int delimPos = (*it).find_first_of('_');
        int tower    = atoi((*it).substr(0, delimPos).c_str());

        std::string remainder = (*it).substr(delimPos+1);
        int len1      = remainder.size();
        int delimPos1 = remainder.find_first_of("_");

        int layer     = atoi(remainder.substr(0, delimPos1).c_str());
        int view      = atoi(remainder.substr(delimPos1+1, len1-delimPos1-1).c_str());
        int plane     = 2*layer + view;
        
        log << MSG::DEBUG << "Tower " << tower << " Layer " << layer << " View " << view << endreq;
        std::vector<int>& curList = m_layerList[tower];
        curList.push_back(plane);                
    }
}

void TkrFailureModeSvc::processTowerList() {
     // Purpose and Method: process the jobOptions input lists of towers
    //                     
   
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theTowers = m_towerListProperty.value( );

    if (theTowers.empty()) return;

    log << MSG::DEBUG << "Towers to kill: " << endreq;

    m_failureModes = m_failureModes || 1 << TOWER;

    std::vector<std::string>::const_iterator it;
    std::vector<std::string>::const_iterator itend = theTowers.end( );
    for (it = theTowers.begin(); it != itend; it++) {
        int tower = atoi((*it).c_str());
        log << MSG::DEBUG << "Tower " << tower << endreq;
        m_towerList.push_back(tower);
    }
}

bool TkrFailureModeSvc::towerFailed(int tower) {
    bool found = false;
    // Search to see if this event id is among the list of ids we want to pause on
    int *loc = std::find(m_towerList.begin(), m_towerList.end(), tower);                
    return (loc != m_towerList.end());
}


bool TkrFailureModeSvc::isFailed(int tower, int layer, int view) {
    // Purpose and Method: check whether given id is in any of the identified lists
    //                     

    bool found = false;

    // check tower first

    if(found = towerFailed(tower)) return found;
    
    // are we looking any further?
    if (layer==-1) return found;

    return found = layerFailed(tower, layer, view);
}

bool TkrFailureModeSvc::layerFailed(int tower, int layer, int view) {
    // Purpose and Method: look for the given id in the tower list
    //                        
    if (m_layerList.empty()) return false;

    // just a number, no physical relation to anything!
    int plane = 2*layer + view; 

    
    std::vector<int> &layerList = m_layerList[tower];
    
    // Search to see if this (tower,layer) is among the list
    int *loc = std::find(layerList.begin(), layerList.end(), plane);                
    
    return (loc != layerList.end());
}


