// Implementation file for TkrFailureModeSvc which handles the failure mode testing
// for the Tkr.
// 
//
// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrFailureModeSvc.cxx,v 1.7 2003/02/07 23:33:08 lsrea Exp $
//
// Author: L. Rochester (after Richard Dubois)


#include "src/TkrFailureModeSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "idents/TowerId.h"

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
    m_visitor = 0;
}

StatusCode  TkrFailureModeSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrFailureModeSvc == riid) {
        *ppvIF = dynamic_cast<ITkrFailureModeSvc*> (this);
        return StatusCode::SUCCESS;
    } else if (IID_ITkrFailureModeSvcCalib == riid) {
        *ppvIF = dynamic_cast<ITkrFailureModeSvcCalib*> (this);
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
    std::cout << "pointer to log file " << &log << std::endl;

    if(m_visitor==0) {
        m_visitor = new BadVisitorFM;
        m_visitor->setService(this);
    }
      
    // Call super-class
    Service::initialize ();
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
        
    status = doInit();
        
    return StatusCode::SUCCESS;
}

StatusCode TkrFailureModeSvc::doInit()
{
    m_failureModes = 0;
    processTowerList();
    processLayerList();

    return StatusCode::SUCCESS;
}

StatusCode TkrFailureModeSvc::update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot)
{
    MsgStream log(msgSvc(), name());
    m_visitor->setLog(&log);
    log << MSG::INFO << "updater called" << endreq;
    StatusCode sc = doInit();
    pDead->traverse(m_visitor);
    pHot->traverse(m_visitor);
    
    
	// taken from bad strips code
    /*int i, j;
    std::cout << "class TkrBadStripsSvc bad strip lists: " << std::endl;
    for(i=0;i<NELEMENTS;i++) {
        const stripCol* v = getBadStrips(i);
        int size = v->size();
        int tower = i/(NLAYERS*NVIEWS);
        int layer = i%(NLAYERS*NVIEWS)/2;
        int view  = i%2;
        if (size) {
            std::cout << " index " << i <<" tower " << tower << " layer "
                << layer << " view " << view <<" size " << size << std::endl << " strips " ;
            for (j=0;j<size;j++) {
                int strip = (*v)[j].getStripNumber();
                std::cout << strip<< " " ;
            }
            std::cout << std::endl;
        }
    }
	*/
    
    //std::cout << fillStream(std::cout) << "*" << std::endl << "*" << "The end!!" << std::endl;
    
    return sc;
    
}


StatusCode TkrFailureModeSvc::finalize () {return StatusCode::SUCCESS;}

void TkrFailureModeSvc::processLayerList() {
    // Purpose and Method: process the jobOptions input lists of (tower,layer) pairs
    //                     
    
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theLayers = m_layerListProperty.value( );
    if (theLayers.empty()) return;

    log << MSG::DEBUG;
    if (log.isActive() ) {
        log << "Layers to kill ";
    }
    log << endreq;

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
        
        log << MSG::DEBUG;
        if (log.isActive() ) {
            log << "Tower " << tower << " Layer " << layer << " View " << view;
        }
        log << endreq;
        std::vector<int>& curList = m_layerList[tower];
        curList.push_back(plane);
        m_failureModes = m_failureModes || 1 << LAYER_SHIFT;
    }
}

void TkrFailureModeSvc::processTowerList() {
     // Purpose and Method: process the jobOptions input lists of towers
    //                     
   
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theTowers = m_towerListProperty.value( );

    if (theTowers.empty()) return;

    log << MSG::DEBUG;
    if (log.isActive () ) {
        log << "Towers to kill: ";
    }
    log << endreq;

    std::vector<std::string>::const_iterator it;
    std::vector<std::string>::const_iterator itend = theTowers.end( );
    for (it = theTowers.begin(); it != itend; it++) {
        int tower = atoi((*it).c_str());
        log << MSG::DEBUG;
        if (log.isActive() ) {
            log << "Tower " << tower;
        }
        log << endreq;
        m_towerList.push_back(tower);
        m_failureModes = m_failureModes || 1 << TOWER_SHIFT;
    }
}

bool TkrFailureModeSvc::towerFailed(int tower) const {
    bool found = false;
    // Search to see if this event id is among the list of ids we want to pause on
    std::vector<int>::const_iterator loc = std::find(m_towerList.begin(), m_towerList.end(), tower);                
    return (loc != m_towerList.end());
}


bool TkrFailureModeSvc::isFailed(int tower, int layer, int view) const {
    // Purpose and Method: check whether given id is in any of the identified lists
    //                     

    bool found = false;

    // check tower first

    if(found = towerFailed(tower)) return found;
    
    // are we looking any further?
    if (layer==-1) return found;

    return found = layerFailed(tower, layer, view);
}

bool TkrFailureModeSvc::layerFailed(int tower, int layer, int view)  const {
    // Purpose and Method: look for the given id in the tower list
    //                        
    if (m_layerList.empty()) return false;

    // just a number, no physical relation to anything!
    int plane = 2*layer + view; 

    // because layerFailed is const, we have to go thru extra shenanigans to access the 
    // vector of failed planes.  (because map[] adds an element if it isn't there!
   
    LayerMap::const_iterator iList = m_layerList.find(tower);
   
    if (iList==m_layerList.end()) { return false; }

    const std::vector<int> &layerList = iList->second;
    std::vector<int>::const_iterator loc = std::find(layerList.begin(), layerList.end(), plane);  
    
    return (loc != layerList.end());
}

std::vector<int>& TkrFailureModeSvc::getLayers(int tower) 
{
    return m_layerList[tower];
}

std::vector<int>& TkrFailureModeSvc::getTowers()
{
    return m_towerList;
}

CalibData::eVisitorRet BadVisitorFM::badTower(unsigned int row, unsigned int col,
                                            int badness)
{
    *m_log << MSG::DEBUG;
    if ((*m_log).isActive()) {
        *m_log << "BadVisitor::badTower called with args " << endreq
        << "row = " << row << " col = "<< col 
        << " badness = " << badness;
    }
    *m_log << endreq;

    int tower = idents::TowerId(row, col).id();
    std::vector<int>& towerList = m_pFailureMode->getTowers();
    std::vector<int>::iterator p = std::find(towerList.begin(), towerList.end(), tower);
    if (p==towerList.end()) { towerList.push_back(tower); }

    return CalibData::CONT;  
}

CalibData::eVisitorRet BadVisitorFM::badPlane(unsigned int row, 
                                            unsigned int col, 
                                            unsigned int tray, bool top,
                                            int badness, bool allBad,
                                            const CalibData::StripCol& strips)
{
    *m_log << MSG::DEBUG;
    if ((*m_log).isActive()) {
        *m_log << "BadVisitor::badPlane called with args" << endreq
        << "row = " << row << ", col = " << col << ", tray = "
        << tray << endreq
        << "top = " << top << ", badness = " 
        << badness << " allBad = " << allBad << endreq
        << "Strip collection contains " << strips.size()
        << " strips. ";
    }
    *m_log << endreq;
    
    if (allBad) {
     
        int tower = idents::TowerId(row, col).id();
        int layer = top ? tray : tray-1;
        int view  = layer%2 ? 1-top : top;
        int plane = 2*layer + view;
              
        std::vector<int> curList = m_pFailureMode->getLayers(tower);
        std::vector<int>::iterator p = std::find(curList.begin(), curList.end(), plane);
        if (p==curList.end()) {curList.push_back(plane);}       
    }
    
    return CalibData::CONT;
}

