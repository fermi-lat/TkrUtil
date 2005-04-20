// Implementation file for TkrFailureModeSvc which handles the failure mode testing
// for the Tkr.
// 
//
// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrFailureModeSvc.cxx,v 1.17 2005/03/01 00:57:46 lsrea Exp $
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
    m_existsList = false;
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

    // Open the message log
    MsgStream log( msgSvc(), name() );

    StatusCode  sc = StatusCode::SUCCESS;

    Service::initialize();

    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrGeometrySvc" 
            << endreq;
        return sc;
    }

    if(m_visitor==0) {
        m_visitor = new BadVisitorFM;
        m_visitor->setService(this);
        m_visitor->setService(m_tkrGeom);
    }

    // Bind all of the properties for this service
    if ( (sc = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }

    sc = doInit();

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
    StatusCode sc = StatusCode::SUCCESS;
    m_visitor->setLog(&log);
    log << MSG::INFO << "Updater called" << endreq;
    //StatusCode sc = doInit();
    if (!m_existsList) {
        if (pDead) pDead->traverse(m_visitor);
        if (pHot)  pHot->traverse(m_visitor);
        m_failureModes = 0;
        if(m_towerList.size()) m_failureModes |= 1<<TOWER_SHIFT;
        if(m_layerList.size()) m_failureModes |= 1<<LAYER_SHIFT;
    } else {
        log << MSG::INFO << "No update done -- layer and tower list is being used instead" << endreq;
    }
    /*
    std::cout << "layerList size " << m_layerList.size() << std::endl;
    int tower;
    for (tower=0;tower<1;++tower) {
        LayerMap::const_iterator iList = m_layerList.find(tower);
        if (iList==m_layerList.end()) continue;

        const std::vector<int> &layerList = iList->second;
        int ii;
        for (ii=0;ii<layerList.size();++ii) {
            std::cout << layerList[ii] <<std::endl;
        }
    }
    */

    return sc;

}

StatusCode TkrFailureModeSvc::finalize () 
{
    delete m_visitor;
    return StatusCode::SUCCESS;
}

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
        m_failureModes = m_failureModes | 1 << LAYER_SHIFT;
        m_existsList = true;

        /*
        std::cout << "debug output from processLayerList" << std::endl;
        int itower;
        for (itower=0;itower<1;++itower) {
            LayerMap::const_iterator iList = m_layerList.find(itower);
            if (iList==m_layerList.end()) continue;

            const std::vector<int> &layerList = iList->second;
            int vsize = layerList.size();
            std::vector<int>::const_iterator loc;
            for (loc=layerList.begin(); loc!=layerList.end();++loc) {
                int ii;
                for (ii=0;ii<vsize;++ii) {
                    std::cout << layerList[ii] <<std::endl;
                }
            }
        }
        */
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
        m_failureModes = m_failureModes | 1 << TOWER_SHIFT;
        m_existsList = true;
    }
}

bool TkrFailureModeSvc::towerFailed(int tower) const {
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

    int tower = idents::TowerId(col, row).id();
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

        int tower = idents::TowerId(col, row).id();

        int layer, view;
        m_tkrGeom->trayToLayer(tray, top, layer, view);
        int plane = m_tkrGeom->trayToPlane(tray, top);

        std::vector<int>& curList = m_pFailureMode->getLayers(tower);
        unsigned int i = 0;
        bool found = false;
        for(; i<curList.size();++i) {
            if (curList[i]==plane) {
                found = true;
                break;
            }
        }
        if (!found) {curList.push_back(plane);}       
    }

return CalibData::CONT;
}

