
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "src/TkrGeometrySvc.h"

#include "idents/TowerId.h"

#include <iostream>
#include <algorithm>

static const SvcFactory<TkrGeometrySvc> s_factory;
const ISvcFactory& TkrGeometrySvcFactory = s_factory;


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrGeometrySvc::TkrGeometrySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   
    return; 
}

StatusCode TkrGeometrySvc::initialize()
{
    // Purpose: load up constants from GlastDetSvc and do some calcs
    // Inputs:  none
    // Output:  TkrGeometrySvc statics initialized
    
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    double siWaferActiveSide;

    m_nviews = 2;

    if (service("GlastDetSvc", m_pDetSvc).isSuccess() &&
        
        m_pDetSvc->getNumericConstByName("xNum", &m_numX).isSuccess() &&
        m_pDetSvc->getNumericConstByName("xNum", &m_numY).isSuccess() &&    
        m_pDetSvc->getNumericConstByName("nWaferAcross", &m_nWaferAcross).isSuccess() &&   
        m_pDetSvc->getNumericConstByName("numTrays", &m_nlayers).isSuccess() &&   
        m_pDetSvc->getNumericConstByName("towerPitch", &m_towerPitch).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiThick", &m_siThickness).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiWaferSide", &m_siWaferSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName(
        "SiWaferActiveSide", &siWaferActiveSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName("stripPerWafer", &m_ladderNStrips).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ladderGap", &m_ladderGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("numSuperGlast", &m_nSuperGlast).isSuccess() &&
        m_pDetSvc->getNumericConstByName("numNoCnvTrays", &m_nNoConverter).isSuccess() )
    {
        sc = StatusCode::SUCCESS;
    } else {
        log << MSG::ERROR << "Failed to get geometry constants" << endreq;
        return StatusCode::FAILURE;
    }
    
    // finish up the constants
    m_nlayers--;
    m_nNoConverter--;
    m_trayWidth = m_nWaferAcross*m_siWaferSide +(m_nWaferAcross-1)*m_ladderGap;
    m_siDeadDistance = 0.5*(m_siWaferSide - siWaferActiveSide);
    m_siStripPitch = siWaferActiveSide/m_ladderNStrips;
    m_siResolution = m_siStripPitch/sqrt(12.);
    
    
    // fill up the m_volId arrays, used for providing volId prefixes
    
    for(int tower=0;tower<m_numX*m_numY;tower++) {
        idents::VolumeIdentifier vId;
        vId.append(0);
        idents::TowerId t(tower);
        vId.append(t.iy());
        vId.append(t.ix());
        vId.append(1);
        m_volId_tower[tower].init(0,0);
        m_volId_tower[tower].append(vId);
    }
    
    int bilayer;
    for(bilayer=0;bilayer<m_nlayers;bilayer++) {
        for (int view=0; view<2; view++) {
            int tray;
            int botTop;            
            layerToTray(bilayer, view, tray, botTop);

            idents::VolumeIdentifier vId;
            vId.append(tray);
            vId.append(view);
            vId.append(botTop);
            // seems that the old silicon plane no longer exists, 
            // only wafers now
            // add in ladder, wafer--this is all fragile
            vId.append(0); vId.append(0); 
            
            m_volId_layer[bilayer][view].init(0,0);
            m_volId_layer[bilayer][view].append(vId);
        }
    }   
    
    // the minimum "trayHeight" (actually tray pitch)
    sc = getMinTrayHeight(m_trayHeight);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to get minTrayHeight"<< endreq;
        return sc;
    }

    // fill the m_layerZ arrays
    sc = fillLayerZ();
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to fill layerZ"<< endreq;
        return sc;
    }
    return sc;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, 
                                            double stripId)
{
    // Purpose: return the global position
    // Method:  pass it on to the detector service
    // Inputs:  (tower, bilayer, view) and strip number (can be fractional)
    // Return:  global position

    idents::VolumeIdentifier volId;
    volId.append(m_volId_tower[tower]);
    volId.append(m_volId_layer[layer][view]);
    return m_pDetSvc->getStripPosition(volId, stripId);
}

void TkrGeometrySvc::trayToLayer(int tray, int botTop, int& layer, int& view)
{
    // Purpose: calculate layer and view from tray and botTop
    // Method:  pass it on to the detector service
    
    m_pDetSvc->trayToLayer(tray, botTop, layer, view);
}

void TkrGeometrySvc::layerToTray(int layer, int view, int& tray, int& botTop) 
{   
    // Purpose: calculate tray and botTop from layer and view
    // Method:  pass it on to the detector service
    
    m_pDetSvc->layerToTray(layer, view, tray, botTop);
}


// queryInterface

StatusCode  TkrGeometrySvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrGeometrySvc == riid) {
        *ppvIF = dynamic_cast<ITkrGeometrySvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service

const IID&  TkrGeometrySvc::type () const {
    return IID_ITkrGeometrySvc;
}


StatusCode TkrGeometrySvc::getMinTrayHeight(double &trayHeight) 
{
    // Purpose: fills the m_trayHeight
    // Method:  cycles through the layers and measure differences in z
    //          measures from top of a layer to top of next lowest layer
    // Inputs:  uses list of volIds of layers
    // Outputs: passes back the minimum tray height
    // Caveats:

    StatusCode sc = StatusCode::SUCCESS;
    HepTransform3D T1, T2;
    trayHeight = 10000.0;
    
    // Geometry service knows about trays, recon wants layers
    for (int tray = 2; tray<=m_nlayers; tray++) {

        idents::VolumeIdentifier volId1, volId2;
        int layer1, layer2;
        int view1, view2;
        
        volId1.append(m_volId_tower[0]);
        volId2.append(m_volId_tower[0]);

        int botTop = 0;
        trayToLayer(tray,   botTop, layer1, view1);    // bottom of each tray
        trayToLayer(tray-1, botTop, layer2, view2);
        
        volId1.append(m_volId_layer[layer1][view1]);
        volId2.append(m_volId_layer[layer2][view2]);

        std::string foo1 = volId1.name();
        std::string foo2 = volId1.name();

        sc = m_pDetSvc->getTransform3DByID(volId1, &T1);
        sc = sc && m_pDetSvc->getTransform3DByID(volId2, &T2);

        double z1 = (T1.getTranslation()).z();
        double z2 = (T2.getTranslation()).z();
        double trayPitch = z1 - z2;

        trayHeight = std::min(trayPitch, trayHeight);  
    }
    return sc;
}

StatusCode TkrGeometrySvc::fillLayerZ() 
{
    // Purpose: Fills the m_layerZ arrays with z positions
    // Method:  cycles through the layers and recovers z from the geometry
    // Inputs:  list of volIds of layers
    // Outputs: fills m_layerZ
    // Caveats:

    StatusCode sc = StatusCode::SUCCESS;
    HepTransform3D T;
    
    for (int bilayer=0; bilayer<m_nlayers; bilayer++) {
        for (int view=0; view<2; view++) {
            idents::VolumeIdentifier volId;
            volId.append(m_volId_tower[0]);
            volId.append(m_volId_layer[bilayer][view]);
            std::string foo = volId.name();
            sc = sc && m_pDetSvc->getTransform3DByID(volId,&T);
            m_layerZ[bilayer][view] = (T.getTranslation()).z();
           
        }
    }
    return sc;
}

double TkrGeometrySvc::getReconLayerZ(int layer, int view)
{
    // Purpose: returns the z for a given plane and view
    // Method:  accesses m_layerZ
    // Inputs   layer and view
    // Outputs: z position
    // Caveats:

    int digiLayer = reverseLayerNumber(layer);
    switch (view) {
    case 0:
        return m_layerZ[digiLayer][0];
        break;
    case 1:
        return m_layerZ[digiLayer][1];
        break;
    default:
        return 0.5*(m_layerZ[digiLayer][0] + m_layerZ[digiLayer][1]);
    }
}
    
double TkrGeometrySvc::getReconLayerZ(int layer) 
{
    // Purpose: returns the average z for a given layer
    // Method:  call getReconLayerZ with 2nd argument
    // Inputs   layer
    // Outputs: z position
    // Caveats:
    return getReconLayerZ(layer, 2);
}
