
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/TkrGeometrySvc.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"

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
        //m_pDetSvc->getNumericConstByName("numTrays", &m_nlayers).isSuccess() &&   
        m_pDetSvc->getNumericConstByName("towerPitch", &m_towerPitch).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiThick", &m_siThickness).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiWaferSide", &m_siWaferSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName(
        "SiWaferActiveSide", &siWaferActiveSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName("stripPerWafer", &m_ladderNStrips).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ladderGap", &m_ladderGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap).isSuccess() //&&
        //m_pDetSvc->getNumericConstByName("numSuperGlast", &m_nSuperGlast).isSuccess() &&
        //m_pDetSvc->getNumericConstByName("numNoCnvTrays", &m_nNoConverter).isSuccess() 
        )
    {
        sc = StatusCode::SUCCESS;
    } else {
        log << MSG::ERROR << "Failed to get geometry constants" << endreq;
        return StatusCode::FAILURE;
    }
    
    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved G4PropagationTool" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve G4PropagationTool" << endreq;
        }
        
    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    }   
    
    sc = fillPropagatorInfo();
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to fill rad len arrays"<< endreq;
        return sc;
    }

    //m_nlayers      = getNumType(ALL);
    //m_nSuperGlast  = getNumType(SUPER);
    //m_nNoConverter = getNumType(NOCONV);
    //m_nRegular     = getNumType(REGULAR);
    
    // finish up the constants
    //m_nlayers--;
    //m_nNoConverter--;
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
    for(bilayer=0;bilayer<numLayers();bilayer++) {
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

    // Get propagator from the service
    
    IPropagatorSvc* propagatorSvc = 0;
    sc = service("GlastPropagatorSvc", propagatorSvc, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "GlastPropagatorSvc is required for this algorithm." 
            << endreq;
        return sc;
    }
    m_KalParticle = propagatorSvc->getPropagator();
   
       // Get the failure mode service 
    m_tkrFail = 0;
    if( service( "TkrFailureModeSvc", m_tkrFail, false).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrFailureModeSvc" << endreq;
        log << MSG::INFO << "Will assume it is not required"    << endreq;
    }

    // Get the alignment service 
    m_tkrAlign = 0;
    if( service( "TkrAlignmentSvc", m_tkrAlign, false).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrAlignmentSvc" << endreq;
        log << MSG::INFO << "Will assume it is not required"    << endreq;
    }

    m_badStrips = 0;
    if( service( "TkrBadStripsSvc", m_badStrips, false).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrBadStripsSvc" << endreq;
        log << MSG::INFO << "Will assume it is not required"    << endreq;
    }
    return StatusCode::SUCCESS;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, 
                                            double stripId) const
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


void TkrGeometrySvc::trayToLayer(int tray, int botTop, 
                                       int& layer, int& view) const
{
    // Purpose: calculate layer and view from tray and botTop
    // Method: use knowledge of the structure of the Tracker
    
    int plane = 2*tray + botTop - 1;
    layer = (plane+10)/2 - 5; // homemade "floor"
    view = ((layer%2==0) ? botTop : (1 - botTop));
    return;
}

void TkrGeometrySvc::layerToTray(int layer, int view, 
                                       int& tray, int& botTop) const
{   
    // Purpose: calculate tray and botTop from layer and view.
    // Method:  use knowledge of the structure of the Tracker
    
    int plane = (2*layer) + (((layer % 2) == 0) ? (1 - view) : (view));
    tray = (plane+1)/2;
    botTop = (1 - (plane % 2));
}


void TkrGeometrySvc::planeToLayer(int plane, 
                                       int& layer, int& view) const
{
    // Purpose: calculate tray and botTop from plane
    // Method:  use knowledge of the structure of the Tracker
        layer = plane/2;
        int element = (plane+3)%4;
        view = element/2;
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
    for (int tray = 2; tray<=numLayers(); tray++) {

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
    
    for (int bilayer=0; bilayer<numLayers(); bilayer++) {
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

StatusCode TkrGeometrySvc::fillPropagatorInfo()
{
    // Purpose: Fills the m_radLen arrays with radiation lengths
    // Method:  runs G4PropagationTool and accesses step info
    // Inputs:  
    // Outputs: filled arrays; stored in digi order 
    // Caveats:

    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());


    Point startPoint = Point(40., 40., 700.);
    Vector startDir  = Vector(0., 0., -1.);

    int i, ind;

    for(i=0; i<NLAYERS; ++i) {
        m_radLenConv[i] = 0.;
        m_radLenRest[i] = 0.;
    }
    for (ind=0; ind<NTYPES; ++ind) {
        m_numLayers[ind]     = 0;
        m_aveRadLenConv[ind] = 0.;
        m_aveRadLenRest[ind] = 0.;
    }

    IPropagator* track = m_G4PropTool;
    track->setStepStart(startPoint, startDir);
    track->step(699.);
    
    int numSteps = track->getNumberSteps();
    int istep = 0;
    idents::VolumeIdentifier id;
    double radlen;

    bool startCount = false;
    for (; istep<numSteps; ++istep) {
        Point stepPoint = track->getStepPosition(istep);
        id = track->getStepVolumeId(istep);
        std::string idName = id.name();
        double arclen = track->getStepArcLen(istep);
        radlen = track->getStepRadLength(istep);
        /*
        std::cout << "step " << istep << " rl " << radlen 
            << " arclen " << arclen 
            << " volid " << id.name()
            << " pos " << stepPoint <<   std::endl;
        */
        //check if in a layer
        if(id[0]==0 && id[3]==1 && id.size()>6) {
            int tray = id[4];
            int item = id[6];
            if (item==2) {
                startCount = true;
            }
            if (!startCount) continue;
            // only the "top" silicon in each layer is in the same tray as layer
            // the others are one tray up
            int layer = (item==2 || item==6 || item==0) ? tray-1 : tray;
            if (item==2) { // converter
                m_radLenConv[layer] = radlen;
            } else {
                m_radLenRest[layer] += radlen;
            }
        }
    }
    
    // choose layer types based on thickness of converter:
    //   under 1% -> NoConverter
    //   over 10% -> SuperGlast
    //   else        Standard

    log << MSG::DEBUG;

    for (i=0; i<NLAYERS; ++i) {
        ind = (int) getDigiLayerType(i);
        radlen = m_radLenConv[i];
        m_numLayers[ind]++;
        m_aveRadLenConv[ind] += radlen;
        m_aveRadLenRest[ind] += m_radLenRest[i];
        m_numLayers[(int)ALL]++;
        m_aveRadLenConv[(int)ALL] += radlen;
        m_aveRadLenRest[(int)ALL] += m_radLenRest[i];

        if(log.isActive()) {
        log << " digiLayer " << i << " Conv " << m_radLenConv[i] << " Rest " 
            << m_radLenRest[i] << endreq;
        }
    }

    log << endreq;

// reverseLayerNumber() starts working here!


    log << MSG::INFO << endreq;
    for (ind=0;ind<NTYPES;++ind) {
        m_aveRadLenConv[ind] = m_numLayers[ind] ? m_aveRadLenConv[ind]/m_numLayers[ind] : 0;
        m_aveRadLenRest[ind] = m_numLayers[ind] ? m_aveRadLenRest[ind]/m_numLayers[ind] : 0;
        log << MSG::INFO  << "Type " << ind << " numLayers " << m_numLayers[ind] 
            << " aveConv " << m_aveRadLenConv[ind] << " aveRest " << m_aveRadLenRest[ind]
            << endreq;
    }  
    log<< MSG::INFO << endreq;
    
    return sc;
}


double TkrGeometrySvc::getReconLayerZ(int layer, int view) const
{
    // Purpose: returns the z for a given plane and view
    // Method:  accesses m_layerZ
    // Inputs   layer and optional view (if absent, return average)
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
   
convType TkrGeometrySvc::getReconLayerType(int layer) const
{
    // Purpose: returns the converter type for a layer
    // Method:  reverses layer, and calls getDigiLayerType()
    // Inputs   layer
    // Outputs: layer type
    // Caveats:
    
    return getDigiLayerType(reverseLayerNumber(layer));
}

convType TkrGeometrySvc::getDigiLayerType(int digiLayer) const
{
    // Purpose: returns the converter type for a layer
    // Method:  tests radlen for this layer
    // Inputs   layer
    // Outputs: layer type
    // Caveats:

    double radlen = m_radLenConv[digiLayer];
    convType type;
    if (radlen<0.01)      { type = NOCONV;}
    else if (radlen>0.10) { type = SUPER;}
    else                  { type = STANDARD;}

    return type;
}
