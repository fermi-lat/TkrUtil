
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
        m_pDetSvc->getNumericConstByName("towerPitch", &m_towerPitch).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiThick", &m_siThickness).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiWaferSide", &m_siWaferSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName(
        "SiWaferActiveSide", &siWaferActiveSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName("stripPerWafer", &m_ladderNStrips).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ladderGap", &m_ladderGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap).isSuccess() 
        )
    {
        sc = StatusCode::SUCCESS;
    } else {
        log << MSG::ERROR << "Failed to get geometry constants" << endreq;
        return StatusCode::FAILURE;
    }

    // finish the constants
    m_trayWidth = m_nWaferAcross*m_siWaferSide +(m_nWaferAcross-1)*m_ladderGap;
    m_siDeadDistance = 0.5*(m_siWaferSide - siWaferActiveSide);
    m_siStripPitch = siWaferActiveSide/m_ladderNStrips;
    m_siResolution = m_siStripPitch/sqrt(12.);
    
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

    // fill up the m_volId tower arrays, used for providing volId prefixes
    // these are the Id's for the 15 towers, with just the first 3 fields

    int tower;
    for(tower=0;tower<m_numX*m_numY;++tower) {
        idents::VolumeIdentifier vId;
        vId.append(0);               // in Tower
        idents::TowerId t(tower);  
        vId.append(t.iy());          // yTower
        vId.append(t.ix());          // xTower
        m_volId_tower[tower].init(0,0);  
        m_volId_tower[tower].append(vId);
    }
 
    // fill the towerType array, and find the test tower
    m_testTower = -1;
    m_xLim[0] = 1000; m_xLim[1] = -1;
    m_yLim[0] = 1000; m_yLim[1] = -1;
    HepTransform3D T;
    StatusCode foundTower = StatusCode::FAILURE;
    for (tower=0;tower<m_numX*m_numY;++tower) {
        idents::TowerId t = idents::TowerId(tower);
        m_towerType[tower] = 0;
        idents::VolumeIdentifier volId;
        volId.init(0,0);
        volId.append(m_volId_tower[tower]);
        volId.append(1); // TKR

        int tray;
        int botTop;
        layerToTray(0, 0,tray, botTop);
        volId.append(tray);
        volId.append(0);
        volId.append(botTop);
        volId.append(0); volId.append(0); // ladder and wafer
        sc = m_pDetSvc->getTransform3DByID(volId,&T);
        if (sc.isFailure()) {
            // start by marking missing towers
            m_towerType[tower] = -1;
            continue;
        } else {
            foundTower = StatusCode::SUCCESS;
            // find the lowest and highest towers in each direction
            m_xLim[0] = std::min(m_xLim[0], t.ix());
            m_xLim[1] = std::max(m_xLim[1], t.ix());
            m_yLim[0] = std::min(m_yLim[0], t.iy());
            m_yLim[1] = std::max(m_yLim[1], t.iy());
        }
           
        // set test tower to first tower actually present
        if(m_testTower<0) m_testTower = tower;  
    }
    if( foundTower.isFailure()) {
        // for now, just fail; might be more clever later
        log << MSG::ERROR << "Failed to find any tower... check geometry!"<< endreq;
        return sc;
    }

    // use the list of existing towers to generate the tower type of each tower.
    // tower type is number of edges not touching another tower (0-4);
    int ix, iy;
    idents::TowerId tempTower;
    for (ix=0; ix<m_numX; ++ix) {
        for (iy=0; iy<m_numY; ++iy) {
            int nExposed = 0;
            idents::TowerId t = idents::TowerId(ix, iy);
            tower = t.id();
            if (m_towerType[tower]==-1) continue;
            // tower exists, check the 4 sides
            if (ix==0) ++nExposed;
            else {
                tempTower = idents::TowerId(ix-1,iy);
                if (m_towerType[tempTower.id()]==-1) ++nExposed;
            }
            if (ix==m_numX-1) ++nExposed; 
            else {
                tempTower = idents::TowerId(ix+1,iy);
                if (m_towerType[tempTower.id()]==-1) ++nExposed;
            }
            if(iy==0) ++nExposed;
            else {
                tempTower = idents::TowerId(ix,iy-1);
                if (m_towerType[tempTower.id()]==-1) ++nExposed;
            }
            if (iy==m_numY-1) ++nExposed;
            else {
                tempTower = idents::TowerId(ix,iy+1);
                if (m_towerType[tempTower.id()]==-1) ++nExposed;
            }
            m_towerType[tower] = nExposed;
        }
    }

    sc = fillPropagatorInfo(); 
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to fill rad len arrays"<< endreq;
        return sc;
    }
   
    int bilayer;
    int view;
    for(bilayer=0;bilayer<numLayers();++bilayer) {
        for (view=0; view<2; ++view) {
            int tray;
            int botTop;            
            layerToTray(bilayer, view, tray, botTop);

            idents::VolumeIdentifier vId;
            vId.append(tray);
            vId.append(view);
            vId.append(botTop);
            vId.append(0); vId.append(0); // ladder and wafer
            
            m_volId_layer[bilayer][view].init(0,0);
            m_volId_layer[bilayer][view].append(vId);
        }
    }  

    // fill the m_layerZ arrays
    sc = fillLayerZ();
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to fill layerZ"<< endreq;
        return sc;
    }
   
    // the minimum "trayHeight" (actually tray pitch)
    // uses layerZ info, so must follow call to fillLayerZ()
    sc = getMinTrayHeight(m_trayHeight);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to get minTrayHeight"<< endreq;
        return sc;
    }

	// Get cal info necessary for tracker recon

	sc = getCalInfo();
    if (sc.isFailure()) return sc;



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
    if( service( "TkrFailureModeSvc", m_tkrFail, true).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrFailureModeSvc" << endreq;
        log << "Will assume it is not required"    << endreq;
    }

    // Get the alignment service 
    m_tkrAlign = 0;
    if( service( "TkrAlignmentSvc", m_tkrAlign, true).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrAlignmentSvc" << endreq;
        log << "Will assume it is not required"    << endreq;
    }

    m_badStrips = 0;
    if( service( "TkrBadStripsSvc", m_badStrips, true).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrBadStripsSvc" << endreq;
        log << "Will assume it is not required"    << endreq;
    }

    m_tkrSplits = 0;
    if( service( "TkrSplitsSvc", m_tkrSplits, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't set up TkrSplitsSvc" << endreq;
        return StatusCode::FAILURE;
    }

    m_tkrToT    = 0;
    if( service( "TkrToTSvc", m_tkrToT, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't set up TkrToTSvc" << endreq;
        return StatusCode::FAILURE;
    }

    log << MSG::INFO << "TkrGeometrySvc successfully initialized" << endreq;
    return StatusCode::SUCCESS;

}

StatusCode TkrGeometrySvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "TkrGeometrySvc finalize called" << endreq;
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
    volId.append(1); // TKR
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


StatusCode TkrGeometrySvc::getMinTrayHeight(double& trayHeight) 
{
    // Purpose: fills the m_trayHeight
    // Method:  cycle through zlayer arrays
    // Inputs:  
    // Outputs: passes back the minimum tray height
    // Caveats:

    StatusCode sc = StatusCode::SUCCESS;
    //HepTransform3D T1, T2;
    trayHeight = 10000.0;
    
    // Geometry service knows about trays, recon wants layers

    for (int layer = 1; layer<numLayers(); ++layer) {
        double trayPitch = getReconLayerZ(layer-1) - getReconLayerZ(layer);
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
    
    for (int bilayer=0; bilayer<numLayers(); ++bilayer) {
        for (int view=0; view<2; ++view) {
            idents::VolumeIdentifier volId;
            volId.append(m_volId_tower[m_testTower]);
            volId.append(1); // TKR
            volId.append(m_volId_layer[bilayer][view]);
            if ((sc=m_pDetSvc->getTransform3DByID(volId,&T)).isFailure()) break;
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
    
    double stayClear;

    if(m_pDetSvc->getNumericConstByName("TKRVertStayClear", &stayClear).isFailure()) {
        log << MSG::DEBUG ;
        if(log.isActive()) log << "Couldn't find TKRVertStayClear";
        log << endreq;
        return StatusCode::FAILURE;
    }

    // get the z coordinate of the top silicon in the bottom tray 
    //         (should be pretty safe!)
    // don't make any assumptions about the view in the bottom tray
    idents::VolumeIdentifier bottom;
    bottom = m_volId_tower[m_testTower]; // testTower
    bottom.append(1); // TKR
    bottom.append(0);                // tray 0
    idents::VolumeIdentifier idBot;
    HepTransform3D botTransform;
    for (int view = 0; view<2; ++view) {
        idBot = bottom;
        idBot.append(view);          // try both views
        idBot.append(1);             // top silicon (*most* bottom trays have one!)
        idBot.append(0); idBot.append(0);  // wafer
        //std::cout << "view " << view << " idBot " << idBot.name() << std::endl;
        if(sc = m_pDetSvc->getTransform3DByID(idBot, &botTransform).isSuccess()) {
            break;
        }
    }
    if (sc.isFailure()) {
        log << MSG::DEBUG ;
        if(log.isActive()) log << "Couldn't find bottom tray silicon";
        log << endreq;
        return sc; }

    double zBot;
    zBot = (botTransform.getTranslation()).z();

    // this is always above the tracker
    double propTop = stayClear+zBot;

    // need to do the swim through the test tower
    idents::TowerId tTower = idents::TowerId(m_testTower);
    int xTower = tTower.ix();
    int yTower = tTower.iy();
    double xStart = (xTower - 0.5*(numXTowers()-1))*towerPitch() + 40.;
    double yStart = (yTower - 0.5*(numYTowers()-1))*towerPitch() + 40.;
    Point startPoint = Point(xStart, yStart, propTop);

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
    // if the top were at the top of the stayclear, stayClear would be enuf
    // 100 mm gets you from the middle of the bottom closeout to below the tracker
    double propRange = stayClear+100.;
    double propBot = propTop - propRange;
    track->step(propRange);
    log << MSG::INFO  << "Propagator goes from "<< propTop << " to " << propBot << std::endl;
    
    int numSteps = track->getNumberSteps();
    int istep = 0;
    idents::VolumeIdentifier id;
    double radlen;
    idents::VolumeIdentifier prefix = m_pDetSvc->getIDPrefix();

    bool startCount = false;
    for (; istep<numSteps; ++istep) {
        Point stepPoint = track->getStepPosition(istep);
        id = track->getStepVolumeId(istep);
        id.prepend(prefix);
        std::string idName = id.name();
        //double arclen = track->getStepArcLen(istep);
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
            // the "bottom" silicon is one tray up
            int layer = (item==2 || item==6 || item==0) ? tray-1 : tray;
            if (item==2) { // converter
                m_radLenConv[layer] = radlen;
            } else {
                m_radLenRest[layer] += radlen;
            }
        }
    }

    log << MSG::DEBUG;

    for (i=0; i<NLAYERS; ++i) {
        convType type = getDigiLayerType(i);
        if (type==ABSENT) continue;
        ind = (int) type;
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

    // choose layer types based on thickness of converter:
    //   no r.l. anywhere in layer -> Absent
    //   under 1% -> NoConverter
    //   over 10% -> SuperGlast
    //   else        Standard

    double radlen = m_radLenConv[digiLayer];
    convType type;
    if (radlen<0.01) {type = (m_radLenRest[digiLayer]> 0) ? NOCONV : ABSENT;}
    else {type = (radlen>0.10) ? SUPER : STANDARD; }

    return type;
}

StatusCode TkrGeometrySvc::getCalInfo()
{
    MsgStream log(msgSvc(), name());

	double csiLength, csiWidth, csiHeight;
    double cellHorPitch, cellVertPitch;
    int nCsiPerLayer, CALnLayer;

    if (m_pDetSvc->getNumericConstByName("CsILength",  &csiLength).isFailure() ||
        m_pDetSvc->getNumericConstByName("CsIHeight",  &csiHeight).isFailure() ||
        m_pDetSvc->getNumericConstByName("CsIWidth",   &csiWidth ).isFailure() ||
        m_pDetSvc->getNumericConstByName("CALnLayer",  &CALnLayer).isFailure() ||
        m_pDetSvc->getNumericConstByName("nCsIPerLayer",  &nCsiPerLayer).isFailure() ||
        m_pDetSvc->getNumericConstByName("cellHorPitch",  &cellHorPitch).isFailure() ||
        m_pDetSvc->getNumericConstByName("cellVertPitch", &cellVertPitch).isFailure())
    {return StatusCode::FAILURE;}

    //get the top and bottom of the CAL crystals

    int count;

    // top layer of the Cal
	// Which id is legal depends on whether the cal geometry is plain or segvols.
	// for the moment, I'll just check a few until I get a good one, and if I don't then
	// I bail... this is a bit of a kludge... seems to work though

    idents::VolumeIdentifier topLayerId;
    topLayerId.append(m_volId_tower[m_testTower]);
    topLayerId.append(0);  // CAL
    topLayerId.append(0);  // layer
    topLayerId.append(0);  // x view
 
    StatusCode sc;
    HepTransform3D transfTop;
    for (count=0;count<3;++count) {
        topLayerId.append(0);
        if((sc = m_pDetSvc->getTransform3DByID(topLayerId,&transfTop)).isSuccess()) break;
    }
	if(sc.isFailure()) {
		log << MSG::ERROR << "Couldn't get Id for layer 0 of CAL" << endreq;
		return sc;
	}
    Vector vecTop = transfTop.getTranslation();
    m_calZTop = vecTop.z()+ 0.5*csiHeight;
    
    m_calZBot = m_calZTop - (CALnLayer-1)*cellVertPitch - csiHeight;

    // get the maximum horizontal dimension of the crystals in a layer
    double calWidth1 = (nCsiPerLayer-1)*cellHorPitch + csiWidth;
    double ModWidth  = std::max(calWidth1, csiLength);
    m_calXWidth = (m_numX-1)*m_towerPitch + ModWidth;
    m_calYWidth = (m_numY-1)*m_towerPitch + ModWidth;

    return StatusCode::SUCCESS;
}

double TkrGeometrySvc::getLATLimit(int view, limitType type) const
{
    int towerNum = getLimitingTower(view, type);
    int nTowers  = (view==0 ? m_numX     : m_numY );
    // lower edge
    double limit = (-0.5*nTowers + towerNum)*m_towerPitch;
    // upper edge
    if (type==HIGH) { limit += m_towerPitch; }
    
    return limit;
}

bool TkrGeometrySvc::isInActiveLAT(Point pos) const
{
    double x = pos.x();
    double y = pos.y();
    return (x>getLATLimit(0,LOW) && x<getLATLimit(0,HIGH) 
        && y>getLATLimit(1,LOW) && y<getLATLimit(1,HIGH) ? true : false );
}