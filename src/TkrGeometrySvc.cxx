
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/TkrGeometrySvc.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"

#include "idents/TowerId.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

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

    sc = getConsts();

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
    initializeArrays();

    // getTestTower checks to see if there is a "bottom Tray"
    // this is crucial for numbering of the planes, 
    //    so must be done before anything else!!!
    sc = getTestTower();
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to find test tower"<< endreq;
        return sc;
    }

    sc = getVolumeInfo();
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to find correct volumes for tracker"<< endreq;
        return sc;
    }

    // number of layers is now known...
    // but will get be calculated another way later... the two *should* agree!
    makeTowerIds();

    //sc = getTestTower();
    //if( sc.isFailure()) {
    //    log << MSG::ERROR << "Failed to find test tower"<< endreq;
    //    return sc;
    //}

    sc  = getTowerLimits();

    if( (sc = getTowerLimits()).isFailure()) {
        // for now, just fail; might be more clever later
        log << MSG::ERROR << "Failed to find any tower... check geometry!"<< endreq;
        return sc;
    }

    getTowerType();

    sc = fillPropagatorInfo(); 
    if( sc.isFailure()) {
        log << MSG::ERROR << "Failed to fill rad len arrays"<< endreq;
        return sc;
    }

    makeLayerIds();

    // the minimum "trayHeight" (actually tray pitch)
    // uses planeZ info, so must follow call to fillPlaneZ()
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
    // Get the bad strips service
    m_badStrips = 0;
    if( service( "TkrBadStripsSvc", m_badStrips, true).isFailure() ) {
        log << MSG::INFO << "Couldn't set up TkrBadStripsSvc" << endreq;
        log << "Will assume it is not required"    << endreq;
    }
    // get the splits service
    m_tkrSplits = 0;
    if( service( "TkrSplitsSvc", m_tkrSplits, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't set up TkrSplitsSvc" << endreq;
        return StatusCode::FAILURE;
    }
    // get the ToT service
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
    volId.append(m_volId_layer[layer][view]);
    return m_pDetSvc->getStripPosition(volId, stripId);
}


void TkrGeometrySvc::trayToLayer(int tray, int botTop, 
                                 int& layer, int& view) const
{
    // Purpose: calculate layer and view from tray and botTop
    // Method: use knowledge of the structure of the Tracker

    int plane = trayToPlane(tray, botTop);
    layer = m_planeToLayer[plane];
    view  = m_planeToView[plane];
    return;
}

void TkrGeometrySvc::layerToTray(int layer, int view, 
                                 int& tray, int& botTop) const
{   
    // Purpose: calculate tray and botTop from layer and view.
    // Method:  use knowledge of the structure of the Tracker

    int plane = m_layerToPlane[layer][view];
    if (plane==-1) {
        tray = -1;
        botTop = -1;
        return;
    }
    tray = planeToTray(plane);
    botTop = planeToBotTop(plane);
    return;
}


void TkrGeometrySvc::planeToLayer(int plane, 
                                  int& layer, int& view) const
{
    // Purpose: calculate tray and botTop from plane
    // Method:  use knowledge of the structure of the Tracker
    layer = m_planeToLayer[plane];
    view  = m_planeToView[plane];
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
    trayHeight = 10000.0;

    // Geometry service knows about trays, recon wants layers

    for (int layer = 1; layer<numLayers(); ++layer) {
        double trayPitch = getReconLayerZ(layer-1) - getReconLayerZ(layer);
        trayHeight = std::min(trayPitch, trayHeight);  
    }
    return sc;
}

StatusCode TkrGeometrySvc::getConsts()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    m_nviews = 2;

    double siWaferActiveSide;

    if (service("GlastDetSvc", m_pDetSvc).isSuccess() &&
        m_pDetSvc->getNumericConstByName("xNum", &m_numX).isSuccess() &&
        m_pDetSvc->getNumericConstByName("xNum", &m_numY).isSuccess() &&    
        m_pDetSvc->getNumericConstByName("nWaferAcross", &m_nWaferAcross).isSuccess() &&   
        m_pDetSvc->getNumericConstByName("towerPitch",   &m_towerPitch).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiThick",      &m_siThickness).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiWaferSide",  &m_siWaferSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName("SiWaferActiveSide", 
                                      &siWaferActiveSide).isSuccess() &&
        m_pDetSvc->getNumericConstByName("stripPerWafer", &m_ladderNStrips).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ladderGap",    &m_ladderGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("ssdGap",       &m_ladderInnerGap).isSuccess() &&
        m_pDetSvc->getNumericConstByName("nFeChips",     &m_chipsPerLadder).isSuccess()
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

    return sc;
}

void TkrGeometrySvc::initializeArrays()
{
    MsgStream log(msgSvc(), name());

    int layer, view, plane;
    for (plane=0; plane<NPLANES; ++plane) {
        m_planeToView[plane] = -1;
        m_planeToLayer[plane] = -1;
        m_isTopPlaneInLayer[plane] = false;
    }

    for (layer=0;layer<NLAYERS;++layer) {
        for (view=0;view<NVIEWS;++view) {
            m_layerToPlane[layer][view] = -1; 
        }
    }

    int i, ind;

    for(i=0; i<NLAYERS; ++i) {
        m_radLenConv[i] = 0.;
        m_radLenRest[i] = 0.;
        m_convZ[i]      = 0.;
    }
    for (ind=0; ind<NTYPES; ++ind) {
        m_numLayers[ind]     = 0;
        m_aveRadLenConv[ind] = 0.;
        m_aveRadLenRest[ind] = 0.;
    }
    for (ind=0; ind<m_numX*m_numY; ++ind) {
        m_towerType[ind] = -1;
    }

    // and for now...
    m_bottomTrayNumber = 0;
    m_topTrayNumber    = NLAYERS;
}

StatusCode TkrGeometrySvc::getTowerLimits()
{

    StatusCode sc;

    HepTransform3D T;
    int tower;
    // fill the towerType array
    m_xLim[0] = 1000; m_xLim[1] = -1;
    m_yLim[0] = 1000; m_yLim[1] = -1;
    StatusCode foundTower = StatusCode::FAILURE;
    for (tower=0;tower<m_numX*m_numY;++tower) {
        idents::TowerId t = idents::TowerId(tower);
        idents::VolumeIdentifier volId;
        volId.init(0,0);
        volId.append(m_volId_tower[tower]);
        int tray;
        int botTop;
        layerToTray(0, 0,tray, botTop);
        if(tray==-1) layerToTray(0, 1, tray, botTop);
        // get the right combination for this plane...
        volId.append(tray);
        volId.append(0);
        volId.append(botTop);
        volId.append(0); volId.append(0); // ladder and wafer
        sc = m_pDetSvc->getTransform3DByID(volId,&T);
        if (sc.isSuccess()) {
            foundTower = StatusCode::SUCCESS;
            // find the lowest and highest towers in each direction
            m_xLim[0] = std::min(m_xLim[0], t.ix());
            m_xLim[1] = std::max(m_xLim[1], t.ix());
            m_yLim[0] = std::min(m_yLim[0], t.iy());
            m_yLim[1] = std::max(m_yLim[1], t.iy());
        }          
    }
    return foundTower;
}

void TkrGeometrySvc::getTowerType()
{
    // use the list of existing towers to generate the tower type of each tower.
    // tower type is number of edges not touching another tower (0-4);
    int ix, iy;
    int tower;
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
}

void TkrGeometrySvc::makeTowerIds()
{
    int tower;
    for(tower=0;tower<m_numX*m_numY;++tower) {
        idents::VolumeIdentifier vId;
        vId.append(0);               // in Tower
        idents::TowerId t(tower);  
        vId.append(t.iy());          // yTower
        vId.append(t.ix());          // xTower
        vId.append(1);               // Tracker
        m_volId_tower[tower].init(0,0);  
        m_volId_tower[tower].append(vId);
    }
}

void TkrGeometrySvc::makeLayerIds()
{
    int bilayer, view;
    for(bilayer=0;bilayer<numLayers();++bilayer) {
        for (view=0; view<NVIEWS; ++view) {
            int tray;
            int botTop;            
            m_volId_layer[bilayer][view].init(0,0);
            layerToTray(bilayer, view, tray, botTop);
            if (tray==-1) continue;
            idents::VolumeIdentifier vId;
            vId.append(tray);
            vId.append(view);
            vId.append(botTop);
            vId.append(0); vId.append(0); // ladder and wafer

            m_volId_layer[bilayer][view].append(vId);
        }
    }  
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
    //bottom.append(1); // TKR
    bottom.append(0);                // tray 0
    idents::VolumeIdentifier idBot;
    HepTransform3D botTransform;
    for (int view = 0; view<2; ++view) {
        idBot = bottom;
        idBot.append(view);          // try both views
        idBot.append(1);             // top silicon (*most* bottom trays have one!)
        idBot.append(0); idBot.append(0);  // ladder, wafer
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

    IPropagator* track = m_G4PropTool;
    track->setStepStart(startPoint, startDir);
    // if the top were at the top of the stayclear, stayClear would be enuf
    // 100 mm gets you from the middle of the bottom closeout to below the tracker
    double propRange = stayClear+100.;
    double propBot = propTop - propRange;
    track->step(propRange);
    log << MSG::INFO  << "Propagator goes from "<< propTop << " to " << propBot << std::endl;

    int numSteps = track->getNumberSteps();
    int istep;
    idents::VolumeIdentifier id;
    idents::VolumeIdentifier prefix = m_pDetSvc->getIDPrefix();

    // now do the layer stuff

    double radlen;
    bool startCount = false;
    for (istep=0; istep<numSteps; ++istep) {
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
            if (item==2 || item==1) {
                startCount = true;
            }
            if (!startCount) continue;
            // only the "top" silicon in each layer is in the same tray as layer
            // the "bottom" silicon is one tray up
            int layer = (item==2 || item==6 || item==0) ? tray-1 : tray;
            if (item==2) { // converter
                m_radLenConv[layer] = radlen;
                m_convZ[layer] = stepPoint.z();
            } else {
                m_radLenRest[layer] += radlen;
            }
        }
    }

    log << MSG::DEBUG;
    int numAll = 0;
    int i, ind;
    for (i=0; i<NLAYERS; ++i) {
        convType type = getLayerType(i);
        if (type==ABSENT) continue;
        ind = (int) type;
        radlen = m_radLenConv[i];
        m_numLayers[ind]++;
        m_aveRadLenConv[ind] += radlen;
        m_aveRadLenRest[ind] += m_radLenRest[i];
        numAll++;
        m_aveRadLenConv[(int)ALL] += radlen;
        m_aveRadLenRest[(int)ALL] += m_radLenRest[i];

        if(log.isActive()) {
            log << " digiLayer " << i << " Conv " << m_radLenConv[i] << " Rest " 
                << m_radLenRest[i] << endreq;
        }
    }
    log << endreq;

    if(numAll!=m_numLayers[(int)ALL]) {
        log<< MSG::ERROR << "Inconsistent layer accounting, " 
            << numAll << " vs " <<m_numLayers[(int)ALL] << endreq;
        return StatusCode::FAILURE;
    } else {
        m_numLayers[(int)ALL] = numAll;
    }


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

double TkrGeometrySvc::getLayerZ(int digiLayer, int view) const
{
    // Purpose: returns the z for a given plane and view
    // Method:  accesses m_planeZ
    // Inputs   layer and optional view (if absent, return average)
    // Outputs: z position
    // Caveats:

    switch (view) {
    case 0:
        return m_planeZ[m_layerToPlane[digiLayer][view]];
        break;
    case 1:
        return m_planeZ[m_layerToPlane[digiLayer][view]];
        break;
    default:
        return 0.5*(m_planeZ[m_layerToPlane[digiLayer][0]] 
        + m_planeZ[m_layerToPlane[digiLayer][1]]);
    }
}


double TkrGeometrySvc::getReconLayerZ(int layer, int view) const
{
    // Purpose: returns the z for a given plane and view
    // Method:  accesses m_planeZ
    // Inputs   layer and optional view (if absent, return average)
    // Outputs: z position
    // Caveats:

    return getLayerZ(reverseLayerNumber(layer), view);
}

convType TkrGeometrySvc::getReconLayerType(int layer) const
{
    // Purpose: returns the converter type for a layer
    // Method:  reverses layer, and calls getDigiLayerType()
    // Inputs   layer
    // Outputs: layer type
    // Caveats:

    return getLayerType(reverseLayerNumber(layer));
}

convType TkrGeometrySvc::getLayerType(int digiLayer) const
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

    //get the top and bottom of the CAL crystals... why do I need to do this?

    int count;

    // top layer of the Cal
    // Which id is legal depends on whether the cal geometry is plain or segvols.
    // for the moment, I'll just check a few until I get a good one, and if I don't then
    // I bail... this is a bit of a kludge... seems to work though

    // should check for x or y, and pick the highest
    // also some sensible action if there is no cal, like set it to lowest tkr plane - 1.

    idents::VolumeIdentifier topLayerId;
    topLayerId.append(m_testTowerId);
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

StatusCode TkrGeometrySvc::getVolumeInfo() 
{
    StatusCode sc = StatusCode::SUCCESS;

    bool found = false;

    HepTransform3D T;
    int tower;

    //move  thru the planes in order, find their z postions, and assign them to layers
    int tray, botTop, layer, view;
    int lastTray, lastFace, lastPlane;
    int plane = 0;

    tower = m_testTower;
    idents::VolumeIdentifier vId, vId1, vId2, vIdTest;
    vId.init(0,0);
    vId.append(0);               // in Tower
    idents::TowerId t(tower);  
    vId.append(t.iy());          // yTower
    vId.append(t.ix());          // xTower
    // would be better to have this, but need to check if it will work
    vIdTest = vId;
    vId.append(1);             // tracker
    for(tray=0; tray<NLAYERS+1; ++tray) {
        vId1 = vId;
        vId1.append(tray);
        for(botTop=0;botTop<2;++botTop) { // sequential order: tray, botTop
            for (view=0;view<2;++view) {
                vId2 = vId1;
                vId2.append(view);
                vId2.append(botTop);
                vId2.append(0); vId2.append(0);
                sc = m_pDetSvc->getTransform3DByID(vId2,&T);
                if (sc.isSuccess()) {
                    found = true;
                    m_planeZ[plane] = (T.getTranslation()).z();
                    m_planeToView[plane] = view;
                    lastTray = tray;
                    lastFace = botTop;
                    lastPlane = plane;
                    ++plane; // just count up from the bottom, numbering planes sequentially
                }
            }
        }
    }
    m_topTrayNumber = ( (lastFace==0 && lastTray>0) ? lastTray : -1);
    m_numPlanes     = plane; // it's one more than the last plane number!
    
    //Figure out the layers by finding the adjacent planes.
    //  "Adjacent" means the distance is less than the layerSeparation
    layer = 0;
    double layerSeparation = 10.0;
    for (plane=0;plane<=lastPlane;++plane) {
       view = m_planeToView[plane];
       m_planeToLayer[plane]       = layer;
       m_layerToPlane[layer][view] = plane;
       if (plane==lastPlane) break;
       double thisZ = m_planeZ[plane];
       if(fabs(thisZ-m_planeZ[plane+1])>layerSeparation) {
            ++layer;
            m_isTopPlaneInLayer[plane] = true;
        }
    }
    // special handling for very top layer
    m_isTopPlaneInLayer[lastPlane] = (m_planeToLayer[lastPlane]==m_planeToLayer[lastPlane-1]);

    m_numLayers[ALL]  = ++layer;
    return StatusCode::SUCCESS;
}


StatusCode TkrGeometrySvc::getTestTower() 
{
    // finds the test tower, and checks for bottom tray
    StatusCode sc = StatusCode::SUCCESS;

    bool found = false;

    HepTransform3D T;
    int tower;

    int tray, botTop, view;
    int layer0 = -1;

    for(tower=0;tower<m_numX*m_numY;++tower) {
        idents::VolumeIdentifier vId, vId1, vId2, vId3, vIdTest;
        vId.init(0,0);
        vId.append(0);               // in Tower
        idents::TowerId t(tower);  
        vId.append(t.iy());          // yTower
        vId.append(t.ix());          // xTower
        // would be better to have this, but need to check if it will work
        vIdTest = vId;
        vId.append(1);             // tracker
        for(tray=0; tray<NLAYERS+1; ++tray) { // what does it mean if the 1st layer fails???
            vId1 = vId;
            vId1.append(tray);
            for(view=0;view<2;++view) {
                vId2 = vId1;
                vId2.append(view);
                for(botTop=0;botTop<2;++botTop) {
                    vId3 = vId2;
                    vId3.append(botTop);
                    vId3.append(0); vId3.append(0);
                    sc = m_pDetSvc->getTransform3DByID(vId3,&T);
                    if (sc.isSuccess()) {
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (found) {
            m_testTower    = tower;
            m_testTowerId  = vIdTest;
            sc = StatusCode::SUCCESS;
            m_bottomTrayNumber = ( (tray==0 && botTop==1) ? tray : -1);
            break;
        }       
    }
    return sc;
}

int TkrGeometrySvc::getPlaneSeparation(const idents::TkrId& id1, const idents::TkrId& id2) const
{
    // returns number of planes between two objects specified by TkrIds.
    //   -1 means number cannot be determined

    // Really don't even need tower values. A TkrId with only tray and botTop will have
    //    a plane associated with it.

    int nDiff = -1;
    if (!id1.hasTray() || !id1.hasBotTop() || !id2.hasTray() || !id2.hasBotTop()) return nDiff;

    int plane1 = getPlane(id1);
    int plane2 = getPlane(id2);
    return abs(plane1 - plane2);
}

int TkrGeometrySvc::getPlane(double z) const
{
    if(z >= m_planeZ[m_numPlanes-1]) return m_numPlanes-1;

    int i;
    for(i=0; i<m_numPlanes-1; ++i) {
        if(z <= m_planeZ[i+1]) {
            return (fabs(m_planeZ[i] - z) < fabs(m_planeZ[i+1] - z)) ? i : i+1;
        }
    }
    return -1;  // this will never happen!!
}

// Good a place as any for this one...
unsigned int TkrGeometrySvc::getDefaultClusterStatus() const
{
    using namespace Event;
    int planeOffset = -trayToPlane(0, 0);
    int layerOffset = getLayer(1) - getLayer(0);
    return 
        (layerOffset>0 ? TkrCluster::maskLAYEROFFSET : 0) |
        (planeOffset>0 ? TkrCluster::maskPLANEOFFSET : 0) |
        (TkrCluster::maskVERSION&(TkrCluster::VERSION<<TkrCluster::shiftVERSION));
}

double TkrGeometrySvc::truncateCoord( double x, double pitch, 
                     int numElements, int& elementNumber, bool reverse) const
{
    //Returns an element number and the coordinates in the local element

    double xScaled = x/pitch;
    // this is the correction for odd number of elements (towers in EM, for example)
    double delta = 0.5*(numElements%2);
    double xMod = xScaled + delta;
    int theFloor = (int) floor(xMod);
    elementNumber = theFloor + numElements/2;
    // if it's outside the actual element, assign the closest
    // this will not happen for real hits, but may for extrapolated hits
    // or transformed hits.
    elementNumber = std::max(std::min(elementNumber, numElements-1),0);
    if (reverse) elementNumber = numElements - 1 - elementNumber;
    return pitch*(xMod - theFloor - 0.5);
}

bool TkrGeometrySvc::inTower(int view, const Point p, int& iXTower, int& iYTower, 
                             double& xActiveDist, double& yActiveDist, 
                             double& xGap, double& yGap) const
{
    double twrPitch = towerPitch();
    double numX = numXTowers();
    double numY = numYTowers();
    bool isInTower = true;
    double xTower = truncateCoord(p.x(), twrPitch, numX, iXTower);
    double yTower = truncateCoord(p.y(), twrPitch, numY, iYTower);
    // check if the point is in an inter-tower gap
    int nWafer = nWaferAcross();
    double xPitch, yPitch, xSiGap, ySiGap;
    double deadGap = siDeadDistance();
    if (view==0) {
        xPitch = ladderPitch();
        yPitch = waferPitch();
        xSiGap   = ladderGap();
        ySiGap   = ladderInnerGap();
    } else {
        yPitch = ladderPitch();
        xPitch = waferPitch();
        ySiGap   = ladderGap();
        xSiGap   = ladderInnerGap();
    }

    // if these are negative, track misses active area of tower
    // probably no point in constraining hit in this plane
    xGap = xSiGap + 2*deadGap;
    yGap = ySiGap + 2*deadGap;
    xActiveDist = fabs(xTower) - nWafer*xPitch + xGap; 
    yActiveDist = fabs(yTower) - nWafer*yPitch + yGap;
    if (xActiveDist>0 && yActiveDist>0) { // test for "inside active Tower"
        // look for internal gaps
        double xWafer, yWafer;
        int iXWafer, iYWafer;
        xWafer = truncateCoord(xTower, xPitch, nWafer, iXWafer);
        yWafer = truncateCoord(yTower, yPitch, nWafer, iYWafer);

        double activeWaferSide = siActiveWaferSide();
        xActiveDist = activeWaferSide - fabs(xWafer);
        yActiveDist = activeWaferSide - fabs(yWafer);
    } else {
        // towerGap is big, so no point in keeping 2 versions
        double towerGap    = twrPitch - nWafer*xPitch + xGap;

        // note that if we are on the outside of an edge tower the "gap" 
        // is infinite, but the simple intertower gap is effectively infinite
        // as far as the fit is concerned, so it shouldn't matter.
        isInTower = false;
        if (xActiveDist<0 ) {
            yGap = twrPitch;
            xGap = (yActiveDist<0 ? twrPitch : towerGap);
        } else {
            xGap = twrPitch;
            yGap = towerGap;
        }
    }
    return isInTower;
}
