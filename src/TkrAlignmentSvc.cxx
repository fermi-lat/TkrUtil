
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Service.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "src/TkrAlignmentSvc.h"

#include "idents/TowerId.h"
#include "idents/VolumeIdentifier.h"

#include <fstream>
#include <algorithm>

#include "xml/IFile.h"

static const SvcFactory<TkrAlignmentSvc> s_factory;
const ISvcFactory& TkrAlignmentSvcFactory = s_factory;


// Service parameters which can be set at run time must be declared.
// This should be done in the constructor.

TkrAlignmentSvc::TkrAlignmentSvc(const std::string& name, 
                                 ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the file to get data from   
    declareProperty("simFile", m_simFile="");
    declareProperty("recFile", m_recFile="");
    declareProperty("testMode", m_testMode=0);
    declareProperty("maximumDelta", m_maxDelta=5.0);
    
    return;
}


StatusCode TkrAlignmentSvc::initialize()
{
    // Purpose: reads in (or set) alignment constants
    // Inputs:  None
    // Outputs: Status code (Success/Failure)
    
    MsgStream log(msgSvc(), name());
    // log.setLevel(MSG::DEBUG);

    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    
    m_simFile = "";
    m_recFile = "";
    m_testMode = 0;
    m_fileFlag = 0;
    
    setProperties();

    // we need the detModel and geometry service

    sc = service("TkrGeometrySvc", m_pGeoSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrGeometrySvc" 
            << endreq;
        return sc;
    }

    sc = service("GlastDetSvc", m_pDetSvc, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get GlastDetSvc" 
            << endreq;
        return sc;
    }

    if (getGeometry().isFailure()) {
        log << MSG::ERROR << "Alignment geometry not created" << endreq;
        return StatusCode::FAILURE;
    }
    
    // initialize alignment constants   
    int i;
    AlignmentConsts alConsts(0., 0., 0., 0., 0., 0.);
    for (i=0;i<NELEMENTS; ++i) {
        m_simConsts[i] = alConsts;
        m_recConsts[i] = alConsts;
    }

    if(!m_testMode) {        
        // If there are no alignment files, service will do nothing
        if (m_simFile=="" && m_recFile=="") {        
            log << MSG::INFO << "No alignment was requested." << endreq;
            log << MSG::INFO << "  No alignment will be done." << endreq;
            return sc;
        }
    } else {
        if (doTest().isFailure()) { return StatusCode::FAILURE;}
    }

    // get the data

    m_mode = "sim";
    if(getData(m_simFile).isFailure()) {
        log << MSG::ERROR << "Simulation Alignment constants not read" << endreq;
        return StatusCode::FAILURE;
    }
    m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);


    m_mode = "rec";
    if(getData(m_recFile).isFailure()) {
        log << MSG::ERROR << "Reconstruction Alignment constants not read" << endreq;
        return StatusCode::FAILURE;
    }
    m_fileFlag = m_fileFlag|(1<<REC_SHIFT);

    return sc;
}

StatusCode TkrAlignmentSvc::finalize()
{
    return StatusCode::SUCCESS;
}

StatusCode TkrAlignmentSvc::getGeometry()
{
    using namespace idents;

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    // Get the z for each tray

    VolumeIdentifier baseVol; // prefix for all trays
    baseVol.init(0,0);
    baseVol.append(0);
    baseVol.append(0);
    baseVol.append(0);
    baseVol.append(1);

    int tray, layer, view;
    double zTop, zBot;

    TkrAlignmentGeomVisitor* visitor = new TkrAlignmentGeomVisitor();

    std::cout << std::endl;

    m_pDetSvc->accept(*visitor);

    for (tray =1; tray<m_pGeoSvc->numLayers(); ++tray) {
        m_pGeoSvc->trayToLayer(tray,0,layer,view);
        zBot = m_pGeoSvc->getReconLayerZ(m_pGeoSvc->reverseLayerNumber(layer), view);
        m_pGeoSvc->trayToLayer(tray,1,layer,view);
        zTop = m_pGeoSvc->getReconLayerZ(m_pGeoSvc->reverseLayerNumber(layer), view);
        m_trayZ[tray] = 0.5*(zBot+zTop);
        m_planeZ[tray][0] = zBot - m_trayZ[tray];
        m_planeZ[tray][1] = zTop - m_trayZ[tray];
    }

    double siThickness = m_pGeoSvc->siThickness();

    //bottom tray
    double trayBotHeight = visitor->getTrayBotHeight();
    m_pGeoSvc->trayToLayer(0,1,layer,view);
    zTop = m_pGeoSvc->getReconLayerZ(m_pGeoSvc->reverseLayerNumber(layer), view); 
    m_trayZ[0] = zTop + 0.5*siThickness - 0.5*trayBotHeight;
    m_planeZ[0][1] = 0.5*(trayBotHeight - siThickness);
    m_planeZ[0][0] = - 0.5*trayBotHeight;

    //top tray
    double trayTopHeight = visitor->getTrayTopHeight();
    tray = m_pGeoSvc->numLayers();
    m_pGeoSvc->trayToLayer(tray,0,layer, view);
    zBot = m_pGeoSvc->getReconLayerZ(m_pGeoSvc->reverseLayerNumber(layer), view);
    m_trayZ[tray] = zBot - 0.5*siThickness + 0.5*trayTopHeight;
    m_planeZ[tray][0] = -0.5*(trayTopHeight - siThickness);
    m_planeZ[tray][1] = 0.5*trayTopHeight;
    
    log << MSG::DEBUG ;
    if (log.isActive()) {
        for(tray = 0; tray<m_pGeoSvc->numLayers()+1; ++tray) {
            std::cout << "tray " << tray << " z: " << m_trayZ[tray] << " "
                << m_planeZ[tray][0] << " " << m_planeZ[tray][1] << std::endl;
        }        
        std::cout << std::endl;
          
        for (layer = 0; layer<m_pGeoSvc->numLayers(); ++layer ) {
            for (view = 0; view<2; ++view) {
                std::cout << "rlayer/view z " << layer << " " << view << " "
                    << m_pGeoSvc->getReconLayerZ(m_pGeoSvc->reverseLayerNumber(layer), view) 
                    << std::endl;
            }
        }
    }
    log << endreq;
    
    // clean up
    delete visitor;

    return sc;   
}
    
StatusCode TkrAlignmentSvc::getData(std::string fileName)
{
    MsgStream log(msgSvc(), name());
    
    StatusCode sc = StatusCode::SUCCESS;
    
    if( fileName == "") return sc;
    
    // this method resolves environmental variables in the file name
    xml::IFile::extractEnvVar(&fileName);    
    log << MSG::INFO << "Input file for " << m_mode << " alignment: " << endreq
        << "    " << fileName << endreq;
    
    std::ifstream theFile;
    theFile.open( fileName.c_str());
    m_dataFile = &theFile;
    
    if (!theFile) {
        log << MSG::ERROR 
            << "Alignment file " << fileName << " not found: check jobOptions." 
            << endreq;
        return StatusCode::FAILURE;
    }
    log << MSG::DEBUG;
    if (log.isActive() ) {
        log << "Do " << m_mode << " alignment constants:";
    }
    log << endreq;
    
    readFromFile();
    theFile.close();
    
    log << MSG::DEBUG ;
    if (log.isActive()) {
        for (int tower = 0; tower < NTOWERS; ++tower) {
            bool first = true;
            for (int layer=0; layer<NLAYERS; ++layer) {
                int index = getIndex(tower, layer, 0, 0, 0);
                if (m_mode=="sim") {
                    if ( !m_simConsts[index].isNull() ) {
                        if (first) std::cout << m_mode <<" consts for tower " << tower 
                            << ", view 0, ladder 0, wafer 0" << std::endl;
                        first = false;
                        std::cout << "layer " << layer << " ";
                        m_simConsts[index].fillStream(std::cout);
                    }
                }   else  {             
                    if ( !m_recConsts[index].isNull() ) {
                        if (first) std::cout << m_mode <<" consts for tower " << tower 
                            << ", view 0, ladder 0, wafer 0" << std::endl;
                        first = false;
                        std::cout << "layer " << layer << " ";
                        m_recConsts[index].fillStream(std::cout);
                    }
                }
            }
        }
    }
    log << endreq;
   
    return sc;
}

StatusCode TkrAlignmentSvc::doTest()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "Test mode " << m_testMode << 
        " requested, each hit will be translated by +1 strip in Global X" << endreq
        << "               and +2 strips in Global Y";
    AlignmentConsts alConsts(0.22, 0.44, 0., 0., 0., 0.);

    int i;

    if (m_testMode&(1<<SIM_SHIFT)) {
        m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);
        log << " in Sim" ;
        for (i=0; i<NELEMENTS; ++i) {m_simConsts[i] = alConsts;
        }
    }
    if (m_testMode&(1<<REC_SHIFT)) {
        m_fileFlag = m_fileFlag|(1<<REC_SHIFT);
        log << " Rec";
        for (i=0; i<NELEMENTS; ++i) {m_recConsts[i] = alConsts;
        }
    }
    log << endreq;
    
    return sc;
}    

StatusCode TkrAlignmentSvc::readFromFile()
{    
    // Purpose:
    // Inputs:  File name
    // Outputs: None
    // Dependencies: None
    // Caveats: None
    
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    std::string junk;
    std::string flag;
    
    while(!m_dataFile->eof()) {
        *m_dataFile >> flag ;
        if (m_dataFile->eof()) break;

        if(flag== "//") { // comment line, just skip 
            std::getline(*m_dataFile, junk);
        } else {
            if (flag=="tower") {
                
                *m_dataFile >> m_tower;
                
                if(m_tower>=m_pGeoSvc->numXTowers()*m_pGeoSvc->numYTowers()) {
                    return StatusCode::FAILURE;
                }
                if (m_dataFile->eof()) break;
                
                double a,b,c,d,e,f;       
                *m_dataFile >> a >> b >> c >> d >> e >> f ;
                
                m_towerConsts = AlignmentConsts(0.001*a, 0.001*b, 0.001*c,
                    0.001*d, 0.001*e, 0.001*f);
                
                if (m_towerConsts.isNull()) { continue; }
                
                
                if (fillTrayConsts().isFailure()) {
                    log << MSG::ERROR << "fillTrayConsts failed!" << endreq;
                    return StatusCode::FAILURE;
                } 
            } else {
                std::getline(*m_dataFile, junk);
            }
        }
    }       
    return sc;
}

StatusCode TkrAlignmentSvc::fillTrayConsts()
{

    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    AlignmentConsts& towerConsts = m_towerConsts;

    int nTray = m_pGeoSvc->numLayers()+1;
    
    // could read this in here
    AlignmentConsts thisTray(0., 0., 0., 0., 0., 0.);

    int tray;

    for (tray = 0; tray< nTray; ++tray) {
        double deltaX = thisTray.getDeltaX() + towerConsts.getDeltaX()
            + towerConsts.getRotY()*m_trayZ[tray];
        double deltaY = thisTray.getDeltaY() + towerConsts.getDeltaY()
            - towerConsts.getRotX()*m_trayZ[tray];
        double deltaZ = thisTray.getDeltaZ() + towerConsts.getDeltaZ();
        double rotX   = thisTray.getRotX()   + towerConsts.getRotX();
        double rotY   = thisTray.getRotY()   + towerConsts.getRotY();
        double rotZ   = thisTray.getRotZ()   + towerConsts.getRotZ();
        m_tray = tray;
        m_trayConsts = AlignmentConsts(deltaX, deltaY, deltaZ,
            rotX, rotY, rotZ);
        if(fillLadderConsts().isFailure()) {return StatusCode::FAILURE;}
    }
    return sc;
}

StatusCode TkrAlignmentSvc::fillLadderConsts()
{
    
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    AlignmentConsts& trayConsts = m_trayConsts;
    
    int nPlane = m_pGeoSvc->numViews();
    int nLadder = m_pGeoSvc->nWaferAcross();
    
    // could read this in here
    AlignmentConsts thisPlane(0., 0., 0., 0., 0., 0.);
    
    int plane, layer, view, ladder;
    
    double xPlane, yPlane;
    double trayWidth   = m_pGeoSvc->trayWidth();
    double ladderGap   = m_pGeoSvc->ladderGap();
    // why don't we ask TkrGeometrySvc to calculate this?
    double ladderWidth = (trayWidth - (nLadder-1)*ladderGap)/nLadder;
    double offset      = -0.5*trayWidth + 0.5*ladderWidth;
    
    for (plane=0; plane< nPlane; ++plane) {
        m_pGeoSvc->trayToLayer(m_tray, plane, layer, view);
        if (layer<0 || layer>=m_pGeoSvc->numLayers()) { continue;}
        
        double zPlane = m_planeZ[m_tray][plane];
        for (ladder=0; ladder<nLadder; ++ladder) {
            if (view==0) {
                xPlane = offset + ladder*(ladderWidth + ladderGap);
                yPlane = 0.;
            } else {
                xPlane = 0.;
                yPlane = offset + ladder*(ladderWidth + ladderGap);
            }
            
            double deltaX = thisPlane.getDeltaX() + trayConsts.getDeltaX()
                + trayConsts.getRotY()*zPlane - trayConsts.getRotZ()*yPlane;
            double deltaY = thisPlane.getDeltaY() + trayConsts.getDeltaY()
                - trayConsts.getRotX()*zPlane + trayConsts.getRotZ()*xPlane;
            double deltaZ = thisPlane.getDeltaZ() + trayConsts.getDeltaZ()
                +trayConsts.getRotX()*yPlane - trayConsts.getRotY()*xPlane;
            double rotX = thisPlane.getRotX() + trayConsts.getRotX();
            double rotY = thisPlane.getRotY() + trayConsts.getRotY();
            double rotZ = thisPlane.getRotZ() + trayConsts.getRotZ();
            
            m_ladder = ladder;
            m_botTop = plane;
            m_ladderConsts = AlignmentConsts(deltaX, deltaY, deltaZ,
                rotX, rotY, rotZ);
            
            if(fillWaferConsts().isFailure()) { return StatusCode::FAILURE;}
        }
    }
    return sc;
}

StatusCode TkrAlignmentSvc::fillWaferConsts()
{   
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    AlignmentConsts& ladderConsts = m_ladderConsts;

    int layer, view, wafer;
    int nWafer = m_pGeoSvc->nWaferAcross();
    
    // could read this in here
    AlignmentConsts thisWafer(0., 0., 0., 0., 0., 0.);
    
    m_pGeoSvc->trayToLayer(m_tray, m_botTop, layer, view);

    double trayWidth   = m_pGeoSvc->trayWidth();
    double waferGap    = m_pGeoSvc->ladderInnerGap();
    double waferWidth  = (trayWidth - (nWafer-1)*waferGap + .1)/nWafer;
    double trayWidth1   = nWafer*waferWidth + (nWafer-1)*waferGap;
    double laddergap   = m_pGeoSvc->ladderGap();
    double offset      = -0.5*trayWidth1 + 0.5*waferWidth;

    double xWafer, yWafer;

    for (wafer=0; wafer<nWafer; ++wafer) { 
        if (view==0) {
            xWafer = 0.;
            yWafer = offset + wafer*(waferWidth+waferGap);
        } else {
            xWafer = offset + wafer*(waferWidth+waferGap);
            yWafer = 0.;
        }
        
        double deltaX = thisWafer.getDeltaX() + ladderConsts.getDeltaX()
             - ladderConsts.getRotZ()*yWafer;
        double deltaY = thisWafer.getDeltaY() + ladderConsts.getDeltaY()
            + ladderConsts.getRotZ()*xWafer;
        double deltaZ = thisWafer.getDeltaZ() + ladderConsts.getDeltaZ()
            +ladderConsts.getRotX()*yWafer - ladderConsts.getRotY()*xWafer;
        double rotX = thisWafer.getRotX() + ladderConsts.getRotX();
        double rotY = thisWafer.getRotY() + ladderConsts.getRotY();
        double rotZ = thisWafer.getRotZ() + ladderConsts.getRotZ();

        int index = getIndex(m_tower, layer, view, m_ladder, wafer);
        if (m_mode=="sim") {
            m_simConsts[index] = AlignmentConsts(deltaX, deltaY, deltaZ,
                rotX, rotY, rotZ);
        } else {
            m_recConsts[index] = AlignmentConsts(deltaX, deltaY, deltaZ,
                rotX, rotY, rotZ);
        }
    }
    return sc;
}

int TkrAlignmentSvc::getIndex(int tower, int layer, 
                              int view, int ladder, int wafer) const
{
    // Purpose:  calculate index into array of consts
    // Inputs:   tower, bilayer, view, ladder, wafer
    // Outputs:  index
    
    int index = -1;
    if (layer<0 || layer>=NLAYERS || tower<0 || tower>=NTOWERS ||
        view<0 || view>=NVIEWS || ladder<0 || ladder >= NLADDERS ||
        wafer<0 || wafer>=NWAFERS)
    {return index;}
    // for now, hardwired to be as large as will ever by needed
    return wafer + NWAFERS*(ladder + NLADDERS*(view + NVIEWS*(layer + NLAYERS*tower)));
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, int tower, int layer, 
                                            int view, int ladder, int wafer) const
{
    // Purpose:  return pointer to an set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set
    
    int index = getIndex(tower, layer, view, ladder, wafer);
    
    return getConsts(type, index);
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, 
                                                  idents::VolumeIdentifier id) const
{
    // Purpose:  return pointer to a set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set
    
    int tower  = idents::TowerId(id[2],id[1]).id();
    int tray   = id[4];
    int botTop = id[6];
    int layer, view;
    m_pGeoSvc->trayToLayer(tray, botTop, layer, view);
    view       = id[5];
    int ladder = id[7];
    int wafer  = id[8];

    int index = getIndex(tower, layer, view, ladder, wafer);
    
    return getConsts(type, index);
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, int index) const
{
    // Purpose:  return pointer
    // Inputs:   index
    // Outputs:  pointer to that set of consts
       
    if (index>=0 && index < NELEMENTS) {
        
        if (type==SIM) {
            return &m_simConsts[index];
        } else {
            return &m_recConsts[index];
        }
        
    } else {
        AlignmentConsts* result = 0;
        return result;
    }
}

void TkrAlignmentSvc::moveMCHit(idents::VolumeIdentifier id, HepPoint3D& entry,
                                HepPoint3D& exit) const
{
    // Purpose:     Move an McHit according to alignment constants
    // Inputs:      volId and entry and exit points
    // Output:      modified entry and exit point arguments
    
    //This is called "dir" but what I'm calculating is a direction normalized to dir.z() = 1;
    //  This yields the slopes in x and y as the x() and y() components.
    HepVector3D dir = (exit-entry);
    double delz = dir.z();
	
    // This doesn't work in some flavors of unix, so rewrite below
	//dir = (delz != 0) ? dir/delz : HepVector3D( 0., 0., 1.);

	HepVector3D tempVec = dir; // just being superstitious here
	if (delz!=0) { dir = tempVec/delz;}
	else         { dir = HepVector3D(0., 0., 1.);}    
    
    // calculation is done separately for entry and exit, because transformation depends
    //    in 2nd order (~ microns) on the coordinates. This can be speeded up if
    //    necessary, by basing the calculation on the average coordinates
    
    const AlignmentConsts* alConsts = getConsts(SIM, id);

    int view = id[5];
    
    HepVector3D deltaEntry = getDelta(view, entry, dir, alConsts);
    HepVector3D deltaExit  = getDelta(view, exit,  dir, alConsts);

    // for now, limit delta to 5 mm (modifiable in jobOptions
    // later fix transformation for special cases

    double mag;
    HepVector3D dirDelta;

    mag = std::min(m_maxDelta, deltaEntry.mag());
    dirDelta = deltaEntry.unit();
    deltaEntry = mag*dirDelta;

    mag = std::min(m_maxDelta, deltaExit.mag());
    dirDelta = deltaExit.unit();
    deltaExit = mag*dirDelta;

    entry = entry + deltaEntry;
    exit  = exit  + deltaExit;
    
    //std::cout << "moveMc " << view << " delta " << deltaEntry.x() << " " << deltaEntry.y() 
    //    << " dir " << dir.x() << " "<< dir.y() << std::endl;
    
    return;
}

void TkrAlignmentSvc::moveCluster(int tower, int layer, int view, int ladder,
                                  HepPoint3D& point) const
{
    // Purpose:  Move an cluster according to alignment constants (dx and dy only)
    // Inputs:   tower, layer, view, ladder, and cluster position
    // Output:   modified position argument
    
    int wafer = 0;
    const AlignmentConsts* alConsts = getConsts(REC, tower, layer, view, ladder, wafer);
    // "x" is the only thing we can correct at this stage...
    //     average for all the wafers in a ladder would be better.
    AlignmentConsts alConsts1(alConsts->getDeltaX(), alConsts->getDeltaY()); 
    HepVector3D dir(0., 0., 1);

    view = -1;

    point = point - getDelta(view, point, dir, &alConsts1);
}

void TkrAlignmentSvc::moveReconHit(int /*tower*/, int /*layer*/, int /*view*/, int /*ladder*/,
                                   HepPoint3D& /*point*/, HepVector3D /*dir*/) const
{
    //placeholder for now
}

HepVector3D TkrAlignmentSvc::getDelta(int view, const HepPoint3D& point,
                                      const HepVector3D& dir,
                                      const AlignmentConsts* alConsts) const 
{
    // Purpose:  Calculate translation in x,y for fixed z
    // Inputs:   view, point and direction
    // Output:   vector of translation
    
    double small = 1.e-3;
    
    double point1 = point.x(); 
    double point2 = point.y(); 
    //double pointZ = point.z(); not yet used
    
    double deltaX = alConsts->getDeltaX();
    double deltaY = alConsts->getDeltaY();
    double deltaZ = alConsts->getDeltaZ();
    
    double rotX = alConsts->getRotX();
    double rotY = alConsts->getRotY();
    double rotZ = alConsts->getRotZ();
    
    double dirZ   = dir.z();
    if (fabs(dirZ)<small) dirZ = small;
    double alpha1 = dir.x()/dirZ;
    double alpha2 = dir.y()/dirZ;

    double pointX, pointY;
    double alphaX, alphaY;

    if (view==0) {
        pointX = point1;
        pointY = point2;
        alphaX = alpha1;
        alphaY = alpha2;
    } else {
        pointX = point2;
        pointY = -point1;
        alphaX = alpha2;
        alphaY = -alpha1;
    }
  
    double rotTerm = deltaZ + rotX*pointY - rotY*pointX; 
    
    double deltaPointX = - deltaX + rotZ*pointY
        + alphaX*rotTerm;
    double deltaPointY = - deltaY - rotZ*pointX
        + alphaY*rotTerm;

    if (view==0) {
        return HepVector3D(deltaPointX,  deltaPointY, 0.);
    } else {
        return HepVector3D(deltaPointY, -deltaPointX, 0.);
    }
}


IGeometry::VisitorRet TkrAlignmentGeomVisitor::pushShape(ShapeType /* s */, const UintVector& idvec, 
        std::string name, std::string /* material*/, const DoubleVector& params, 
        VolumeType /*type*/)
{
    if(name=="oneCAL") {
        return AbortSubtree;
    }
    
    if (name=="towerRow") {
        if(idvec[0]>0) return AbortSubtree;
    }
    if (name=="oneTower") { 
        if(idvec[0]>0) return AbortSubtree;
    }
    if (name=="trayBot") {m_trayBotHeight = params[8];}

    if (name=="trayTop") {m_trayTopHeight = params[8];}
    
    bool debug = false;
    if (debug) {
        if (name=="traySuper"||name=="trayTop"||name=="trayBot"
            ||name=="trayReg"||name=="trayNoConv"||name=="oneTKRStack"
            ||name=="oneTKR") { 
            std::cout << name << ": z, height " 
                << params[2] << ", " << params[8] << std::endl;
        }
    }
    if (name=="trayTop") return AbortSubtree;

    return More;
}

void TkrAlignmentGeomVisitor::popShape() { return;}


// queryInterface

StatusCode  TkrAlignmentSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrAlignmentSvc == riid) {
        *ppvIF = dynamic_cast<ITkrAlignmentSvc*> (this);
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
    return StatusCode::SUCCESS;
}

// access the type of this service

const IID&  TkrAlignmentSvc::type () const {
    return IID_ITkrAlignmentSvc;
}
