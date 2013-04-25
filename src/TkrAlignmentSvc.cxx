/**
@file TkrAlignmentSvc.cxx

@brief handles Tkr alignment
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrAlignmentSvc.cxx,v 1.49 2013/04/10 23:15:43 lsrea Exp $
*/

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Service.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "src/TkrAlignmentSvc.h"

#include "idents/TowerId.h"
#include "idents/VolumeIdentifier.h"

#include "facilities/Util.h"

//#ifdef DEFECT_NO_STRINGSTREAM
//#include <strstream>
//#else
#include <sstream>
//#endif

#include <fstream>
#include <algorithm>
#include <string>
//do I really need these?
#include <string.h>
#include <cctype>

//static const SvcFactory<TkrAlignmentSvc> s_factory;
//const ISvcFactory& TkrAlignmentSvcFactory = s_factory;
DECLARE_SERVICE_FACTORY(TkrAlignmentSvc);

// Service parameters which can be set at run time must be declared.
// This should be done in the constructor.

namespace{

    AlignmentConsts nullConsts(0.,0.,0.,0.,0.,0.);

}

TkrAlignmentSvc::TkrAlignmentSvc(const std::string& name, 
                                 ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    m_simScaleVec.assign(6,1.0);
    m_recScaleVec.assign(6,1.0);
   
    // Name of the file to get data from 
    // These constants ordinarily come from the calibration system now
    declareProperty("simFile", m_simFile="");
    declareProperty("recFile", m_recFile="");
    declareProperty("testMode", m_testMode=0);
    declareProperty("maximumDelta", m_maxDelta=5.0);
    declareProperty("printConsts",  m_printConsts=false);
    // these scaling factors are for testing (for use with calibration system)
    declareProperty("simScale",     m_simScale=1.0);
    declareProperty("recScale",     m_recScale=1.0);
    declareProperty("simScaleVec",  m_simScaleVec);
    declareProperty("recScaleVec",  m_recScaleVec);

    return;
}

StatusCode TkrAlignmentSvc::initialize()
{
    // Purpose: reads in alignment constants
    // Inputs:  None
    // Outputs: Status code (Success/Failure)

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();

    m_simFile = "";
    m_recFile = "";
    m_testMode = 0;
    m_fileFlag = 0;

    m_itemCol.clear();
    m_pItem = m_itemCol.begin();

    setProperties();

    // we need the detModel and geometry service

    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrGeometrySvc" 
            << endreq;
        return sc;
    }

    m_nTowers = m_tkrGeom->numXTowers()*m_tkrGeom->numYTowers();
    m_nLayers = m_tkrGeom->numLayers();
    m_nTrays  = m_tkrGeom->numTrays();
    m_nViews  = 2;
    m_nFaces  = 2;
    m_nLadders = m_tkrGeom->nWaferAcross();
    m_nWafers  = m_nLadders;

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
    AlignmentConsts alConsts(0., 0., 0., 0., 0., 0.);

    bool valid = m_simConsts.setDims(m_nTowers, m_nLayers, m_nViews, m_nLadders, m_nWafers);
    if (!valid) {
        log << MSG::ERROR << "simConsts not created" << endreq;
        return StatusCode::FAILURE;
    }
    m_simConsts.setRangeCheck(false);
    m_simConsts.setValue(alConsts);

    valid = m_recConsts.setDims(m_nTowers, m_nLayers, m_nViews, m_nLadders, m_nWafers);
    if (!valid) {
        log << MSG::ERROR << "recConsts not created" << endreq;
        return StatusCode::FAILURE;
    }
    m_recConsts.setRangeCheck(false);
    m_recConsts.setValue(alConsts);

    if(!m_testMode) {        
        // If there are no alignment files, the constructor does not get initial alignment constants.
        // However, constants can be received later through the update() method from the data base.
        if (m_simFile=="" && m_recFile=="") {        
            log << MSG::INFO << "no text style alignment files specified in job options" << endreq;
            return sc;
        }
    } else {
        if (doTest().isFailure()) { return StatusCode::FAILURE;}
    }

    // get the data

    setMode("sim");
    if(getData(m_simFile).isFailure()) {
        log << MSG::ERROR << "Simulation Alignment constants not read" << endreq;
        return StatusCode::FAILURE;
    }
    if (m_simFile!="") m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);

    setMode("rec");
    if(getData(m_recFile).isFailure()) {
        log << MSG::ERROR << "Reconstruction Alignment constants not read" << endreq;
        return StatusCode::FAILURE;
    }
    if (m_recFile!="") m_fileFlag = m_fileFlag|(1<<REC_SHIFT);

    m_useFlags = false;

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

    m_pDetSvc->accept(*visitor);

    double siThickness = m_tkrGeom->siThickness();
    int firstTray = 1;
    int endTray   = m_nTrays-1;

    bool valid = m_faceZ.setDims(m_nTrays, m_nFaces);
    if (!valid) return StatusCode::FAILURE;
    m_faceZ.setRangeCheck(false);
    m_faceZ.setValue(0.0);

    m_trayZ.resize(m_nTrays, 0.0);

    // some shenanigans to accommodate EGRET.
    // LAT has bot bottom and top trays: nTrays == nLayers + 1;
    // EGRET has neither bottom nor top tray: nTrays == nLayers;
    // I think there's an EM with a configuration that isn't covered here,
    //    but I don't care!

    if(m_tkrGeom->getBottomTrayFlag()==1) { // There's a bottom tray; handle it!
        double trayBotHeight = visitor->getTrayBotHeight();
        m_tkrGeom->trayToLayer(0,1,layer,view);
        zTop = m_tkrGeom->getLayerZ(layer, view); 
        m_trayZ[0] = zTop + 0.5*siThickness - 0.5*trayBotHeight;
        m_faceZ(0,0) = 0.5*(trayBotHeight - siThickness);
        m_faceZ(0,1) = - 0.5*trayBotHeight;
    } else { // otherwise include the lowest tray in the loop
        firstTray = 0;
        if(m_tkrGeom->getTopTrayFlag()==1) {endTray--;}
    }

    if(m_tkrGeom->getTopTrayFlag()==1) { // there's a top tray...
        double trayTopHeight = visitor->getTrayTopHeight();
        tray = m_nTrays - 1;
        m_tkrGeom->trayToLayer(tray,0,layer, view);
        zBot = m_tkrGeom->getLayerZ(layer, view);
        m_trayZ[tray] = zBot - 0.5*siThickness + 0.5*trayTopHeight;
        m_faceZ(tray,0) = -0.5*(trayTopHeight - siThickness);
        m_faceZ(tray,0) = 0.5*trayTopHeight;
    } 

    for (tray =firstTray; tray<endTray; ++tray) {
        m_tkrGeom->trayToLayer(tray,0,layer,view);
        zBot = m_tkrGeom->getLayerZ(layer, view);
        m_tkrGeom->trayToLayer(tray,1,layer,view);
        zTop = m_tkrGeom->getLayerZ(layer, view);
        m_trayZ[tray] = 0.5*(zBot+zTop);
        double halfHeight = 0.5*(zTop - zBot);
        m_faceZ(tray, 0) = -halfHeight;
        m_faceZ(tray, 1) =  halfHeight;
    }

    log << MSG::DEBUG ;
    if (log.isActive()) {
        for(tray = 0; tray<m_nTrays; ++tray) {
            std::cout << "tray " << tray << " z: " << m_trayZ[tray] << " "
                << m_faceZ(tray, 0) << " " << m_faceZ(tray, 1) << std::endl;
        }        
        std::cout << std::endl;

        for (layer = 0; layer<m_nLayers; ++layer ) {
            for (view = 0; view<m_nViews; ++view) {
                std::cout << "rlayer/view z " << layer << " " << view << " "
                    << m_tkrGeom->getLayerZ(layer, view) 
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

    int ret = facilities::Util::expandEnvVar(&fileName);

    if(ret>=0) {
        log << MSG::INFO << "Input file for " << m_mode << " alignment: " << endreq
            << "    " << fileName << endreq;
    } else {
        log << MSG::ERROR << "Input filename " << fileName << " not resolved" << endreq;
        return StatusCode::FAILURE;
    }

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

    if (readFromFile().isFailure()) {
        log << MSG::ERROR << "Failed to read input file" << endreq;
        return StatusCode::FAILURE;
    }
    theFile.close();

    log << MSG::DEBUG;
    if (log.isActive()) {
        m_pItem = m_itemCol.begin();
        for (; m_pItem<m_itemCol.end(); ++m_pItem) {
            AlignmentItem item = **m_pItem;
            log << itemType[item.getType()] << " " << item.getNumber() << " "  
                << item.getConsts().getDeltaX() << " " 
                << item.getConsts().getDeltaY() << " " 
                << item.getConsts().getDeltaZ() << " " 
                << item.getConsts().getRotX() << " " 
                << item.getConsts().getRotY() << " " 
                << item.getConsts().getRotZ() << " " 
                << endreq;
        }
    }
    log << endreq;

    if (processConstants().isFailure()) {
        log << MSG::ERROR << "Failed to process constants" << endreq;
        return StatusCode::FAILURE;
    }
  
    ConstsVec* ptrConsts = (m_calibType==SIM ? &m_simConsts : &m_recConsts );

    log << MSG::DEBUG ;
    if (log.isActive()) {
        log << "Debug output for alignment consts follows:" << endreq;
        AlignmentConsts alConsts;
        int tower, layer, view, ladder, wafer;
        bool first = true;
        for (tower = 0; tower < m_nTowers; ++tower) {
            for (layer=0; layer<m_nLayers; ++layer) {
                if(first) {
                    log << m_mode 
                        << " consts for tower " << tower
                        << ", layer " << layer;
                }
                for(view=0;view<2;++view) {
                    if(first) log << endreq << "view " << view ;
                    for(ladder=0;ladder<m_nLadders;++ladder) {
                        if(first) log << endreq << " ladder " << ladder << endreq;
                        for(wafer=0;wafer<m_nWafers;++wafer) {
                            if(first) log << "  wafer " << wafer << ": ";
                            alConsts = 
                                (*ptrConsts)(tower, layer, view, ladder, wafer);
                            if ( !alConsts.isNull() ) {
                                log.stream() << alConsts;
                                log << endreq;
                            }
                        }
                    }
                }
                first = false;
            }
        }
    }
    log << endreq;

    if(m_printConsts) {
        std::cout << "text file, m_calibType = " << m_calibType << std::endl;
        int nElements = ptrConsts->getNumElements();
        int i;
        for (i=0; i<nElements;++i) {
            std::cout << (*ptrConsts)[i] << std::endl;
        }
    }
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

    if (m_testMode&(1<<SIM_SHIFT)) {
        m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);
        log << " in Sim" ;
        m_simConsts.setValue(alConsts);
    }
    if (m_testMode&(1<<REC_SHIFT)) {
        m_fileFlag = m_fileFlag|(1<<REC_SHIFT);
        log << " Rec";
        m_recConsts.setValue(alConsts);
    }
    log << endreq;

    return sc;
}    

StatusCode TkrAlignmentSvc::readFromFile()
{    
    // Purpose: fills the m_alignCol vector with info from the alignment input file
    // Inputs:  
    // Outputs: 
    // Dependencies: None
    // Caveats: None

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    std::string mystring, mystring0;
    std::string flag;

    // clear the list
    m_itemCol.clear();
    m_pItem = m_itemCol.begin();
    int lineNum = 0;

    while(!m_dataFile->eof()) {
        mystring = "";
        getline(*m_dataFile, mystring);
        mystring0 = mystring;
        lineNum++;
        // upper-case the string
        unsigned int i;
        for (i=0;i<mystring.size();++i) {mystring[i] = toupper(mystring[i]);}
        int pos = mystring.find("//");
        if(pos>0) mystring = mystring.substr(0,pos);

        // this stuff apparently no longer needed
        //#ifdef DEFECT_NO_STRINGSTREAM
        //std::istrstream mystream(mystring.c_str());
        //#else
        std::istringstream mystream(mystring);
        //#endif

        //count the tokens
        int tokenCount = 0;
        while(mystream>>mystring) {tokenCount++;}
        if (tokenCount==0) continue;
        // now read back the stream
        mystream.clear();  // to reset the eof flag
        mystream.seekg(0); // to get the pointer back to the beginning

        if(!(mystream>>flag)) break ;

        if(flag.substr(0,2)== "//") { // comment line, just skip 
        } else {
            int iflag;
            int type = -1;
            for (iflag=0; iflag<ntypes; ++iflag) {
                if (flag==itemType[iflag]) {
                    type = iflag;
                    break;
                }
            }
            if (type==-1) {
                std::cout << "Error in type " << flag << std::endl;
                return StatusCode::FAILURE;
            }

            int number;
            if (!(mystream>>number)) break;

            AlignmentConsts consts; // created as null
            if (tokenCount==8) {
                double a,b,c,d,e,f;       
                mystream >> a >> b >> c >> d >> e >> f ;
                consts = AlignmentConsts(
                    0.001*a, 0.001*b, 0.001*c,
                    0.001*d, 0.001*e, 0.001*f
                    );
            } else if(tokenCount!=2) {
                log << MSG::ERROR << "Error in Alignment File: " << endreq 
                    << "Line No. " << lineNum << ": " << mystring0 << endreq
                    << "Please fix this and resubmit the job"
                    << endreq;
                return StatusCode::FAILURE;
            }

            AlignmentItem* pItem = new AlignmentItem(iflag, number,consts); 
            m_itemCol.push_back(pItem);
        }
    }
    return sc;
}

bool TkrAlignmentSvc::getNextItem(aType type, AlignmentItem& item) 
{
    // Purpose: checks on the next item in the alignCol
    // Inputs:  requested type, item to be filled
    // Outputs: true if the item matches the input type, and if so returns item
    // Dependencies: None
    // Caveats: None

    if (m_pItem==m_itemCol.end()) return false;
    AlignmentItem* pItem = *m_pItem;
    if (pItem->getType()!=type) return false;
    ++m_pItem;
    item = *pItem;
    return true;
}

StatusCode TkrAlignmentSvc::processConstants()
{    
    // Purpose: produces the basic alignment constants
    // Inputs:  takes input from alignCol
    // Outputs: fills m_simConsts and/or m_recConsts
    // Dependencies: None
    // Caveats: None

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    AlignmentItem item;
    m_pItem = m_itemCol.begin();

    // for each tower in the list, the chain of fillers is called:
    // tower calls tray, tray calls face, etc.

    // do only the ones that are in the input file
    // the rest will automatically have null constants

    while (getNextItem(TOWER,item)) {
        m_tower = item.getNumber();
        if (!m_tkrGeom->isTower(m_tower)) {
            log << MSG::ERROR << "Tower " << m_tower 
                << " doesn't exist. Check input." << endreq;
        }
        m_towerConsts = item.getConsts();
        if (fillTrayConsts().isFailure()) {
            log << MSG::ERROR << "fillTrayConsts failed!" << endreq;
            return StatusCode::FAILURE;
        } 
    }
    return sc;
}

StatusCode TkrAlignmentSvc::fillTrayConsts()
{
    // Purpose: calculates the local tray constants and passes them to the faces
    // Inputs:  takes input from alignCol
    // Outputs: constants passed to the faces
    // Dependencies: None
    // Caveats: None

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // here we have to do all of them
    // even if a tray is not in the list, the tower constants must be passed down.

    // keep track of which ones have been explicitly called
    std::vector<bool> done(m_nTrays, false);

    AlignmentItem item;

    // do the ones that are in the input file
    while (getNextItem(TRAY,item)) {
        m_tray = item.getNumber();
        if (m_tray>=m_nTrays) {
            log << MSG::ERROR << "Tray " << m_tray
                << " doesn't exist. Check input. " << endreq;
            return StatusCode::FAILURE;
        }
        AlignmentConsts thisTray = item.getConsts();
        calculateTrayConsts(thisTray);
        if(fillFaceConsts().isFailure()) {return StatusCode::FAILURE;}
        done[m_tray] = true;
    }

    // if the tower consts are zero, nothing left to do
    if(m_towerConsts.isNull()) return sc;

    // now do the rest, with null tray constants
    int itry = m_nTrays;
    while (itry--) {
        if (!done[itry]) {
            m_tray = itry;
            calculateTrayConsts(nullConsts);
            if(fillFaceConsts().isFailure()) {return StatusCode::FAILURE;}
        }
    }
    return sc;
}

void TkrAlignmentSvc::calculateTrayConsts( AlignmentConsts& thisTray) const
{
    // Purpose: merges the tower and tray constants
    // Inputs:  tray constants
    // Outputs: merged constants
    // Dependencies: None
    // Caveats: None

    double deltaX = thisTray.getDeltaX() + m_towerConsts.getDeltaX()
        + m_towerConsts.getRotY()*m_trayZ[m_tray];
    double deltaY = thisTray.getDeltaY() + m_towerConsts.getDeltaY()
        - m_towerConsts.getRotX()*m_trayZ[m_tray];
    double deltaZ = thisTray.getDeltaZ() + m_towerConsts.getDeltaZ();
    double rotX   = thisTray.getRotX()   + m_towerConsts.getRotX();
    double rotY   = thisTray.getRotY()   + m_towerConsts.getRotY();
    double rotZ   = thisTray.getRotZ()   + m_towerConsts.getRotZ();
    m_trayConsts = AlignmentConsts(deltaX, deltaY, deltaZ, rotX, rotY, rotZ);
    return;
}

StatusCode TkrAlignmentSvc::fillFaceConsts()
{
    // Purpose: calculates the local face constants and passes them to the ladders
    // Inputs:  takes input from alignCol
    // Outputs: constants passed to the ladders
    // Dependencies: None
    // Caveats: None

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // here we have to do all of them, because the tower consts are non-zero

    std::vector<bool> done(m_nFaces, false);

    AlignmentItem item;
    // do the ones that are in the input file
    while (getNextItem(FACE,item)) {
        m_face = item.getNumber();
        if (m_face>=m_nFaces) {
            log << MSG::ERROR << "Face " << m_face
                << " doesn't exist. Check input." << endreq;
            return StatusCode::FAILURE;
        }

        AlignmentConsts thisFace = item.getConsts();
        calculateFaceConsts(thisFace);
        if(fillLadderConsts().isFailure()) {return StatusCode::FAILURE;}
        done[m_face] = true;
    }

    if (m_trayConsts.isNull()) return sc;

    // now do the rest
    int itry = m_nFaces;
    while (itry--) {
        if (!done[itry]) {
            m_face = itry;
            if(m_tkrGeom->trayToPlane(m_tray, m_face)<0) continue;
            calculateFaceConsts(nullConsts);
            if(fillLadderConsts().isFailure()) {return StatusCode::FAILURE;}
        }
    }
    return sc;
}

void TkrAlignmentSvc::calculateFaceConsts( AlignmentConsts& thisFace) const
{
    // Purpose: merges the tray and face constants
    // Inputs:  face constants
    // Outputs: merged constants
    // Dependencies: None
    // Caveats: None

    double zPlane = m_faceZ(m_tray, m_face);

    double deltaX = thisFace.getDeltaX() + m_trayConsts.getDeltaX()
        + m_trayConsts.getRotY()*zPlane;
    double deltaY = thisFace.getDeltaY() + m_trayConsts.getDeltaY()
        - m_trayConsts.getRotX()*zPlane;
    double deltaZ = thisFace.getDeltaZ() + m_trayConsts.getDeltaZ();
    double rotX   = thisFace.getRotX()   + m_trayConsts.getRotX();
    double rotY   = thisFace.getRotY()   + m_trayConsts.getRotY();
    double rotZ   = thisFace.getRotZ()   + m_trayConsts.getRotZ();

    m_faceConsts = AlignmentConsts(deltaX, deltaY, deltaZ, rotX, rotY, rotZ);
    return;
}


StatusCode TkrAlignmentSvc::fillLadderConsts()
{
    // Purpose: calculates the local ladder constants and passes them to the wafers
    // Inputs:  takes input from alignCol
    // Outputs: constants passed to the wafers
    // Dependencies: None
    // Caveats: None

    StatusCode sc = StatusCode::SUCCESS;    
    MsgStream log(msgSvc(), name());

    std::vector<bool> done(m_nLadders, false);

    AlignmentItem item;
    //int nLadders = m_tkrGeom->nWaferAcross();
    while (getNextItem(LADDER,item)) {
        m_ladder = item.getNumber();
        if (m_ladder>=m_nLadders) {
            log << MSG::ERROR << "Ladder " << m_ladder
                << " doesn't exist. Check input." << endreq;
            return StatusCode::FAILURE;
        }

        AlignmentConsts thisLadder = item.getConsts();
        calculateLadderConsts(thisLadder);
        if(fillWaferConsts().isFailure()) {return StatusCode::FAILURE;}
        done[m_ladder] = true;
    }

    if (m_faceConsts.isNull()) return sc;

    // now do the rest
    int itry = m_nLadders;
    while (itry--) {
        if (!done[itry]) {
            m_ladder = itry;
            calculateLadderConsts(nullConsts);
            if(fillWaferConsts().isFailure()) {return StatusCode::FAILURE;}
        }
    }
    return sc;
}

void TkrAlignmentSvc::calculateLadderConsts(AlignmentConsts& thisLadder) const
{
    // Purpose: merges the tray and ladder constants
    // Inputs:  ladder constants
    // Outputs: merged constants
    // Dependencies: None
    // Caveats: None

    int layer, view;
    m_tkrGeom->trayToLayer(m_tray, m_face, layer, view);

    double trayWidth   = m_tkrGeom->trayWidth();
    double ladderGap   = m_tkrGeom->ladderGap();
    // why don't we ask TkrGeometrySvc to calculate this?
    double ladderWidth = (trayWidth - (m_nLadders-1)*ladderGap)/m_nLadders;
    double offset      = -0.5*trayWidth + 0.5*ladderWidth;
    offset += m_ladder*(ladderWidth + ladderGap);

    double xPlane = offset;
    double yPlane = 0.0;
    if (view!=0) std::swap(xPlane, yPlane);

    double deltaX = thisLadder.getDeltaX() + m_faceConsts.getDeltaX()
        - m_faceConsts.getRotZ()*yPlane;
    double deltaY = thisLadder.getDeltaY() + m_faceConsts.getDeltaY()
        + m_faceConsts.getRotZ()*xPlane;
    double deltaZ = thisLadder.getDeltaZ() + m_faceConsts.getDeltaZ()
        +m_faceConsts.getRotX()*yPlane - m_faceConsts.getRotY()*xPlane;
    double rotX = thisLadder.getRotX() + m_faceConsts.getRotX();
    double rotY = thisLadder.getRotY() + m_faceConsts.getRotY();
    double rotZ = thisLadder.getRotZ() + m_faceConsts.getRotZ();

    m_ladderConsts = AlignmentConsts(deltaX, deltaY, deltaZ, rotX, rotY, rotZ);
    return;
}

StatusCode TkrAlignmentSvc::fillWaferConsts()
{   
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    std::vector<bool> done(m_nWafers, false);

    int layer, view;
    m_tkrGeom->trayToLayer(m_tray, m_face, layer, view);

    int index;
    AlignmentItem item;
    
    ConstsVec* ptrConsts = (m_calibType==SIM ? &m_simConsts : &m_recConsts );
    while (getNextItem(WAFER,item)) {
        m_wafer = item.getNumber();
        if (m_wafer>=m_nWafers) {
            log << MSG::ERROR << "Wafer " << m_wafer
                << " doesn't exist. Check input." << endreq;
            return StatusCode::FAILURE;
        }
        if( !ptrConsts->checkIndex(m_tower, layer, view, m_ladder, m_wafer)) continue;
        AlignmentConsts thisWafer = item.getConsts();

        calculateWaferConsts(thisWafer);
        // here are the final constants!
        index = ptrConsts->getIndex(m_tower, layer, view, m_ladder, m_wafer);
        (*ptrConsts)[index] = m_waferConsts;
        done[m_wafer] = true;
    }

    if (m_ladderConsts.isNull()) return sc;

    // now do the rest
    int itry = m_nWafers;

    while (itry--) {
        if (!done[itry]) {
            m_wafer = itry;
            if( !ptrConsts->checkIndex(m_tower, layer, view, m_ladder, m_wafer)) continue;
            AlignmentConsts thisWafer = item.getConsts();
            calculateWaferConsts(nullConsts);
            // here are the final constants!
            index = ptrConsts->getIndex(m_tower, layer, view, m_ladder, m_wafer);
            (*ptrConsts)[index] = m_waferConsts;
        }
    }
    return sc;
}

void TkrAlignmentSvc::calculateWaferConsts(AlignmentConsts& thisWafer) const
{
    // Purpose: merges the ladder and wafer constants
    // Inputs:  wafer constants
    // Outputs: merged constants
    // Dependencies: None
    // Caveats: None

    double trayWidth   = m_tkrGeom->trayWidth();
    double waferGap    = m_tkrGeom->ladderInnerGap();
    double waferWidth  = (trayWidth - (m_nWafers-1)*waferGap + .1)/m_nWafers;
    double trayWidth1   = m_nWafers*waferWidth + (m_nWafers-1)*waferGap;
    double offset      = -0.5*trayWidth1 + 0.5*waferWidth;
    offset += m_wafer*(waferWidth+waferGap);

    double xWafer = 0.0;
    double yWafer = offset;
    int layer, view;
    m_tkrGeom->trayToLayer(m_tray, m_face, layer, view);
    if (view==1) std::swap(xWafer, yWafer);

    double deltaX = thisWafer.getDeltaX() + m_ladderConsts.getDeltaX()
        - m_ladderConsts.getRotZ()*yWafer;
    double deltaY = thisWafer.getDeltaY() + m_ladderConsts.getDeltaY()
        + m_ladderConsts.getRotZ()*xWafer;
    double deltaZ = thisWafer.getDeltaZ() + m_ladderConsts.getDeltaZ()
        +m_ladderConsts.getRotX()*yWafer - m_ladderConsts.getRotY()*xWafer;
    double rotX = thisWafer.getRotX() + m_ladderConsts.getRotX();
    double rotY = thisWafer.getRotY() + m_ladderConsts.getRotY();
    double rotZ = thisWafer.getRotZ() + m_ladderConsts.getRotZ();

    // apply scaling
    // which scale?
    // The overall scale can be useful, but the individual scales are tricky
    //   because all the constants get mixed in the final corrections
    //   I'll leave them in, but beware!!!!

    double scale;
    std::vector<double> scaleVec;

    if(m_calibType==SIM) {
        scale = m_simScale;
        scaleVec = m_simScaleVec;
    } else {
        scale = m_recScale;
        scaleVec = m_recScaleVec;
    }

    deltaX *= scale*scaleVec[0];
    deltaY *= scale*scaleVec[1];
    deltaZ *= scale*scaleVec[2];
    rotX   *= scale*scaleVec[3];
    rotY   *= scale*scaleVec[4];
    rotZ   *= scale*scaleVec[5];

    m_waferConsts = AlignmentConsts(deltaX, deltaY, deltaZ,
        rotX, rotY, rotZ);

    return;
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(calibType type, int tower, int layer, 
                                                  int view, int ladder, int wafer) const
{
    // Purpose:  return pointer to an set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set

    ConstsVec* ptrConsts = (type==SIM ? &m_simConsts : &m_recConsts );
    int index = ptrConsts->getIndex(tower, layer, view, ladder, wafer);

    return getConsts(type, index);
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(calibType type, 
                                                  idents::VolumeIdentifier id) const
{
    // Purpose:  return pointer to a set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set

    int tower  = idents::TowerId(id[2],id[1]).id();
    int tray   = id[4];
    int botTop = id[6];
    int layer, view;
    m_tkrGeom->trayToLayer(tray, botTop, layer, view);
    view       = id[5];
    int ladder = id[7];
    int wafer  = id[8];

    ConstsVec* ptrConsts = (type==SIM ? &m_simConsts : &m_recConsts );
    int index = ptrConsts->getIndex(tower, layer, view, ladder, wafer);

    return getConsts(type, index);
}

const AlignmentConsts* TkrAlignmentSvc::getConsts(calibType type, int index) const
{
    // Purpose:  return pointer
    // Inputs:   index
    // Outputs:  pointer to that set of consts

    ConstsVec* ptrConsts = (type==SIM ? &m_simConsts : &m_recConsts );

    if (index>=0 && index < ptrConsts->getNumElements()) {
        return &(*ptrConsts)[index];
    } else {
        AlignmentConsts* result = 0;
        return result;
    }
}

void TkrAlignmentSvc::moveMCHit(idents::VolumeIdentifier id, HepPoint3D& entry,
                                HepPoint3D& exit,
                                HepVector3D dir) const
{
    // Purpose:     Move an McHit according to alignment constants
    // Inputs:      volId and entry and exit points
    // Output:      modified entry and exit point arguments

    //This is called "dir" but what I'm calculating is a direction normalized to dir.z() = 1;
    //  This yields the slopes in x and y as the x() and y() components.
    //HepVector3D dir = (exit-entry);

    HepVector3D tempVec;
    if(dir.mag()==0.0) {tempVec = exit - entry;}
    else               {tempVec = dir;}
    double delz = dir.z();

    if (delz!=0) { tempVec = tempVec/delz;}
    else         { tempVec = HepVector3D(0., 0., 1.);}    

    // calculation is done separately for entry and exit, because transformation depends
    //    in 2nd order (~ microns) on the coordinates. This can be speeded up if
    //    necessary, by basing the calculation on the average coordinates

    const AlignmentConsts* alConsts = getConsts(SIM, id);

    int view = id[5];

    HepVector3D deltaEntry = getDelta(view, entry, tempVec, alConsts);
    HepVector3D deltaExit  = getDelta(view, exit,  tempVec, alConsts);

    // for now, limit delta to m_maxDelta (default = 5mm, modifiable in jobOptions
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

HepVector3D TkrAlignmentSvc::deltaReconPoint(const HepPoint3D& point, const HepVector3D& dir, 
                                             int layer, int view, 
                                             unsigned flags,
                                             alignTask task, 
                                             const AlignmentConsts* consts
                                             ) const
{   
    // alignTask can have the values NULLTASK, APPLYCONSTS, pr FINDTOWERCONSTS.
    //   if NULLTASK, a zero vector is returned
    //   if APPLYCONSTS, the constants are supplied by the alignment service
    //   if FINDTOWERCONSTS,  the constants are expected to be passed in
    //      if they're not, a zero vector is returned;
    //      In this mode, transformations are made with respect to local tower 
    //         coordinates, not wafer coordinates.

    HepPoint3D localPoint;
    int nXTower, nYTower;
    idents::VolumeIdentifier volId;

    if ((task==NULLTASK) || (task!=APPLYCONSTS && consts==0)) return HepVector3D(0., 0., 0.);

    if (task==APPLYCONSTS) {
        volId = getGeometryInfo(layer, view, point, localPoint);
    } else {
        localPoint = getTowerCoordinates(point, nXTower, nYTower);
    }

    double small = 1.e-3;
    double dirZ   = dir.z();
    if (fabs(dirZ)<small) dirZ = small;
    double alphaX = dir.x()/dirZ;
    double alphaY = dir.y()/dirZ;

    double deltaPointX, deltaPointY;

    m_faceConsts = nullConsts;
    if(task==FINDTOWERCONSTS) {
        // need to generate the "face" constants to apply to the "face" coordinates
        // first get tray and view
        int tray, face;
        m_tkrGeom->layerToTray(layer, view, tray, face);
        m_tray = tray;
        m_face = face;
        // now the input consts
        m_towerConsts = *consts;
        // this is for testing
        //m_towerConsts = AlignmentConsts(0.0, 0.0, 0.0, 0.0, 0.002, 0.0);
        calculateTrayConsts(nullConsts);
        calculateFaceConsts(nullConsts);
        //std::cout << "tower coord: " << localPoint << std::endl;
        //std::cout << "calculated face consts: " << m_faceConsts << std::endl;
    }

    const AlignmentConsts* alConsts = (task==APPLYCONSTS ? getConsts(REC, volId): &m_faceConsts);

    applyDelta(localPoint.x(), localPoint.y(), alphaX, alphaY, alConsts,
        deltaPointX, deltaPointY, flags);

    // for now, limit delta to m_maxDelta (default is 5 mm, modifiable in jobOptions
    // later fix transformation for special cases
    double mag;
    HepVector3D dirDelta;
    HepVector3D deltaPoint = HepVector3D(deltaPointX, deltaPointY, 0.0);

    mag = std::min(m_maxDelta, deltaPoint.mag());
    dirDelta = deltaPoint.unit();
    deltaPoint = mag*dirDelta;

    return -deltaPoint;
}

//void TkrAlignmentSvc::moveReconPoint(HepPoint3D& point, const HepVector3D& dir, 
//                                     int layer, int view, alignTask task, 
//                                     const AlignmentConsts* consts,
//                                     const unsigned flags) const
//{
//    HepVector3D deltaPoint = deltaReconPoint(point, dir, layer, view, task, consts, flags);
//
//    //now subtract(??) this delta from the global point
//    point += deltaPoint;
//    return;
//}

HepPoint3D TkrAlignmentSvc::getTowerCoordinates(
    const HepPoint3D& globalPoint, int& nXTower, int& nYTower) const
{
    double xTower = m_tkrGeom->truncateCoord(globalPoint.x(), m_tkrGeom->towerPitch(), 
        m_tkrGeom->numXTowers(), nXTower);
    double yTower = m_tkrGeom->truncateCoord(globalPoint.y(), m_tkrGeom->towerPitch(), 
        m_tkrGeom->numYTowers(), nYTower);
    return HepPoint3D(xTower, yTower, globalPoint.z());
}

idents::VolumeIdentifier TkrAlignmentSvc::getGeometryInfo(
    int layer, int view, const HepPoint3D& globalPoint, HepPoint3D& alignmentPoint) const
{
    int nXTower, nYTower, ladder, wafer;

    HepPoint3D towerPoint
        = getTowerCoordinates(globalPoint, nXTower, nYTower);

    double ladderPitch = m_tkrGeom->ladderPitch();
    double waferPitch  = m_tkrGeom->waferPitch();
    int nLadders = m_tkrGeom->nWaferAcross();
    double xLocal, yLocal;

    // local means in the coordinate system of the wafer...
    //  for y-measuring trays, this is rotated 90 degrees
    // 
    // because the wafers are *numbered backwards* for the y-measuring planes,
    //  the transformation is not obvious. The minus signs that you would expect
    //  are *not* there.

    double xTower = towerPoint.x();
    double yTower = towerPoint.y();

    if (view==1) std::swap(xTower, yTower);
    xLocal = m_tkrGeom->truncateCoord(xTower, ladderPitch, nLadders, ladder);
    yLocal = m_tkrGeom->truncateCoord(yTower, waferPitch,  nLadders, wafer);

    int tray, botTop;
    m_tkrGeom->layerToTray(layer, view, tray, botTop);
    idents::VolumeIdentifier id;
    id.init(0,0);
    id.append(0);
    id.append(nYTower);
    id.append(nXTower);
    id.append(1);
    id.append(tray);
    id.append(view);
    id.append(botTop);
    id.append(ladder);
    id.append(wafer);

    if(view==1) std::swap(xLocal, yLocal);
    alignmentPoint = HepPoint3D(xLocal, yLocal, 0.);

    return id;
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
        pointX = -point2;
        pointY = point1;
        alphaX = -alpha2;
        alphaY = alpha1;
    }

    double deltaPointX, deltaPointY;
    applyDelta(pointX, pointY, alphaX, alphaY, alConsts, deltaPointX, deltaPointY);

    if (view==0) {
        return HepVector3D(deltaPointX,  deltaPointY, 0.);
    } else {
        return HepVector3D(deltaPointY, -deltaPointX, 0.);
    }
}

void TkrAlignmentSvc::applyDelta(double pointX, double pointY, 
                                 double alphaX, double alphaY,
                                 const AlignmentConsts* alConsts, 
                                 double& deltaPointX, double& deltaPointY,
                                 unsigned int flags) const
{
    
  /*
    double deltaX = ((flags&XANDY)!=0 || !m_useFlags ? alConsts->getDeltaX() : 0.0);
    double deltaY = ((flags&XANDY)!=0 || !m_useFlags ? alConsts->getDeltaY() : 0.0);
    double deltaZ = ((flags&ANGLE)!=0 || !m_useFlags ? alConsts->getDeltaZ() : 0.0);

    double rotX =   ((flags&ANGLE)!=0 || !m_useFlags ? alConsts->getRotX() : 0.0);
    double rotY =   ((flags&ANGLE)!=0 || !m_useFlags ? alConsts->getRotY() : 0.0);
    double rotZ =   ((flags&ROTZ)!=0  || !m_useFlags ? alConsts->getRotZ() : 0.0);
  */

    double deltaX = alConsts->getDeltaX();
    double deltaY = alConsts->getDeltaY();
    double deltaZ = alConsts->getDeltaZ();

    double rotX =   alConsts->getRotX();
    double rotY =   alConsts->getRotY();
    double rotZ =   alConsts->getRotZ();

	//std::cout << "TkrAlignmentSvc::applyDelta, dx/y/z/ drotx/y/z "
    //    << deltaX << " " << deltaY << " " << deltaZ << " " << rotX << " " << rotY << " " << rotZ << std::endl;


    double rotTerm = deltaZ + rotX*pointY - rotY*pointX; 

    deltaPointX = - deltaX + rotZ*pointY
        + alphaX*rotTerm;
    deltaPointY = - deltaY - rotZ*pointX
        + alphaY*rotTerm;

    //std::cout << "   deltaPointX/Y " << deltaPointX << " " << deltaPointY << std::endl;
}

void TkrAlignmentSvc::update(CalibData::TkrTowerAlignCalib* pTowerAlign, 
                             CalibData::TkrInternalAlignCalib* pInternalAlign) 
{
    // This fills the alignment arrays from the appropriate xml file
    // It is called every time the event falls outside the current validity range
    // At first I thought it would be really clever to modify the existing code
    // so that the same structure would work for both .txt and .xml files.
    // But then I decided that it wouldn't be that clever after all.

    ConstsVec* ptrConsts = (m_calibType==SIM ? &m_simConsts : &m_recConsts );

    MsgStream log(msgSvc(), name());

    int layer, view, index;
    int nElements = ptrConsts->getNumElements();

    //clear out the consts vector
    for(index=0; index<nElements; ++index) {
        (*ptrConsts)[index] = nullConsts;
    }

    StatusCode sc = StatusCode::SUCCESS;

    CLHEP::Hep3Vector disp, rot;
    AlignmentConsts thisTray, thisFace, thisLadder, thisWafer;
    for(m_tower=0; m_tower<m_nTowers; ++m_tower) {
        // get the tower constant
        m_towerConsts = nullConsts;
        if(pTowerAlign) {
            sc = pTowerAlign->getTowerAlign(m_tower, disp, rot);
            m_towerConsts = makeConsts(disp, rot);
            log << MSG::INFO << "inter-tower constants for calib type " << m_calibType << " loaded" << endreq;
        }
        else
            log << MSG::WARNING << "inter-tower constants for calib type " << m_calibType << " not loaded" << endreq;
        // log internal align only once!
        if ( pInternalAlign )
            log << MSG::INFO << "intra-tower constants for calib type " << m_calibType << " loaded" << endreq;
        else
            log << MSG::WARNING << "intra-tower constants for calib type " << m_calibType << " not loaded" << endreq;
        for(m_tray=0; m_tray<m_nTrays; ++m_tray) {
            // then the tray
            thisTray = nullConsts;
            if(pInternalAlign) {
                sc = pInternalAlign->getTrayAlign(m_tower, m_tray, disp, rot);
                thisTray = makeConsts(disp, rot);
            }
            calculateTrayConsts(thisTray);
            for(m_face=0;m_face<m_nFaces;++m_face) {
                // then face, if it exists
                m_tkrGeom->trayToLayer(m_tray, m_face, layer, view);
                if (layer<0) continue;
                thisFace = nullConsts;
                if(pInternalAlign) {
                    sc = pInternalAlign->getFaceAlign(
                        m_tower, m_tray, m_face, disp, rot);
                    thisFace = makeConsts(disp, rot);
                }
                calculateFaceConsts(thisFace);
                for(m_ladder=0;m_ladder<m_nLadders;++m_ladder) {
                    // next, ladder
                    thisLadder = nullConsts;
                    if(pInternalAlign) {
                        sc = pInternalAlign->getLadderAlign(
                            m_tower, m_tray, m_face, m_ladder, disp, rot);
                        thisLadder = makeConsts(disp, rot);
                    }
                    calculateLadderConsts(thisLadder);
                    for(m_wafer=0;m_wafer<m_nWafers;++m_wafer) {
                        // finally, wafer
                        thisWafer = nullConsts;
                        if(pInternalAlign) {
                            sc = pInternalAlign->getWaferAlign(
                                m_tower, m_tray, m_face, m_ladder, m_wafer, disp, rot);
                            thisWafer = makeConsts(disp, rot);
                        }
                        calculateWaferConsts(thisLadder);
                        // store it!
                        int index = ptrConsts->getIndex(
                            m_tower, layer, view, m_ladder, m_wafer);
                        (*ptrConsts)[index] = m_waferConsts;
                    }
                }
            }
        }
    }

    m_fileFlag = m_fileFlag | (1 << (m_calibType==SIM? SIM_SHIFT : REC_SHIFT));

    // print out the consts vector
    if(m_printConsts) {
        std::cout << "calibType = " << m_calibType << std::endl;
        for(index=0; index<nElements; ++index) {
            std::cout << (*ptrConsts)[index] << std::endl;
        }
    }
   
    if(sc.isFailure()||!sc.isFailure()) return;
}

IGeometry::VisitorRet TkrAlignmentGeomVisitor::pushShape(ShapeType /* s */, const UintVector& idvec, 
                                                         std::string name, std::string /* material*/, const DoubleVector& params, 
                                                         VolumeType /*type*/, SenseType /*sense*/)
{
    bool debug = false;

    if(name=="oneCAL") {
        return AbortSubtree;
    }

    if (name=="towerRow") {
        // some of the names are blank if you don't start from topVol = LAT
        if(idvec.size()==0) return More;
        if(idvec[0]>0) return AbortSubtree;
    }
    if (name=="oneTower") { 
        if (debug) {
            std::cout << "Element: " << name << std::endl;
            std::cout << " idvec("<< idvec.size() << " elements):" ;
            for (unsigned int iv=0; iv<idvec.size(); ++iv ) {
                std::cout << idvec[iv] << " ";
            }
            std::cout << std::endl;
        }
        if (idvec.size()==0) return More;
        if(idvec[0]>0) return AbortSubtree;
    }

    // trayBot or emTrayBot
    if (name.find("rayBot")!=std::string::npos) {m_trayBotHeight = params[8];}

    // trayTop or emTrayTop
    if (name.find("rayTop")!=std::string::npos) {m_trayTopHeight = params[8];}

    if (debug) {
        // something like trayBot or emTraySuper or oneTkrStack
        if (name.find("ray")!=std::string::npos || name.find("neTKR")!=std::string::npos) { 
            std::cout << name << ": z, height " 
                << params[2] << ", " << params[8] << std::endl;
        }
    }
    if (name.find("rayTop")!=std::string::npos) return AbortSubtree;

    return More;
}

void TkrAlignmentGeomVisitor::popShape() { return;}


// queryInterface

StatusCode  TkrAlignmentSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
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

const InterfaceID&  TkrAlignmentSvc::type () const {
    return IID_ITkrAlignmentSvc;
}
