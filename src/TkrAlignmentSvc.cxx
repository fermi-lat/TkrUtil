
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "src/TkrAlignmentSvc.h"

#include "idents/TowerId.h"

#include <fstream>
#include <algorithm>
#include <iostream>

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
    declareProperty("simFile", m_simFile);
    declareProperty("recFile", m_recFile);
    declareProperty("testMode", m_testMode);
    
    return;
}

StatusCode TkrAlignmentSvc::initialize()
{
    // Purpose: reads in (or set) alignment constants
    // Inputs:  None
    // Outputs: Status code (Success/Failure)
    
    MsgStream log(msgSvc(), name());
    log.setLevel(MSG::DEBUG);
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    
    m_simFile = "";
    m_recFile = "";
    m_testMode = 0;
    m_fileFlag = 0;
    
    setProperties();
    
    // initialize files
    
    int i;
    AlignmentConsts alConsts(0., 0., 0., 0., 0., 0.);
    for (i=0;i<NELEMENTS; i++) {
        m_simConsts[i] = alConsts;
        m_recConsts[i] = alConsts;
    }
    
    if(!m_testMode) {
        
        // If there is no bad strips file, service will do nothing
        if (m_simFile=="" && m_recFile=="") {        
            log << MSG::INFO << "No alignment was requested." << endreq;
            log << MSG::INFO << "  No alignment will be done." << endreq;
            return sc;
        }
        
        // this method resolves environmental variables in the file name
        if( m_simFile != "") {
            xml::IFile::extractEnvVar(&m_simFile);    
            log << MSG::INFO << "Input file for simulation alignment: " 
                << m_simFile << endreq;
        }
        
        if( m_recFile != "") {
            xml::IFile::extractEnvVar(&m_recFile);    
            log << MSG::INFO << "Input file for reconstruction alignment: " 
                << m_recFile << endreq;
        }
        
        // read in the sim file, if requested
        if (m_simFile!="") {
            std::ifstream simFile;
            simFile.open( m_simFile.c_str());
            
            if (!simFile) {
                log << MSG::ERROR << "  Simulation Alignment file not found: check jobOptions." << endreq;
                return StatusCode::FAILURE;
            }
            log << MSG::DEBUG << "simulation alignment constants:" << endreq;
            readFromFile(&simFile, m_simConsts);
            simFile.close();
            m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);
        }
        
        // read in the rec file, if requested.
        if (m_recFile!="") {
            std::ifstream recFile;
            recFile.open( m_recFile.c_str());
            
            if (!recFile) {
                log << MSG::ERROR << "  Simulation Alignment file not found: check jobOptions." << endreq;
                return StatusCode::FAILURE;
            }
            log << MSG::DEBUG << "reconstruction alignment constants:" << endreq;
            readFromFile(&recFile, m_recConsts);
            recFile.close();
            m_fileFlag = m_fileFlag|(1<<REC_SHIFT);
            
        }
    } else {
        log << MSG::INFO << "Test mode " << m_testMode << 
            " requested, each hit will be translated by +1 strip in Global X" << endreq
            << "               and +2 strips in Global Y";
        AlignmentConsts alConsts(0.22, 0.44, 0., 0., 0., 0.);
        if (m_testMode&(1<<SIM_SHIFT)) {
            m_fileFlag = m_fileFlag|(1<<SIM_SHIFT);
            log << " in Sim" ;
            for (i=0; i<NELEMENTS; i++) {m_simConsts[i] = alConsts;
            }
        }
        if (m_testMode&(1<<REC_SHIFT)) {
            m_fileFlag = m_fileFlag|(1<<REC_SHIFT);
            log << " Rec";
            for (i=0; i<NELEMENTS; i++) {m_recConsts[i] = alConsts;
            }
        }
        log << endreq;
    }
       
    return sc;
}

StatusCode TkrAlignmentSvc::finalize()
{
    return StatusCode::SUCCESS;
}

void TkrAlignmentSvc::readFromFile(std::ifstream* file, AlignmentConsts* aConsts)
{    
    // Purpose:
    // Inputs:  File name
    // Outputs: None
    // Dependencies: None
    // Caveats: None
    
    bool read = true;           // for testing
    
    std::string junk;
    
    // format of file:
    
    while(read && !file->eof()) {
        int tower, layer, view, ladder, wafer;
        int sequence;
        
        
        *file >> tower ;
        if(tower==-1) { // comment line, just skip 
            std::getline(*file, junk);
            continue;
        }
        
        if (file->eof()) break;
        *file >> layer >> view >> ladder >> wafer;
        
        sequence = getIndex(tower, layer, view, ladder, wafer);
        
        std::cout << "Alignment Consts[ "<< sequence << "] for t,l,v,l,w " << tower << " " 
            << layer << " " << view << " " << ladder << " " << wafer << std::endl;
        
        double a,b,c,d,e,f;
        
        *file >> a >> b >> c >> d >> e >> f ;
        
        std::cout << "   alConsts: " << a << " " << b << " " << c << " " 
            << d << " " << e << " " << f << std::endl;
        
        if(sequence>=0 && sequence<NELEMENTS) {
            aConsts[sequence] = AlignmentConsts(a,b,c,d,e,f);
        }       
    }  
    return;
}

int TkrAlignmentSvc::getIndex(int tower, int layer, 
                              int view, int ladder, int wafer) const
{
    // Purpose:  calculate index into array of consts
    // Inputs:   tower, bilayer, view, ladder, wafer
    // Outputs:  index
    
    int index = -1;
    if (layer<0 || layer>=NLAYERS || tower<0 || tower>=NTOWERS ||
        view<0 || view>=NVIEWS || ladder<0 || ladder >= NLADDERS || wafer<0 || wafer>=NWAFERS)
    {return index;}
    // for now, hardwired to be as large as will ever by needed
    return wafer + NWAFERS*(ladder + NLADDERS*(view + NVIEWS*(layer + NLAYERS*tower)));
}

AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, int tower, int layer, 
                                            int view, int ladder, int wafer) const
{
    // Purpose:  return pointer to an set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set
    
    int index = getIndex(tower, layer, view, ladder, wafer);
    
    return getConsts(type, index);
}

AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, idents::VolumeIdentifier id) const
{
    // Purpose:  return pointer to a set of alignment consts
    // Inputs:   const type, tower, layer, view, ladder, wafer
    // Outputs:  pointer to that set
    
    int tower  = idents::TowerId(id[2],id[1]).id();
    int tray   = id[4];
    int botTop = id[6];
    int layer  = tray+botTop-1;
    int view   = id[5];
    int ladder = id[7];
    int wafer  = id[8];

    int index = getIndex(tower, layer, view, ladder, wafer);
    
    return getConsts(type, index);
}

AlignmentConsts* TkrAlignmentSvc::getConsts(constType type, int index) const
{
    // Purpose:  return pointer
    // Inputs:   index
    // Outputs:  pointer to that set of consts
       
    if (index>=0 && index < NELEMENTS) {
        
        if (type==SIM) {
            return const_cast<AlignmentConsts*>(&m_simConsts[index]);
        } else {
            return const_cast<AlignmentConsts*>(&m_recConsts[index]);
        }
        
    } else {
        AlignmentConsts* result = 0;
        return 0;
    }
}

void TkrAlignmentSvc::moveMCHit(idents::VolumeIdentifier id, HepPoint3D& entry,
                                HepPoint3D& exit) const
{
    // Purpose:     Move an McHit according to alignment constants
    // Inputs:      volId and entry and exit points
    // Output:      modified entry and exit point arguments
    
    HepVector3D dir = (exit-entry);
    double delz = dir.z();
    
    dir = (delz != 0) ? dir/delz : HepVector3D( 0., 0., 1.);
    
    // calculation is done separately for entry and exit, because transformation depends
    //    in 2nd order (~ microns) on the coordinates. This can be speeded up if
    //    necessary, by basing the calculation on the average coordinates
    
    AlignmentConsts* alConsts = getConsts(SIM, id);

    int view = id[5];
    
    HepVector3D deltaEntry = getDelta(view, entry, dir, alConsts);
    HepVector3D deltaExit  = getDelta(view, exit,  dir, alConsts);
    
    entry = entry + deltaEntry;
    exit  = exit  + deltaExit;
    
    return;
}

void TkrAlignmentSvc::moveCluster(int tower, int layer, int view, int ladder,
                                  HepPoint3D& point) const
{
    // Purpose:  Move an cluster according to alignment constants (dx and dy only)
    // Inputs:   tower, layer, view, ladder, and cluster position
    // Output:   modified position argument
    
    int wafer = 0;
    AlignmentConsts* alConsts = getConsts(REC, tower, layer, view, ladder, wafer);
    // "x" is the only thing we can correct at this stage...
    //     average for all the wafers in a ladder would be better.
    AlignmentConsts alConsts1(alConsts->getDeltaX(), alConsts->getDeltaY()); 
    HepVector3D dir(0., 0., 1);

    view = -1;

    point = point - getDelta(view, point, dir, &alConsts1);
}

HepVector3D TkrAlignmentSvc::getDelta(int view, const HepPoint3D& point,
                                      const HepVector3D& dir,
                                      const AlignmentConsts* alConsts) const 
{
    // Purpose:  Calculate translation in x,y for fixed z
    // Inputs:   view, point and direction
    // Output:   vector of translation
    
    
    double pointX = point.x(); 
    double pointY = point.y(); 
    double pointZ = point.z();
    
    double deltaX = alConsts->getDeltaX();
    double deltaY = alConsts->getDeltaY();
    double deltaZ = alConsts->getDeltaZ();
    
    double deltaRotX = alConsts->getDeltaRotX();
    double deltaRotY = alConsts->getDeltaRotY();
    double deltaRotZ = alConsts->getDeltaRotZ();
    
    double alphaX = dir.x();
    double alphaY = dir.y();
    
    double rotTerm = deltaZ + deltaRotY*pointX - deltaRotX*pointY; 
    
    double deltaPointXLocal = - deltaX - deltaRotZ*pointY
        + alphaX*rotTerm;
    double deltaPointYLocal = - deltaY + deltaRotZ*pointX
        + alphaY*rotTerm;
    
    if(view==1) {
        return HepVector3D(deltaPointYLocal, -deltaPointXLocal, 0.);
    } else {
        return HepVector3D(deltaPointXLocal,  deltaPointYLocal, 0.);
    }
}


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
