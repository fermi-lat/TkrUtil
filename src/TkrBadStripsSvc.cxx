/** 
 @file TkrBadStripsSvc.cxx

 @brief Maintains lists of bad strips, and provides access methods 

 First version 3-Jun-2001
  @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrBadStripsSvc.cxx,v 1.9 2003/12/03 19:12:26 lsrea Exp $
*/


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "src/TkrBadStripsSvc.h"
#include "idents/TowerId.h"

#include "GaudiKernel/SmartDataPtr.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "facilities/Util.h"

static const SvcFactory<TkrBadStripsSvc> s_factory;
const ISvcFactory& TkrBadStripsSvcFactory = s_factory;


// Service parameters which can be set at run time must be declared.
// This should be done in the constructor.

TkrBadStripsSvc::TkrBadStripsSvc(const std::string& name, 
                                 ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the file to get data from   
    declareProperty("badStripsFile", m_badStripsFile="");
    declareProperty("killDigi", m_killDigi=false );
    m_visitor = 0;
       
    return;
}

StatusCode TkrBadStripsSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if(m_visitor==0) {
        m_visitor = new BadVisitor;
        m_visitor->setService(this);
    }

    Service::initialize();
    
    m_badStripsFile = "";
    m_empty = true;
    
    setProperties();
    
    for (int i=0;i<NELEMENTS;i++) { m_stripsCol[i].clear();}
    
    // commented code in this routine was the original attempt to implement
    //  the bad strips as a vector of vectors.
    
    // int size = 0;
    // makeCol(size); //make sure that Collection is sensibly initialized

    // this method resolves environmental variables in the file name
    if (m_badStripsFile!="") {
        //xml::IFile::extractEnvVar(&m_badStripsFile); 
        int ret =  facilities::Util::expandEnvVar(&m_badStripsFile);
        if (ret>=0) {
            log << MSG::INFO << "Input file for bad strips: " 
                << m_badStripsFile << endreq;
        } else {
            log << MSG::ERROR << "Input filename " << m_badStripsFile << " not resolved" << endreq;
            return StatusCode::FAILURE;
        }
    }

    sc = doInit();
    
    return sc;
}


StatusCode TkrBadStripsSvc::doInit()
{
    // Purpose: reads in bad strips and constructs in-memory bad strip vectors
    // Inputs:  None
    // Outputs: Status code (Success/Failure)
    
    MsgStream log(msgSvc(), name());
    
    StatusCode sc = StatusCode::SUCCESS;
    
    // If there is no bad strips file, service will do nothing
    if (m_badStripsFile=="") {        
        log << MSG::INFO << "No bad strips file was requested." << endreq;
        return sc;
    }
        
    // open bad strips file
    std::ifstream file;
    file.open( m_badStripsFile.c_str());
    
    if (!file) {
        log << MSG::ERROR << "  File not found: check jobOptions." << endreq;
        return StatusCode::FAILURE;
    }
    
    //int size = m_ntowers*m_nlayers*m_nviews;    
    //makeCol(size);
    
    readFromFile(&file);
    
    file.close();
    
    return sc;
}

StatusCode TkrBadStripsSvc::finalize()
{
    delete m_visitor;
    return StatusCode::SUCCESS;
}

/*void TkrBadStripsSvc::makeCol(const int size)
{
//    m_stripsCol.assign(size);
//  return;
}
*/

StatusCode TkrBadStripsSvc::update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot)
{
    MsgStream log(msgSvc(), name());
    m_visitor->setLog(&log);
    log << MSG::INFO << "updater called " << endreq;
    StatusCode sc = doInit();
    pDead->traverse(m_visitor);
    pHot->traverse(m_visitor);
    m_empty = m_empty && m_visitor->isEmpty();
    
    return sc;
    
}

void TkrBadStripsSvc::readFromFile(std::ifstream* file)
{    
    // Purpose: read bad strips from file and make in-memory vectors
    // Inputs:  File name
    // Outputs: None
    // Dependencies: None
    // Caveats: None
    
    bool read = true;           // for testing
    bool makestrips = true;     // for testing
    
    int nStrips = 0;
    std::string junk;
    
    // format of file:
    //
    // -1 at beginning of line is a comment
    //
    // -2 at beginning of line sets the tag (default is 1)
    // 
    // for each layer with bad strips:
    // tower# sequence#  [strip#] [strip#] ...  -1
    // all whitespace is ignored, so data for a layer may span lines
    //
    // 
    
    while(read && !file->eof()) {
        int tower;
        int sequence;
        
        int tag = 1;
        
        *file >> tower ;
        if(tower==-1) { // comment line, just skip 
            std::getline(*file, junk);
            continue;
        }
        if (tower==-2) { // set the tag
            *file >> tag;
            std::getline(*file, junk);
            continue;
        }
        
        if (file->eof()) break;
        *file >> sequence;
        
        // kludge until the geometry supplies this info
        // converts layer (0...35) to bilayer (digi format) and view
        int layer = sequence/2;
        int element = (sequence+3)%4;
        int view = element/2;
        
        stripCol* v;
        // my private use of getBadStrips requires non-const pointer
        // to build the vector of bad strips...
        // but public uses should return const pointer, so...
        
        if (makestrips) v = const_cast<stripCol*> 
            (getBadStrips(tower, layer, 
            static_cast<idents::GlastAxis::axis>(view)));
        int strip = -1;
        *file >> strip;
        while (strip>=0) {
            if (makestrips) {
                addStrip(v, TaggedStrip(strip, tag));
                *file >> strip;
                nStrips++;
            }
            
            // sort strips in ascending order of strip number 
            // after each line is read in
            if (makestrips) std::sort(v->begin(), v->end());           
        }  
        if (!v->empty()) {m_empty = false;}       
    }
    return;
}

int TkrBadStripsSvc::getIndex(int tower, int layer, 
                              idents::GlastAxis::axis axis) const
{
    // Purpose:  calculate index into array of vectors
    // Inputs:   tower, bilayer, axis
    // Outputs:  index
    
    int view; 
    // this is to decouple the store from the current definition of axes
    // not that it will ever change
    int index = -1;
    if (axis==idents::GlastAxis::X)  {view = 0;}
    else if (axis==idents::GlastAxis::Y) {view = 1;}
    else {return index;}
    if (layer<0 || layer>=NLAYERS || tower<0 || tower>=NTOWERS) {return index;} 
    // for now, hardwired to be as large as will ever by needed
    return view + NVIEWS*(layer + NLAYERS*tower);
}

void TkrBadStripsSvc::addStrip(stripCol* v, TaggedStrip taggedStrip) 
{
    // Purpose: add a bad strip to the list, already tagged bad
    // Inputs:  strip number
    // Outputs: None
    
    v->push_back(taggedStrip);
    return;
}

const stripCol* TkrBadStripsSvc::getBadStrips(int tower, int layer, 
                                              idents::GlastAxis::axis axis) const
{
    // Purpose:  return pointer to a bad strip vector
    // Inputs:   tower, layer, axis
    // Outputs:  pointer to that vector
    
    int index = getIndex(tower, layer, axis);
    
    return getBadStrips(index);
}


const stripCol* TkrBadStripsSvc::getBadStrips(int index) const
{
    // Purpose:  return pointer to a bad strip vector
    // Inputs:   index
    // Outputs:  pointer to that vector
    
    // original code... maybe some day...
    // int ind = (m_stripsCol.size()==0) ? m_stripsCol.size() : index;
    
    if (index>=0 && index < NELEMENTS) {return &m_stripsCol[index];}
    else                         {return 0;}
}


bool TkrBadStripsSvc::isBadStrip(int tower, int layer, 
                                 idents::GlastAxis::axis axis, 
                                 int strip) const 
{
    // Purpose: determine if a given strip is bad
    // Inputs:  tower, bilayer, axis, strip#
    // Output:  true if strip is in the list ( that is, is bad)
    
    const stripCol* v = getBadStrips(tower, layer, axis);
    return isBadStrip(v, strip);
}

bool TkrBadStripsSvc::isBadStrip(const stripCol* v, int strip) const
{
    // Purpose: determine if a given strip is bad
    // Inputs:  index, strip#
    // Output:  true if strip is in the list
    
    bool isBad = false;
    stripCon_it it;
    for (it=v->begin(); it!=v->end(); it++) {
        if ( (*it).getStripNumber()==strip ) {
            isBad = true;
            break;
        }
    }
    return isBad;
    
    // this was the orginal code, before the tag was variable
    // might be useful again some day
    //stripCol_it it = std::find(v->begin(), v->end(), tagBad(strip));
    //return (it!=v->end());
}

bool TkrBadStripsSvc::empty() const
{
    return m_empty;
}

// queryInterface

StatusCode  TkrBadStripsSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrBadStripsSvc == riid) {
        *ppvIF = dynamic_cast<ITkrBadStripsSvc*> (this);
    } else if(IID_ITkrBadStripsSvcCalib == riid) {
        *ppvIF = dynamic_cast<ITkrBadStripsSvcCalib*> (this);
    } else {
        return Service::queryInterface (riid, ppvIF);
    }
    return StatusCode::SUCCESS;
}

// access the type of this service

const IID&  TkrBadStripsSvc::type () const {
    return IID_ITkrBadStripsSvc;
}


CalibData::eVisitorRet BadVisitor::badTower(unsigned int /* row */, 
                                            unsigned int /* col */,
                                            int /* badness */) {
    return CalibData::CONT;    
}

CalibData::eVisitorRet BadVisitor::badPlane(unsigned int row, 
                                            unsigned int col, 
                                            unsigned int tray, bool top,
                                            int badness, bool allBad,
                                            const CalibData::StripCol& strips)
{
    *m_log << MSG::DEBUG;
    if((*m_log).isActive()) {
        *m_log << "BadVisitor::badPlane called with args" << endreq
        << "row = " << row << ", col = " << col << ", tray = "
        << tray << endreq    
        << "top = " << top << ", badness = " 
        << badness << " allBad = " << allBad << endreq
        << "Strip collection contains " << strips.size()
        << " strips. ";
    }
    *m_log << endreq;
    
    if (!allBad) { 
        unsigned int i;
        int tower = idents::TowerId(row, col).id();
        int layer = top ? tray : tray-1;
        int view  = layer%2 ? 1-top : top;
        int tag = 1;
        idents::GlastAxis::axis iview = view ? idents::GlastAxis::Y : idents::GlastAxis::X;

        TkrBadStripsSvc* pBad = dynamic_cast<TkrBadStripsSvc*>(m_pBadStrips);
        
        stripCol* v;
        v = const_cast<stripCol*> (pBad->getBadStrips(tower, layer, iview));
        if(v==0) {
            std::cout << " error in plane specification! Plane will be skipped!" << std::endl;
            return CalibData::CONT;
        }
        
        for (i=0;i<strips.size();i++) {
            int strip = strips[i];
            pBad->addStrip(v, TaggedStrip(strip, tag));
            m_nStrips++;
        }

        // removing duplicates
        std::sort(v->begin(), v->end());     
        stripCol::iterator p = std::unique(v->begin(), v->end());
        v->erase(p, v->end());
    }
        
    return CalibData::CONT;
}
