/** 
@file TkrBadStripsSvc.cxx

@brief Maintains lists of bad strips, and provides access methods 

First version 3-Jun-2001
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrBadStripsSvc.cxx,v 1.23 2006/03/21 01:15:48 usher Exp $
*/


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "src/TkrBadStripsSvc.h"
#include "idents/TowerId.h"
#include "TkrUtil/ITkrMakeClustersTool.h"

#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IToolSvc.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "facilities/Util.h"

static const SvcFactory<TkrBadStripsSvc> s_factory;
const ISvcFactory& TkrBadStripsSvcFactory = s_factory;

namespace {
    const std::string mType[2] = {"SIM", "REC"};
}

// Service parameters which can be set at run time must be declared.
// This should be done in the constructor.

TkrBadStripsSvc::TkrBadStripsSvc(const std::string& name, 
                                 ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the file to get data from   
    declareProperty("badStripsFile", m_commonBadStripsFile="");
    declareProperty("simBadStripsFile", m_badStripsFile[0]="");
    declareProperty("recBadStripsFile", m_badStripsFile[1]="");
    //declareProperty("killDigi", m_killDigi=false );
    m_visitor = 0;

    return;
}

StatusCode TkrBadStripsSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    Service::initialize();

    m_tkrGeom = 0;
    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get TkrGeometrySvc" 
            << endreq;
        return sc;
    }

    if(m_visitor==0) {
        m_visitor = new BadVisitor;
        m_visitor->setService(this);
        m_visitor->setService(m_tkrGeom);
    }

    m_commonBadStripsFile = "";
    m_badStripsFile[0] = "";
    m_badStripsFile[1] = "";

    setProperties();

    if(m_commonBadStripsFile!="") {
        m_badStripsFile[0] = m_commonBadStripsFile;
        m_badStripsFile[1] = m_commonBadStripsFile;
    }

    int nTowers = m_tkrGeom->numXTowers()*m_tkrGeom->numYTowers();
    int nLayers = m_tkrGeom->numLayers();
    int nViews  = 2;

    calibType type;
    int i;
    for (i=SIM; i<NCALIBTYPES; ++i) {
        type = (calibType) i;
        SetCalibType(type);
        m_empty[type] = true;
        bool valid = m_stripsCol[type].setDims(nTowers, nLayers, nViews);
        if (!valid) return StatusCode::FAILURE;
        stripCol null(0);
        m_stripsCol[type].setValue(null);

        // this method resolves environmental variables in the file name
        if (m_badStripsFile[m_calibType]!="") {
            int ret =  facilities::Util::expandEnvVar(&m_badStripsFile[m_calibType]);
            if (ret>=0) {
                log << MSG::INFO << "Input file for bad strips(type = " << mType[m_calibType] << ") " 
                    << m_badStripsFile[m_calibType] << endreq;
            } else {
                log << MSG::ERROR << "Input filename " << m_badStripsFile[m_calibType] 
                << " not resolved" << endreq;
                return StatusCode::FAILURE;
            }
        }

        m_pBadDigi[type] = 0;
        m_pBadClus[type] = 0;
        m_pBadMap[type]  = 0;

        sc = doInit();
        if (sc.isFailure()) return sc;
    }

    return sc;
}

StatusCode TkrBadStripsSvc::doInit()
{
    // Purpose: reads in bad strips and constructs in-memory bad strip vectors
    // Inputs:  None
    // Outputs: Status code (Success/Failure)

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    Event::TkrClusterCol* pClusters = new Event::TkrClusterCol(0);
    Event::TkrIdClusterMap* pMap = new Event::TkrIdClusterMap;
    setBadClusterCol(pClusters);
    setBadIdClusterMap(pMap);

    // If there is no bad strips file, service will do nothing at init time
    if (m_badStripsFile[m_calibType]=="") {        
        log << MSG::INFO << "No bad strips file was requested for type " 
            << mType[m_calibType] << "." << endreq;
        return sc;
    }

    // open bad strips file
    std::ifstream file;
    file.open( m_badStripsFile[m_calibType].c_str());

    if (!file) {
        log << MSG::ERROR << "  File not found: check jobOptions." << endreq;
        return StatusCode::FAILURE;
    }

    //int size = m_ntowers*m_nlayers*m_nviews;    
    //makeCol(size);

    readFromFile(&file);
    sc = generateBadClusters();

    file.close();

    return sc;
}

StatusCode TkrBadStripsSvc::finalize()
{
    delete m_visitor;
    return StatusCode::SUCCESS;
}

StatusCode TkrBadStripsSvc::update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot)
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    m_visitor->setLog(&log);
    log << MSG::INFO << "Updater called " << endreq;
    if (m_badStripsFile[m_calibType]=="") {
        if( pDead) pDead->traverse(m_visitor);
        if (pHot)  pHot->traverse(m_visitor);
        m_empty[m_calibType] = m_empty[m_calibType] && m_visitor->isEmpty();
    } else {
        log << MSG::INFO 
            << "No update done -- badStripsFile is being used instead" << endreq;
    }

    // generate bad clusters if there's been an update
    // this is probably the place to do it.
    sc = generateBadClusters();
    return sc;
}
StatusCode TkrBadStripsSvc::generateBadClusters()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    if(m_empty) { return sc; }

    IToolSvc* toolSvc = 0;
    if ((sc=service("ToolSvc",toolSvc, true)).isFailure() ){
        log << MSG::ERROR << "Couldn't fine ToolSvc" << endreq;
        return sc;
    }
    // not sure why this is needed yet...
    //sc = StatusCode::SUCCESS;

    Event::TkrDigiCol* pDigis = new Event::TkrDigiCol;
    makeBadDigiCol(pDigis); // deletes old one if it exists
    if (pDigis) {
        Event::TkrClusterCol* pClusters = new Event::TkrClusterCol(0);
        Event::TkrIdClusterMap* pMap = new Event::TkrIdClusterMap;

        ITkrMakeClustersTool* pMakeClusters;
        if ((sc=toolSvc->retrieveTool("TkrMakeClustersTool", pMakeClusters)).isFailure()) {
            log << MSG::ERROR << "Couldn't retrieve TkrMakeClusterTool" << endreq;
            return sc;
        }

        if((sc=pMakeClusters->makeClusters(pClusters, pMap, pDigis, BADCLUSTERS)).isFailure()) {
            log << MSG::ERROR << "makeClusters failed" << endreq;
            return sc;
        }

        setBadClusterCol(pClusters);
        //int size = pMap->size();
        setBadIdClusterMap(pMap);
    }

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
    // tower# plane#  [strip#] [strip#] ...  -1
    // all whitespace is ignored, so data for a layer may span lines
    //
    // 

    while(read && !file->eof()) {
        int tower;
        int plane;

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
        *file >> plane;

        // kludge until the geometry supplies this info
        // converts layer (0...35) to bilayer (digi format) and view

        //int layer = plane/2;
        //int element = (plane+3)%4;
        //int view = element/2;

        int layer, view;
        m_tkrGeom->planeToLayer(plane, layer, view);

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
        if (!v->empty()) {m_empty[m_calibType] = false;}       
    }
    return;
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

    //int index = getIndex(tower, layer, axis);

    return &m_stripsCol[m_calibType](tower, layer, axis);
}


bool TkrBadStripsSvc::isBadStrip(int tower, int layer, 
                                 idents::GlastAxis::axis axis, 
                                 int strip) const 
{
    // Purpose: determine if a given strip is bad
    // Inputs:  tower, bilayer, axis, strip#
    // Output:  true if strip is in the list ( that is, is bad)

    const stripCol* v = &m_stripsCol[m_calibType](tower, layer, axis);
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

StatusCode TkrBadStripsSvc::makeBadDigiCol(Event::TkrDigiCol* pDigis) 
{ 
    MsgStream log( msgSvc(), name() );

    // loop over all the bad strips and make an empty digi for each plane
    if (m_pBadDigi[m_calibType]) { delete m_pBadDigi[m_calibType]; }

    int nTowers = m_stripsCol[0].getDim(0);
    int nLayers = m_stripsCol[0].getDim(1);
    int nViews  = m_stripsCol[0].getDim(2);

    m_pBadDigi[m_calibType] = pDigis;
    int tower, layer, view;
    for (tower=0; tower<nTowers; ++tower) {
        for (layer=0; layer<nLayers; ++layer) {
            for (view=0; view<nViews; ++view) {
                idents::GlastAxis::axis axis = ( view ? idents::GlastAxis::X : idents::GlastAxis::Y);
                const stripCol* strips = getBadStrips(tower, layer, axis);
                if (strips==0 || strips->size()==0) continue;
                int tot[2] = { 0, 0};
                Event::TkrDigi* digi = new Event::TkrDigi(layer, axis, idents::TowerId(tower), tot);
                m_pBadDigi[m_calibType]->push_back(digi);
            }
        }
    }
    log << MSG::INFO << m_pBadDigi[m_calibType]->size() <<" digis created for bad cluster processing" << endreq;

    return StatusCode::SUCCESS;
}

// queryInterface

StatusCode  TkrBadStripsSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
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

const InterfaceID&  TkrBadStripsSvc::type () const {
    return IID_ITkrBadStripsSvc;
}

CalibData::eVisitorRet BadVisitor::badTower
(unsigned int /* row */, 
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
         int tower = idents::TowerId(col, row).id();
         //int layer = top ? tray : tray-1;
         //int view  = layer%2 ? 1-top : top;
         int layer, view;
         m_tkrGeom->trayToLayer(tray, top, layer, view);
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
