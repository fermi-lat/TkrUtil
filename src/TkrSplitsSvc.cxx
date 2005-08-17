/*
@file TkrSplitsSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.cxx,v 1.14 2005/08/17 18:38:49 lsrea Exp $

*/

#include "src/TkrSplitsSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "GaudiKernel/SmartDataPtr.h"
//#include "CalibData/CalibModel.h"
//#include "CalibData/Tkr/TkrSplitsCalib.h"

#include "idents/TowerId.h"
#include "idents/TkrId.h"

#include <algorithm>
#include <fstream>
#include "facilities/Util.h"
#include "xmlBase/IFile.h"

// declare the service factories for the TkrSplitsSvc
static SvcFactory<TkrSplitsSvc> a_factory;
const ISvcFactory& TkrSplitsSvcFactory = a_factory; 

namespace {
    // can be removed when geometry is iterfaced here
    const int NSTRIPS = 64;
    const int NCHIPS  = 24;
    const int defaultSplit = (NCHIPS/2)*NSTRIPS - 1;
}

TkrSplitsSvc::TkrSplitsSvc(const std::string& name,ISvcLocator* svc) 
: Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // declare the properties

    declareProperty("splitsFile",       m_splitsFile="");
    declareProperty("defaultMaxStrips", m_defaultMaxStrips=64);
    declareProperty("maxStripsFile",    m_maxStripsFile="");
    declareProperty("cableBufferSize",  m_cableBuffer=128);
}

StatusCode  TkrSplitsSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrSplitsSvc == riid) {
        *ppvIF = dynamic_cast<ITkrSplitsSvc*> (this);
        return StatusCode::SUCCESS;
    } else {
        return Service::queryInterface (riid, ppvIF);
    }
}

const IID&  TkrSplitsSvc::type () const {
    return IID_ITkrSplitsSvc;
}

StatusCode TkrSplitsSvc::initialize () 
{
    // Purpose and Method: Initialize the lists of dead units
    //                     
    // Inputs: None        
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode  sc = StatusCode::SUCCESS;

    Service::initialize();

    // Open the message log
    MsgStream log( msgSvc(), name() );

    Service::initialize();

    m_tkrGeom = 0;
    if( service( "TkrGeometrySvc", m_tkrGeom, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't retrieve TkrGeometrySvc" << endreq;
        return StatusCode::FAILURE;
    }

    sc = service("CalibDataSvc", m_pCalibDataSvc, true);

    if ( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Could not get IDataProviderSvc interface of CalibDataSvc" 
            << endreq;
        return sc;
    }

    // Bind all of the properties for this service
    if ( (sc = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }

    m_pSplits = 0;

    // back door for quick testing... this overrides the database access.
    if (m_splitsFile!="") {
        int ret =  facilities::Util::expandEnvVar(&m_splitsFile);
        if (ret>=0) {
            log << MSG::INFO << "Input file for splits: " 
                << m_splitsFile << endreq;
        } else {
            log << MSG::ERROR << "Input filename " << m_splitsFile << " not resolved" << endreq;
            return StatusCode::FAILURE;
        }
    }

    if (m_maxStripsFile!="") {
        int ret =  facilities::Util::expandEnvVar(&m_maxStripsFile);
        if (ret>=0) {
            log << MSG::INFO << "Input file for read controller buffer sizes: " 
                << m_maxStripsFile << endreq;
        } else {
            log << MSG::ERROR << "Input filename " << m_maxStripsFile << " not resolved" << endreq;
            return StatusCode::FAILURE;
        }
    }

    sc = doInit();

    return sc;
}

int TkrSplitsSvc::getSplitPoint(int tower, int layer, int view) const 
{
    return getLastC0Strip(tower, layer, view);
}

int TkrSplitsSvc::getEnd(int tower, int layer, int view, int strip) const 
{
    return ( (strip<=getLastC0Strip(tower, layer, view) ? 0 : 1) );
}

int TkrSplitsSvc::getLastC0Strip(int tower, int layer, int view) const
{
    // 
    MsgStream log( msgSvc(), name() );

    int tray, face;
    m_tkrGeom->layerToTray(layer, view, tray, face);
    if (m_splitsFile!="") {
        return m_splits[tower][tray][face];
    } else if (!m_pSplits) {
        return defaultSplit;
    } else {
        bool isTop = (face==1);
        idents::TowerId twr(tower);
        int towerX = twr.ix();
        int towerY = twr.iy();
        idents::TkrId thisPlane(towerY, towerX, tray, isTop);
        CalibData::RangeBase* pPlane = m_pSplits->getChannel(thisPlane);
        CalibData::TkrSplit* pSplit = dynamic_cast<CalibData::TkrSplit*>(pPlane);
        //int highChip = pSplit->getHigh();
        /*
        log << MSG::INFO 
        << "Tower " << tower << " Tray " << tray << " botTop " << botTop 
        << "High chip = " << highChip << endreq;
        */
        return pSplit->getHigh()*NSTRIPS - 1;
    }
}

int TkrSplitsSvc::getMaxStrips(int tower, int layer, int view, int end) const
{
    // for now, just return m_defaultStrips... we need something to handle contouring
    if(m_maxStripsFile=="") {
        return m_defaultMaxStrips;
    } else {
        int tray, face;
        m_tkrGeom->layerToTray(layer, view, tray, face);
        return m_maxStrips[tower][tray][face][end];
    }
}

void TkrSplitsSvc::update(CalibData::TkrSplitsCalib* pSplits)
{
    MsgStream log( msgSvc(), name() );

    m_pSplits = pSplits;

    log << MSG::INFO << "Splits pointer updated" << endreq;
}

StatusCode TkrSplitsSvc::doInit()
{
    // Open the message log
    MsgStream log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;

    // Keep this for a while

    int tower, tray, face;
    for(tower=0;tower<NTOWERS;++tower) {
        for (tray=0;tray<NTRAYS;++tray) {
            for (face=0;face<NFACES;++face) {
                m_splits[tower][face][face] = defaultSplit;
            }
        }
    }

    if (m_splitsFile!="") {
        xmlBase::IFile myFile(m_splitsFile.c_str());

        char buffer[8];
        for (tower=0;tower<NTOWERS;++tower) {
            sprintf(buffer,"Tower%d",tower);
            if(myFile.contains(buffer,"splits")) {
                std::vector<int> splits = myFile.getIntVector(buffer, "splits");
                int size = splits.size();
                if(size!=NTRAYS*NFACES) {
                    log << MSG::ERROR << buffer << ": splits vector size: " << splits.size()                      
                        << ", should be " << NLAYERS*NVIEWS << endreq;
                    return StatusCode::FAILURE;
                }
                int index = 0;
                for (tray=0; tray<NTRAYS; ++tray) {
                    for (face=0; face<NFACES; ++face,++index) {
                    int split = splits[index];
                    if(split<-1 || split>(NCHIPS-1)) {
                        log << MSG::ERROR << buffer << ": invalid split value " 
                            << split << endreq;
                        return StatusCode::FAILURE;
                    }
                    m_splits[tower][tray][face] = (split+1)*NSTRIPS -1;
                }
                }
            }
        }
    }

    int end;
    for(tower=0;tower<NTOWERS;++tower) {
        for (tray=0;tray<NTRAYS;++tray) {
            for (face=0;face<NFACES;++face) {
                for (end=0; end<2; ++end) {
                    m_maxStrips[tower][face][face][end] = m_defaultMaxStrips;
                }
            }
        }
    }

    if (m_maxStripsFile!="") {
        xmlBase::IFile myFile(m_maxStripsFile.c_str());

        char buffer[8];
        for (tower=0;tower<NTOWERS;++tower) {
            sprintf(buffer,"Tower%d",tower);
            if(myFile.contains(buffer,"maxStrips")) {
                std::vector<int> maxStrips = myFile.getIntVector(buffer, "maxStrips");
                int size = maxStrips.size();
                if(size!=NPLANES*2) {
                    log << MSG::ERROR << buffer << ": maxStrips vector size: " << maxStrips.size()                      
                        << ", should be " << NLAYERS*NVIEWS << endreq;
                    return StatusCode::FAILURE;
                }
                int plane, end;
                for (plane=0; plane<NPLANES; ++plane) {
                    for (end=0; end<2; ++end) {
                        tray = m_tkrGeom->planeToTray(plane);
                        face = m_tkrGeom->planeToBotTop(plane);
                        int maxStr = maxStrips[2*plane+end];
                        if(maxStr<0 || maxStr>64) {
                            log << MSG::ERROR << buffer << ": invalid maxStrips value " 
                                << maxStr << endreq;
                            return StatusCode::FAILURE;
                        }
                        m_maxStrips[tower][tray][face][end] = maxStr;
                    }
                }
            }
        }
    }

    return sc;
}

StatusCode TkrSplitsSvc::finalize( ) {
    MsgStream log(msgSvc(), name());  
    return StatusCode::SUCCESS;
}
