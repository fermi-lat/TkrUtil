/*
@file TkrSplitsSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.cxx,v 1.6 2004/08/24 23:45:45 lsrea Exp $

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
#include "xml/IFile.h"

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

    declareProperty("splitsFile", m_splitsFile="");
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

    // Open the message log
    MsgStream log( msgSvc(), name() );

    Service::initialize();
      
    m_pGeoSvc = 0;
    if( service( "TkrGeometrySvc", m_pGeoSvc, true).isFailure() ) {
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
    sc = doInit();
    
    return sc;
}
 
int TkrSplitsSvc::getSplitPoint(int tower, int layer, int view) const 
{
    return getLastC0Strip(tower, layer, view);
}

int TkrSplitsSvc::getEnd(int tower, int layer, 
                         int view, int strip) const 
{
    return ( (strip<=getLastC0Strip(tower, layer, view) ? 0 : 1) );
}
        
int TkrSplitsSvc::getLastC0Strip(int tower, int layer, int view) const
{
    // 
    MsgStream log( msgSvc(), name() );

    if (m_splitsFile!="") {
        return m_splits[tower][layer][view];
    } else if (!m_pSplits) {
        return defaultSplit;
    } else {
        int tray, botTop;
        m_pGeoSvc->layerToTray(layer, view, tray, botTop);
        bool isTop = (botTop==1);
        idents::TowerId twr(tower);
        int towerX = twr.ix();
        int towerY = twr.iy();
        idents::TkrId thisPlane(towerY, towerX, tray, isTop);
        CalibData::RangeBase* pPlane = m_pSplits->getChannel(thisPlane);
        CalibData::TkrSplit* pSplit = dynamic_cast<CalibData::TkrSplit*>(pPlane);
        int highChip = pSplit->getHigh();
        /*
        log << MSG::INFO 
            << "Tower " << tower << " Tray " << tray << " botTop " << botTop 
            << "High chip = " << highChip << endreq;
        */

        return pSplit->getHigh()*NSTRIPS - 1;
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

    int tower, layer, view;
    for(tower=0;tower<NTOWERS;++tower) {
        for (layer=0;layer<NLAYERS;++layer) {
            for (view=0;view<NVIEWS;++view) {
                m_splits[tower][layer][view] = defaultSplit;
            }
        }
    }
   
    if (m_splitsFile!="") {
        xml::IFile myFile(m_splitsFile.c_str());

        char buffer[8];
        for (tower=0;tower<NTOWERS;++tower) {
            sprintf(buffer,"Tower%d",tower);
            if(myFile.contains(buffer,"splits")) {
                std::vector<int> splits = myFile.getIntVector(buffer, "splits");
                int size = splits.size();
                if(size!=NLAYERS*NVIEWS) {
                    log << MSG::ERROR << buffer << ": splits vector size: " << splits.size()                      
                        << ", should be " << NLAYERS*NVIEWS << endreq;
                    return StatusCode::FAILURE;
                }
                int plane;
                for (plane=0; plane<NLAYERS*NVIEWS; ++plane) {
                    m_pGeoSvc->planeToLayer(plane, layer, view);

                    int split = splits[plane];
                    if(split<-1 || split>(NCHIPS-1)) {
                        log << MSG::ERROR << buffer << ": invalid split value " 
                            << split << endreq;
                        return StatusCode::FAILURE;
                    }
                    m_splits[tower][layer][view] = (split+1)*NSTRIPS -1;
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
