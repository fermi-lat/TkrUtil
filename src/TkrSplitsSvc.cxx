/*
@file TkrSplitsSvc.cxx

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.cxx,v 1.2 2004/03/12 05:49:22 lsrea Exp $

*/

#include "src/TkrSplitsSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "idents/TowerId.h"

#include <algorithm>
#include <fstream>
#include "facilities/Util.h"
#include "xml/IFile.h"

// declare the service factories for the TkrSplitsSvc
static SvcFactory<TkrSplitsSvc> a_factory;
const ISvcFactory& TkrSplitsSvcFactory = a_factory; 

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

    declareProperty("splitsFile", m_splitsFile= "");
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
    
    StatusCode  status = StatusCode::SUCCESS;

    // Open the message log
    MsgStream log( msgSvc(), name() );
      
    // Call super-class
    Service::initialize ();

    m_geoSvc = 0;
    if( service( "TkrGeometrySvc", m_geoSvc, true).isFailure() ) {
        log << MSG::ERROR << "Couldn't retrieve TkrGeometrySvc" << endreq;
        return StatusCode::FAILURE;
    }
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
        

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

    status = doInit();
    
    return status;
}

StatusCode TkrSplitsSvc::doInit()
{
    // Open the message log
    MsgStream log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;

    // test of getting TkrGeometrySvc from inside TkrSplitsSvc... It works!
    //int stripsPerLadder  = m_geoSvc->ladderNStrips();

    // can be removed when geometry is iterfaced here
    const int NSTRIPS = 64;
    const int NCHIPS  = 24;
    const int defaultSplit = (NCHIPS/2)*NSTRIPS - 1;

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
                    // fix this kludge by adding a method to TkrGeometry,
                    //    or maybe TkrTowerID??
                    layer = plane/2;
                    int element = (plane+3)%4;
                    view = element/2;

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
