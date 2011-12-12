/** @file TkrHitTruncationTool.h
*/

/**
* @class TkrHitTruncationTool
*
* @brief This tool analyzes the digis to infer truncation
*        
* File and Version Information:
*      $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrUtil/src/TkrHitTruncationTool.h,v 1.1 2011/03/26 22:32:12 lsrea Exp $
*/


#ifndef TkrHitTruncationTOOL_H
#define TkrHitTruncationTOOL_H

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"

#include "TkrUtil/ITkrHitTruncationTool.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

namespace {
}

class TkrHitTruncationTool : public AlgTool, virtual public ITkrHitTruncationTool,
                             virtual public IIncidentListener
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrHitTruncationTool(const std::string& type, const std::string& name, 
        const IInterface* parent);
    ~TkrHitTruncationTool() {}

    StatusCode initialize();
    StatusCode analyzeDigis();
    StatusCode finalize();

    double getDistanceToTruncation(int tower, int tray, int face, int view, 
        double localX);
    double getDistanceToTruncation(idents::TkrId id, Vector towerPos);
    double getDistanceToTruncation(int tower, int plane, Vector towerPos);

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeom;
    /// splits service
    ITkrSplitsSvc*      m_splitsSvc;
     /// Pointer to the Gaudi data provider service
    IDataProviderSvc*   m_dataSvc;
    ///
    IGlastDetSvc*       m_detSvc;

    bool m_newEvent;

    Event::TkrTruncationInfo::TkrTruncationMap* m_truncMap;

    /// this is called by the incident service at the beginning of an event
    void handle(const Incident& inc);

};

//static ToolFactory<TkrHitTruncationTool> s_factory;
//const IToolFactory& TkrHitTruncationToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrHitTruncationTool);

#endif
