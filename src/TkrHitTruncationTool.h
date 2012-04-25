/** @file TkrHitTruncationTool.h
*/

/**
* @class TkrHitTruncationTool
*
* @brief This tool analyzes the digis to infer truncation
*        
* File and Version Information:
*      $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrHitTruncationTool.h,v 1.3 2012/01/20 19:22:56 lsrea Exp $
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
#include "TkrUtil/ITkrDiagnosticTool.h"
#include "TkrUtil/ITkrMapTool.h"
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
    void addEmptyDigis();
    void removeEmptyDigis();
    StatusCode trimDigis();
    void removeEmptyTruncs();
    void setTrimCount(int trimCount) { m_trimCount = trimCount; }

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeom;
    /// splits service
    ITkrSplitsSvc*      m_splitsSvc;
     /// Pointer to the Gaudi data provider service
    IDataProviderSvc*   m_dataSvc;
    ///
    IGlastDetSvc*       m_detSvc;
    ITkrDiagnosticTool* m_diagTool;
    ITkrMapTool*        m_mapTool;

    bool m_newEvent;

    bool m_trimDigis;
    int  m_trimCount;

    Event::TkrTruncationInfo::TkrTruncationMap* m_truncMap;
    Event::TkrDigiCol*                          m_digiCol;
    Event::TkrTruncationInfo*                   m_truncationInfo;

    /// this is called by the incident service at the beginning of an event
    void handle(const Incident& inc);
    bool m_trimDigis;
    int  m_trimCount;

    void doRCLoop();
    StatusCode setPointers();
    StatusCode newEvent();

};

//static ToolFactory<TkrHitTruncationTool> s_factory;
//const IToolFactory& TkrHitTruncationToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrHitTruncationTool);

#endif
