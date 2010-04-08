/**
 * @class ITkrDiagnosticTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrDiagnosticTool.h,v 1.1 2009/10/16 18:47:30 lsrea Exp $
 */
#ifndef ITkrDiagnosticTool_h
#define ITkrDiagnosticTool_h

namespace {
 
    // set the fields in the indices
    const int viewMult  = 1;
    const int layerMult = 2;
    const int towerMult = 1000;
}

#include "GaudiKernel/IAlgTool.h"

#include "TkrUtil/TkrTowerBits.h"

static const InterfaceID IID_ITkrDiagnosticTool("ITkrDiagnosticTool", 1 , 0);

class ITkrDiagnosticTool : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrDiagnosticTool; }

    virtual StatusCode getTkrVector(unsigned short& tkrVector) = 0;
    virtual StatusCode getTkrDiagnosticData() = 0;
    virtual StatusCode calculateTkrVector(unsigned short& towerBits) = 0;
        
    virtual void setTriggerInfo(int tower, int gtcc, int gtrc) = 0;
    virtual bool isSetTrigger(int tower, int plane, int end) = 0;
    virtual void clearTriggerInfo() = 0;
};

#endif
