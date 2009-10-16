/**
 * @class ITkrDiagnosticTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#ifndef ITkrDiagnosticTool_h
#define ITkrDiagnosticTool_h

namespace {
    // assign the ends to the 8 cables
    const int endArray[8] = {0, 1, 1, 0, 1, 0, 0, 1};

    // set the fields in the indices
    const int planeMult = 2;
    const int viewMult  = 1;
    const int layerMult = 2;
    const int geoMult   = 2;
    const int towerMult = 1000;
    const int gtccMult  = 10;

    const int nRc       = 9;
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
