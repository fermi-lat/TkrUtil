/** @file ITkrSplitsSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrSplitsSvc.h,v 1.4 2005/04/11 22:52:01 lsrea Exp $
*/

#ifndef ITkrSplitsSvc_H
#define ITkrSplitsSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/TkrSplitsCalib.h"
#include "Event/Digi/TkrDigi.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrSplitsSvc("ITkrSplitsSvc", 1 , 1);

/** @class ITkrSplitsSvc
* @brief Interface class for TkrSplitsSvc
*
* Author:  L. Rochester
*
*/

namespace {
    // this codes the correspondence between physical space and cable space
    enum {NRANGE=2, NVIEW=2, NEND=2};
    const int cableIndex[NRANGE][NVIEW][NEND] = {3,2,0,1,6,7,5,4};
    }


class ITkrSplitsSvc : virtual public IInterface {

public:
    static const InterfaceID& interfaceID() { return IID_ITkrSplitsSvc; }

    /// get the controller for this strip
    virtual int  getEnd(int tower, int layer, int view, int strip) 
        const = 0;
    /// get the split point for this plane (lastC0Strip)
    virtual int  getSplitPoint(int tower, int layer, int view) 
        const = 0;
    /// update to latest pointer when calibration changes
    virtual void update(CalibData::TkrSplitsCalib* pSplits) = 0;
    /// max number of hits
    virtual int  getMaxStrips(int tower, int layer, int view, int end=0)
        const = 0;
    virtual int  getCableBufferSize() const = 0;
    /// cable number for a given plane and end
    virtual int getCableIndex(int bilayer, int view, int end) const = 0;

};

#endif // ITkrSplitsSvc_H

