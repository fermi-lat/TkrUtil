/** @file ITkrSplitsSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrSplitsSvc.h,v 1.2 2004/03/13 19:40:37 lsrea Exp $
*/

#ifndef ITkrSplitsSvc_H
#define ITkrSplitsSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/TkrSplitsCalib.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrSplitsSvc("ITkrSplitsSvc", 1 , 0);

/** @class ITkrSplitsSvc
* @brief Interface class for TkrSplitsSvc
*
* Author:  L. Rochester
*
*/

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

};

#endif // ITkrSplitsSvc_H

