/** @file ITkrToTSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrToTSvc.h,v 1.1 2004/03/13 19:40:37 lsrea Exp $
*/

#ifndef ITkrToTSvc_H
#define ITkrToTSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"
#include "idents/VolumeIdentifier.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrToTSvc("ITkrToTSvc", 1 , 0);

/** @class ITkrToTSvc
* @brief Interface class for TkrSplitsSvc
*
* Author:  L. Rochester
*
*/

class ITkrToTSvc : virtual public IInterface {

public:
    static const InterfaceID& interfaceID() { return IID_ITkrToTSvc; }

    /// interface methods here
    virtual double getGain      (const int tower, const int layer, const int view,
        const int strip) const = 0;
    virtual double getThreshold (const int tower, const int layer, const int view,
        const int strip) const = 0;
};

#endif // ITkrToTSvc_H

