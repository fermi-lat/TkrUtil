/** @file ITkrSplitsSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrSplitsSvc.h,v 1.4 2003/04/11 23:27:15 lsrea Exp $
*/

#ifndef ITkrSplitsSvc_H
#define ITkrSplitsSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"

#include "TkrUtil/ITkrGeometrySvc.h"

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

    /// get the list of enabled failure mode conditions
    virtual int  getEnd(int tower, int layer, int view, int strip) const =0;
    virtual int  getSplitPoint(int tower, int layer, int view) const = 0;
};

#endif // ITkrSplitsSvc_H

