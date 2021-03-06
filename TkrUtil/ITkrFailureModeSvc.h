/** @file ITkrFailureModeSvc.h
@brief Abstract interface to TkrFailureModeSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrFailureModeSvc.h,v 1.5 2004/08/24 23:45:44 lsrea Exp $
*/



#ifndef ITkrFailureModeSvc_H
#define ITkrFailureModeSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"

#include "idents/TkrId.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrFailureModeSvc("ITkrFailureModeSvc", 3, 0);

/** @class ITkrFailureModeSvc
* @brief Interface class for TkrFailureModeSvc
*
* Author:  L. Rochester (after R. Dubois)
*
*/

class ITkrFailureModeSvc : virtual public IInterface {

public:
    static const InterfaceID& interfaceID() { return IID_ITkrFailureModeSvc; }

    /// get the list of enabled failure mode conditions
    virtual int getFailureConditions() const =0;

    /// look for object in list of failed objects
    virtual bool isFailed(int towerId, int layer = -1, int view = -1) const = 0;
    virtual bool isFailed(const idents::TkrId& tkrId ) const = 0;
    virtual bool empty() const = 0;
};

#endif // ITkrFailureModeSvc_H

