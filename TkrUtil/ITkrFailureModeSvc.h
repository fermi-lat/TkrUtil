#ifndef ITkrFailureModeSvc_H
#define ITkrFailureModeSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrFailureModeSvc("ITkrFailureModeSvc", 1 , 0);

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
    virtual int getFailureConditions()=0;

    /// look for object in list of failed objects
    virtual bool isFailed(int towerId, int layer = -1, int view = -1) = 0;
};

#endif // ITkrFailureModeSvc_H

