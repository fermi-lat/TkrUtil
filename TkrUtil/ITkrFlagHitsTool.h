/** @file ITkrFlagHitsTool.h
@brief Abstract interface for methods to classify track hits
@author Leon Rochester
$Header$
*/

#ifndef _H_ITkrFlagHitsTool
#define _H_ITkrFlagHitsTool

#include "GaudiKernel/IAlgTool.h"

#include "geometry/Point.h"
#include "idents/TkrId.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrFlagHitsTool("ITkrFlagHitsTool", 0 , 0); 

/** @class ITkrFlagHitsTool
* @brief Abstract interface for methods to query TkrClusters
*/

class   ITkrFlagHitsTool : virtual public IAlgTool {
public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrFlagHitsTool; }

    //virtual int flagHits(Point pos, int tower, int layer, int view, 
    //    double xTower, double yTower,
    //    double xError, double yError, double hitWidth, double sigmaCut,
    //    double *paramList, unsigned int& status_bits) const = 0;

    virtual int flagHits(idents::TkrId tkrId, 
        Event::TkrTrackParams inParams, double zIn, 
        double minError, double maxError, double nSigma, 
        Event::TkrTrackParams& outParams,  
        unsigned int& status_bits) const = 0;
};

#endif  // _H_ITkrFlagHitsTool
