/** @file ITkrToTSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrToTSvc.h,v 1.4 2004/12/16 23:28:30 usher Exp $
*/

#ifndef ITkrToTSvc_H
#define ITkrToTSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"
#include "idents/VolumeIdentifier.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrToTSvc("ITkrToTSvc", 3 , 0);

/** @class ITkrToTSvc
* @brief Interface class for TkrToTSvc
*
* Author:  L. Rochester
*
*/

class ITkrToTSvc : virtual public IInterface {

public:
    static const InterfaceID& interfaceID() { return IID_ITkrToTSvc; }

    /// interface methods here
    virtual double getGain       (int tower, int layer, int view, int strip) const = 0;
    virtual double getThreshold  (int tower, int layer, int view, int strip) const = 0;
    virtual double getGain2      (int tower, int layer, int view, int strip) const = 0;
    virtual double getQuality    (int tower, int layer, int view, int strip) const = 0;
    virtual double getMuonFactor (int tower, int layer, int view, int strip) const = 0;
    virtual double getCountsPerMicrosecond () const = 0;
    virtual double getFCPerMip()  const = 0;
    virtual double getMevPerMip() const = 0;
    virtual int    getMaxToT() const = 0;
    
    virtual double getCharge(double ToT, int tower, int layer, int view, int strip) const = 0;
    virtual double getMipsFromToT(double ToT, int tower, int layer, int view, int strip) const = 0;
    virtual double getMipsFromCharge(double charge, int tower, int layer, int view, int strip) const = 0;
};

#endif // ITkrToTSvc_H
