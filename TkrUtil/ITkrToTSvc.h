/** @file ITkrToTSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /home/cvs/SLAC/TkrUtil/TkrUtil/ITkrToTSvc.h,v 1.3.2.1 2004/12/14 02:57:14 lsrea Exp $
*/

#ifndef ITkrToTSvc_H
#define ITkrToTSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"
#include "idents/VolumeIdentifier.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrToTSvc("ITkrToTSvc", 3 , 0);

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
    virtual double getGain2     (const int tower, const int layer, const int view,
        const int strip) const = 0;
    virtual double getQuality   (const int tower, const int layer, const int view,
        const int strip) const = 0;
    virtual double getCountsPerMicrosecond () const = 0;
    virtual double getMevPerMip() const = 0;
    virtual double getFCPerMip() const = 0;
    virtual int    getMaxToT() const = 0;
    
    virtual double getCharge(double ToT, int tower, int layer, int view, int strip) const = 0;
    virtual double getMipsFromToT(double ToT, int tower, int layer, int view, int strip) const = 0;
    virtual double getMipsFromCharge(double charge, int tower, int layer, int view, int strip) const = 0;
};

#endif // ITkrToTSvc_H
