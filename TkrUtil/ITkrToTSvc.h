/** @file ITkrToTSvc.h
@brief Abstract interface to TkrSplitsSvc (q.v.)
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrToTSvc.h,v 1.9 2005/12/20 02:35:57 lsrea Exp $
*/

#ifndef ITkrToTSvc_H
#define ITkrToTSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"
#include "idents/VolumeIdentifier.h"
#include "idents/TkrId.h"

#include "CalibData/Tkr/TkrTot.h"
#include "CalibData/Tkr/TkrScale.h"


// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrToTSvc("ITkrToTSvc", 4 , 0);

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
    virtual double getQuad       (int tower, int layer, int view, int strip) const = 0;
    virtual double getQuality    (int tower, int layer, int view, int strip) const = 0;
    virtual double getMuonScale  (int tower, int layer, int view, int strip) const = 0;
    virtual double getCountsPerMicrosecond () const = 0;
    virtual double getFCPerMip()  const = 0;
    virtual double getMevPerMip() const = 0;
    virtual int    getMaxToT() const = 0;
    
    virtual double getCharge(double rawToT, int tower, 
        int layer, int view, int strip) const = 0;
    virtual int    getRawToT(double eDep, int tower, 
        int layer, int view, int strip) const = 0;
    virtual double getMipsFromToT(double rawToT, int tower, 
        int layer, int view, int strip) const = 0;
    virtual double getMipsFromCharge(double charge) const = 0;

        /// update to latest pointer when calibration changes
    virtual void update(CalibData::TkrTotCol* pToT) = 0;
    virtual void update(CalibData::TkrScaleCol* pScale) = 0;

};

#endif // ITkrToTSvc_H
