/** @file ITkrFailureModeSvcCalib.h
@brief Abstract interface to TkrFailureModeSvc, used by TkrCalibAlg
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrAlignmentSvc.h,v 1.6 2003/04/08 22:20:05 lsrea Exp $
*/



#ifndef ITkrFailureModeSvcCalib_H
#define ITkrFailureModeSvcCalib_H 1

// Include files
#include "GaudiKernel/IInterface.h"

#include "TkrUtil/ITkrFailureModeSvc.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/BadStrips.h"

#include <map>


// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_ITkrFailureModeSvcCalib("ITkrFailureModeSvcCalib", 2 , 0);

typedef std::map <int, std::vector<int> > LayerMap;


/** @class ITkrFailureModeSvcCalib
* @brief An interface class for TkrFailureModeSvc
*
* Author:  L. Rochester (after R. Dubois)
*
*/

class ITkrFailureModeSvcCalib : virtual public IInterface, virtual public ITkrFailureModeSvc 
{

public:
    static const InterfaceID& interfaceID() { return IID_ITkrFailureModeSvcCalib; }

    virtual StatusCode update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot) = 0;

    virtual std::vector<int>& getLayers(int tower) = 0;
    virtual std::vector<int>& getTowers() = 0;
};

#endif // ITkrFailureModeSvcCalib_H

