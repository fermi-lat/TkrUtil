/** @file ITkrBadStripsSvcCalib.h
@brief Abstract Interface to TkrBadStripsSvc, used by TkrCalibAlg
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrBadStripsSvcCalib.h,v 1.2 2003/04/11 23:27:15 lsrea Exp $
*/



#ifndef __ITKRBADSTRIPSSVCCALIB_H
#define __ITKRBADSTRIPSSVCCALIB_H 1

#include "GaudiKernel/IInterface.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/BadStrips.h"


//----------------------------------------------
//
//   ITkrBadStripsSvcCalib
//
//   Tracker BadStrips Service. This is the interface
//   for the calibration database.
//----------------------------------------------
//             Leon Rochester, 3-June-2001
//----------------------------------------------

static const InterfaceID IID_ITkrBadStripsSvcCalib("ITkrBadStripsSvcCalib", 2 , 0); 



class ITkrBadStripsSvcCalib //: public virtual ITkrBadStripsSvc
{
public:

    //! Constructor of this form must be provided

	static const InterfaceID& interfaceID() { return IID_ITkrBadStripsSvcCalib; }

	virtual StatusCode update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot) = 0;
 };

#endif
