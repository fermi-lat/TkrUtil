/*
@file TkrSplitsSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.h,v 1.3 2004/03/13 19:40:37 lsrea Exp $

*/
#ifndef TkrSplitsSvc_H
#define TkrSplitsSvc_H 1

// Include files
#include "TkrUtil/ITkrSplitsSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"

/** @class TkrSplitsSvc
* @brief Service to retrieve the splits constants
* the TKR.
*
* Author:  L. Rochester
*
*/

class TkrSplitsSvc : public Service, virtual public ITkrSplitsSvc  {

public:

    TkrSplitsSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    //StatusCode execute();
    StatusCode finalize();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrSplitsSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

    /// get the last C0 strip for this layer
    int getSplitPoint(int tower, int layer, int view) const;

    /// tell which end this strip belongs to
    int getEnd(int tower, int layer, int view, int strip) const;

    /// update the pointer
    void update(CalibData::TkrSplitsCalib* pSplits);

private:
    /// internal init method
    StatusCode doInit();
    /// get the constant
    int getLastC0Strip(int tower, int layer, int view) const;
   /// pointer to data provider svc
    IDataProviderSvc* m_pCalibDataSvc;
    /// pointer to the geometry
    ITkrGeometrySvc* m_geoSvc;
    /// pointer to the calibration data
    CalibData::TkrSplitsCalib* m_pSplits;
};

#endif // TkrSplitsSvc_H

