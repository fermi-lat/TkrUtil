/*
@file TkrSplitsSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.h,v 1.7 2005/04/11 22:52:02 lsrea Exp $

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

    enum {NTOWERS=16, NTRAYS =19, NLAYERS=18, NVIEWS=2, NFACES=2, NPLANES=36};
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

    /// get the last C0 strip for this plane
    int getSplitPoint(int tower, int layer, int view) const;

    /// tell which end this strip belongs to
    int getEnd(int tower, int tray, int layer, int strip) const;

    /// update the pointer
    void update(CalibData::TkrSplitsCalib* pSplits);

    /// get max hits
    int getMaxStrips(int tower, int layer, int view, int end) const;

    /// get the cable buffer size
    int getCableBufferSize() const { return m_cableBuffer; }

    int getCableIndex(int layer, int view, int end) const {
        return cableIndex[layer%2][view][end];
    }

private:
    /// internal init method
    StatusCode doInit();
    /// get the constant
    int getLastC0Strip(int tower, int layer, int view) const;
    /// pointer to data provider svc
    IDataProviderSvc* m_pCalibDataSvc;
    /// pointer to the geometry
    ITkrGeometrySvc* m_tkrGeom;
    /// pointer to the calibration data
    CalibData::TkrSplitsCalib* m_pSplits;
    /// name of the input file, if present
    std::string m_splitsFile;
    /// array containing splits, for use as a quick test
    int m_splits[NTOWERS][NTRAYS][NFACES];
    /// default maxStrips
    int m_defaultMaxStrips;
    /// File containing the max strips information
    std::string m_maxStripsFile;
    /// full maxStrips
    int m_maxStrips[NTOWERS][NTRAYS][NFACES][2];
    /// cable buffer size
    int m_cableBuffer;
};

#endif // TkrSplitsSvc_H

