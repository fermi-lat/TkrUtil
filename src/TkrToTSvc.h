/*
@file TkrToTSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrToTSvc.h,v 1.2 2004/03/12 05:49:22 lsrea Exp $

*/
#ifndef TkrToTSvc_H
#define TkrToTSvc_H 1

// Include files
#include "TkrUtil/ITkrToTSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"

/** @class TkrToTSvc
* @brief Service to store and compare to a list of desired failure modes in
* the TKR.
*
* Author:  L. Rochester (after R.Dubois)
*
*/

class TkrToTSvc : public Service, virtual public ITkrToTSvc  {

public:

    enum {NTOWERS=16, NLAYERS=18, NVIEWS=2, NSTRIPS=1536};

    TkrToTSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode finalize();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrToTSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

    double getGain(const int tower, const int layer, const int view, 
        const int strip) const 
    {
        if (tower>-1 && tower <NTOWERS && layer>-1 && layer<NLAYERS
            && view>-1 && view<NVIEWS && strip>-1 && strip<NSTRIPS) 
        {
            return (double) m_ToTGain[tower][layer][view][strip];
        }else {
            return -1.;
        }
    }
    double getThreshold(const int tower, const int layer, const int view, 
        const int strip) const
    {
        if (tower>-1 && tower <NTOWERS && layer>-1 && layer<NLAYERS
            && view>-1 && view<NVIEWS && strip>-1 && strip<NSTRIPS) 
        {
            return (double) m_ToTThreshold[tower][layer][view][strip];
        }else {
            return -1.;
        }
    }

private:
    /// internal init method
    StatusCode doInit();
    /// name of file containing splits
    std::string m_ToTFile;
    /// default Gain
    double m_defaultGain;
    /// default Threshold
    double m_defaultThreshold;
    /// array of gains, in microseconds/fC
    float m_ToTGain      [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Thresholds, in microseconds = extrapolation to zero charge
    float m_ToTThreshold [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// pointer to geometry service
    ITkrGeometrySvc* m_geoSvc;
};


#endif // TkrToTSvc_H

