/*
@file TkrToTSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrToTSvc.h,v 1.3 2004/04/10 05:57:02 lsrea Exp $

*/
#ifndef TkrToTSvc_H
#define TkrToTSvc_H 1

// Include files
#include "TkrUtil/ITkrToTSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"
#include "idents/TowerId.h"

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
    double getGain2(const int tower, const int layer, const int view, 
        const int strip) const 
    {
        if (tower>-1 && tower <NTOWERS && layer>-1 && layer<NLAYERS
            && view>-1 && view<NVIEWS && strip>-1 && strip<NSTRIPS) 
        {
            return (double) m_ToTGain2[tower][layer][view][strip];
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
    double getQuality(const int tower, const int layer, const int view, 
        const int strip) const
    {
        if (tower>-1 && tower <NTOWERS && layer>-1 && layer<NLAYERS
            && view>-1 && view<NVIEWS && strip>-1 && strip<NSTRIPS) 
        {
            return (double) m_ToTQuality[tower][layer][view][strip];
        }else {
            return -1.;
        }
    }
    double getCountsPerMicrosecond() const { return m_countsPerMicrosecond;}

private:
    /// internal init method
    StatusCode doInit();
    /// mode: currently "default" or "EM"
    std::string m_mode;
    /// name of file containing splits
    std::string m_ToTFile;
    /// default Gain
    double m_defaultGain;
    /// default quadratic term;
    double m_defaultGain2;
    /// default Threshold
    double m_defaultThreshold;
    /// default quality factor
    double m_defaultQuality;
    /// ToT counts per microsecond
    double m_countsPerMicrosecond;
    /// array of gains, in microseconds/fC
    float m_ToTGain      [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quadratic terms, in microseconds/fC**2
    float m_ToTGain2     [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Thresholds, in microseconds = extrapolation to zero charge
    float m_ToTThreshold [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quality factors
    float m_ToTQuality   [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// pointer to geometry service
    ITkrGeometrySvc* m_geoSvc;
};


#endif // TkrToTSvc_H
