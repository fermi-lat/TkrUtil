/*
@file TkrToTSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrToTSvc.h,v 1.4 2004/05/10 23:58:51 lsrea Exp $

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
    double getMevPerMip() const { return m_mevPerMip; }
    double getFCPerMip() const { return m_fCPerMip; }
    int    getMaxToT() const { return m_maxToT; }

    double getCharge(double ToT, int tower, int layer, int view, int strip) const;
    double getMipsFromToT(double ToT, int tower, int layer, int view, int strip) const;
    double getMipsFromCharge(double charge, int tower, int layer, int view, int strip) const;

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
    /// default muon correction factor
    double m_defaultMuonFactor;
    /// ToT counts per microsecond
    double m_countsPerMicrosecond;
    /// Energy deposited by a mini particle traversing a silicon plane
    double m_mevPerMip;
    /// Charge deposited by a mini particle traversing a silicon plane
    double m_fCPerMip;
    /// array of gains, in microseconds/fC
    int    m_maxToT;
    float m_ToTGain      [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quadratic terms, in microseconds/fC**2
    float m_ToTGain2     [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Thresholds, in microseconds = extrapolation to zero charge
    float m_ToTThreshold [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quality factors
    float m_ToTQuality   [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Muon correction normalizations, should be order(1)
    float m_ToTMuonFactor [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// pointer to geometry service
    ITkrGeometrySvc* m_tkrGeom;
};


#endif // TkrToTSvc_H

