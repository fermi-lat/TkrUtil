/*
@file TkrToTSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrUtil/src/TkrToTSvc.h,v 1.14 2006/11/02 19:34:48 lsrea Exp $

*/
#ifndef TkrToTSvc_H
#define TkrToTSvc_H 1

// Include files
#include "TkrUtil/ITkrToTSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"
#include "idents/TowerId.h"
#include "idents/TkrId.h"

/** @class TkrToTSvc
* @brief Service to store and compare to a list of desired failure modes in
* the TKR.
*
* Author:  L. Rochester (after R.Dubois)
*
*/

enum calibType {GAIN=0, QUAD, THRESHOLD, SCALE, QUALITY};

class TkrToTSvc : public Service, virtual public ITkrToTSvc  {

public:

    TkrToTSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode finalize();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrToTSvc::interfaceID(); 
    }

    /// return the service type
    const InterfaceID& type() const;

    double getGain(int tower, int layer, int view, int strip) const;
    double getQuad(int tower, int layer, int view, int strip) const; 
    double getThreshold(int tower, int layer, int view, int strip) const; 
    double getQuality(int tower, int layer, int view, int strip) const; 
    double getMuonScale(int tower, int layer, int view, int strip) const;
    void   doWarning(const std::string type, const idents::TkrId id, int string) const;

    double getCountsPerMicrosecond() const { return m_countsPerMicrosecond;}
    double getMevPerMip() const { return m_mevPerMip; }
    double getFCPerMip() const { return m_fCPerMip; }
    int    getMaxToT() const { return m_maxToT; }

    double getCharge(double rawToT, int tower, int layer, int view, int strip) const;
    double getMipsFromToT(double rawToT, int tower, int layer, int view, int strip) const;
    double getMipsFromCharge(double charge) const;
    int    getRawToT(double eDep, int tower, int layer, int view, int strip) const;

    /// update the pointer
    void update(CalibData::TkrTotCol* pToT) { m_pToT = pToT; }
    /// update the other pointer!
    void update(CalibData::TkrScaleCol* pScale) {m_pScale = pScale; }

private:
    /// internal init method
    StatusCode doInit();

    /// check index
    bool valid(int tower, int layer, int view, int strip) const
    {
        return (tower>-1 && tower <m_nTowers && layer>-1 && layer<m_nLayers
            && view>-1 && view<m_nViews && strip>-1 && strip<m_nStrips);
    }

    void getConsts(idents::TkrId id, int strip, 
        double& threshold, double& gain, double& quad, double& muonScale) const;

    idents::TkrId getId(int tower, int layer, int view) const;

    /// mode: currently "default" or "EM"
    std::string m_mode;
    /// default Gain
    double m_defaultGain;
    /// default quadratic term;
    double m_defaultQuad;
    /// default Threshold
    double m_defaultThreshold;
    /// default quality factor
    double m_defaultQuality;
    /// default muon correction factor
    double m_defaultMuonScale;
    /// ToT counts per microsecond
    double m_countsPerMicrosecond;
    /// Energy deposited by a mini particle traversing a silicon plane
    double m_mevPerMip;
    /// Charge deposited by a mini particle traversing a silicon plane
    double m_fCPerMip;
    int    m_maxToT;
    /// flag to use default values if not in calibration file (otherwise use zeros)
    bool   m_useDefaultIfMissing;

    // Philippe's correction
    double m_linCorr;
    double m_quadCorr;

    /// pointer to geometry service
    ITkrGeometrySvc* m_tkrGeom;
    /// pointer to ToT consts
    CalibData::TkrTotCol* m_pToT;
    /// pointer to muonScale consts
    CalibData::TkrScaleCol* m_pScale;
    /// flag for using one tower for all the consts
    bool m_useSingleTowerConsts;
    /// tower to duplicate
    int m_baseTower;

    /// some useful constants
    int m_nTowers;
    int m_nLayers;
    int m_nViews;
    int m_nStrips;

    mutable int m_callCount;
};


#endif // TkrToTSvc_H

