/*
@file TkrSplitsSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrSplitsSvc.h,v 1.1 2004/03/10 18:35:03 lsrea Exp $

*/
#ifndef TkrSplitsSvc_H
//#define TkrSplitsSvc_H 1

// Include files
#include "TkrUtil/ITkrSplitsSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"

/** @class TkrSplitsSvc
* @brief Service to store and compare to a list of desired failure modes in
* the TKR.
*
* Author:  L. Rochester (after R.Dubois)
*
*/

class TkrSplitsSvc : public Service, virtual public ITkrSplitsSvc  {

public:

    enum {NTOWERS=16, NLAYERS=18, NVIEWS=2};

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

    /// set the pointer to the geometry
    ///void setTkrGeometrySvc(TkrGeometrySvc* p_geoSvc) { m_geoSvc = p_geoSvc; }

    /// get the last C0 strip for this layer
    int getSplitPoint(int tower, int layer, int view) const {
        return m_splits[tower][layer][view];
    }

    /// tell which end this strip belongs to
    int getEnd(int tower, int layer, int view, int strip) const {
        return (strip<=m_splits[tower][layer][view]? 0 : 1);
    }

private:
    /// internal init method
    StatusCode doInit();
    /// name of file containing splits
    std::string m_splitsFile;
    /// array of splits
    int m_splits[NTOWERS][NLAYERS][NVIEWS];
    /// pointer to the geometry
    ITkrGeometrySvc* m_geoSvc;

};


#endif // TkrSplitsSvc_H

