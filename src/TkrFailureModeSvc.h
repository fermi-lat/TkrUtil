#ifndef TkrFailureModeSvc_H
#define TkrFailureModeSvc_H 1

// Include files
#include "TkrUtil/ITkrFailureModeSvcCalib.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"

#include <vector>
#include <map>

/** @class TkrFailureModeSvc
* @brief Service to store and compare to a list of desired failure modes in
* the TKR.
*
* Author:  L. Rochester (after R.Dubois)
*
*/

enum {TOWER_SHIFT = 0, LAYER_SHIFT = 1};

namespace {
    /**
    @class BadVisitorFM

    Minimal class derived from CalibData::BadStripsVisitor to
    check out BadStrips visitor interface.
    */
    class BadVisitorFM : public CalibData::BadStripsVisitor {
    public:
        BadVisitorFM(MsgStream* log=0) : m_log(log){}

        void setLog(MsgStream* log) {m_log = log;}

        virtual CalibData::eVisitorRet badTower(unsigned int row, unsigned int col,
            int badness);

        virtual CalibData::eVisitorRet badPlane(unsigned int row, unsigned int col, 
            unsigned int tray, bool top,
            int badness, bool allBad,
            const CalibData::StripCol& strips);

        void setService(ITkrFailureModeSvcCalib* pFailureMode) {m_pFailureMode = pFailureMode;}
        void setService(ITkrGeometrySvc*         pGeoSvc)      {m_pGeoSvc      = pGeoSvc;}

    private:
        MsgStream* m_log;
        ITkrFailureModeSvcCalib* m_pFailureMode;
        ITkrGeometrySvc*         m_pGeoSvc;
    };
}

class TkrFailureModeSvc : public Service, virtual public ITkrFailureModeSvcCalib {

public:

    TkrFailureModeSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    //StatusCode execute();
    StatusCode finalize();

    /// get the list of enabled failure mode conditions
    int getFailureConditions() const {return m_failureModes;}

    bool empty() const { return m_failureModes==0;}
    void setFailureModes( int modeBit ) { m_failureModes = m_failureModes | modeBit;}


    StatusCode update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot);


    /// Find out if object is marked Failed
    bool isFailed(int towerId, int layer, int view) const;

    std::vector<int>& getLayers( int tower);
    std::vector<int>& getTowers ();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrFailureModeSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

private:

    /// look for tower in list of dead towers
    bool towerFailed(int tower) const;

    /// look for layer in list of dead layers
    bool layerFailed(int tower, int layer, int view) const;

    /// process the input list of towers
    void processTowerList();

    /// process the input list of layer
    void processLayerList();

    /// do common init for update
    StatusCode doInit();

private:

    /// List of towers from jobOptions
    StringArrayProperty m_towerListProperty;

    /// List of layers from jobOptions
    StringArrayProperty m_layerListProperty;

    /// bitmap of failure modes
    int m_failureModes;
    /// tells whether a list was read in
    bool m_existsList;

    /// vector of towers to fail
    std::vector<int> m_towerList;

    /// vector of layers to fail
    //    std::map <int, std::vector<int> > m_layerList;

    BadVisitorFM* m_visitor;


    /// vector of layers to fail
    LayerMap  m_layerList;

    ITkrGeometrySvc*     m_pGeoSvc;
    };


#endif // TkrFailureModeSvc_H

