#ifndef TkrFailureModeSvc_H
#define TkrFailureModeSvc_H 1

// Include files
#include "TkrUtil/ITkrFailureModeSvc.h"
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

class TkrFailureModeSvc : public Service, virtual public ITkrFailureModeSvc {

public:

    TkrFailureModeSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    /// get the list of enabled failure mode conditions
    int getFailureConditions() {return m_failureModes;};

    enum {TOWER, LAYER};

    /// Find out if object is marked Failed
    bool isFailed(int towerId, int layer, int view);

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrFailureModeSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

private:

    /// look for tower in list of dead towers
    bool towerFailed(int tower);

    /// look for layer in list of dead layers
    bool layerFailed(int tower, int layer, int view);

    /// process the input list of towers
    void processTowerList();

    /// process the input list of layer
    void processLayerList();

private:

    /// List of towers from jobOptions
    StringArrayProperty m_towerListProperty;

    /// List of layers from jobOptions
    StringArrayProperty m_layerListProperty;

    /// bitmap of failure modes
    int m_failureModes;

    /// vector of towers to fail
    std::vector<int> m_towerList;

    /// vector of layers to fail
    std::map <int, std::vector<int> > m_layerList;

};

#endif // TkrFailureModeSvc_H

