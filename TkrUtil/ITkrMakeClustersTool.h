/** @file ITkrMakeClustersTool.h
* @brief Abstract interface for methods to make TkrClusters

 @author Leon Rochester

 $Header$
*/

#ifndef _H_ITkrMakeClustersTool
#define _H_ITkrMakeClustersTool

#include "GaudiKernel/IAlgTool.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Digi/TkrDigi.h"
#include "TkrUtil/ITkrBadStripsSvc.h"


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrMakeClustersTool("ITkrMakeClustersTool", 0 , 0); 

/** @class ITkrMakeClustersTool
* @brief Abstract interface for methods to query TkrClusters
*/

class   ITkrMakeClustersTool : virtual public IAlgTool {
public:
 
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrMakeClustersTool; }

    virtual StatusCode makeClusters(
        Event::TkrClusterCol* pClus, 
        Event::TkrIdClusterMap* clusMap,
        Event::TkrDigiCol* pTkrDigiCol,
        ITkrBadStripsSvc::clusterType clType=ITkrBadStripsSvc::STANDARDCLUSTERS) = 0;

};

#endif  // _H_ITkrMakeClustersTool
