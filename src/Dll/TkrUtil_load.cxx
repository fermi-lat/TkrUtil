/** 
* @file TkrUtil_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/Dll/TkrUtil_load.cxx,v 1.13 2004/09/12 03:53:07 lsrea Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(TkrUtil) {
    DECLARE_SERVICE( TkrFailureModeSvc   );
    DECLARE_SERVICE( TkrBadStripsSvc     );
    DECLARE_SERVICE( TkrGeometrySvc      );
    DECLARE_SERVICE( TkrAlignmentSvc     );
    DECLARE_SERVICE( TkrSplitsSvc        );
    DECLARE_SERVICE( TkrToTSvc           );

    DECLARE_TOOL(    TkrQueryClustersTool);
    DECLARE_TOOL(    TkrMakeClustersTool );

    DECLARE_ALGORITHM ( TkrCalibAlg      );
 } 



