/** 
* @file TkrUtil_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/Dll/TkrUtil_load.cxx,v 1.12 2004/03/13 19:40:37 lsrea Exp $
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
    //DECLARE_TOOL(    TkrMeritTool        ); gone...

    DECLARE_ALGORITHM ( TkrCalibAlg      );
    //This is for test only, not part of TkrUtil package...
    //DECLARE_ALGORITHM ( EvtClock         );
} 



