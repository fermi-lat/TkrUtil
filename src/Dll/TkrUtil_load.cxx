/** 
* @file GlastDigi_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/Dll/TkrUtil_load.cxx,v 1.6 2003/01/27 00:37:46 lsrea Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(TkrUtil) {
    DECLARE_SERVICE( TkrFailureModeSvc   );
    DECLARE_SERVICE( TkrBadStripsSvc     );
    DECLARE_SERVICE( TkrGeometrySvc      );
    DECLARE_SERVICE( TkrAlignmentSvc     );

    DECLARE_TOOL(    TkrQueryClustersTool);
    DECLARE_TOOL(    TkrMeritTool        );

    DECLARE_ALGORITHM ( TkrCalibAlg      );
    //This is for test only, not part of TkrUtil package...
    //DECLARE_ALGORITHM ( EvtClock         );
} 



