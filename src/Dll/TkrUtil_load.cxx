/** 
* @file GlastDigi_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/Dll/TkrUtil_load.cxx,v 1.2 2003/01/06 18:22:47 lsrea Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();


DECLARE_FACTORY_ENTRIES(TkrUtil) {
    DECLARE_SERVICE( TkrFailureModeSvc );
    DECLARE_SERVICE( TkrBadStripsSvc   );
    DECLARE_SERVICE( TkrGeometrySvc    );
} 



