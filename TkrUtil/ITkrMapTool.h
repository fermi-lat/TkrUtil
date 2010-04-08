/**
 * @class ITkrMapTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrMapTool.h,v 1.1 2009/10/16 18:47:30 lsrea Exp $
 */
#ifndef ITkrMapTool_h
#define ITkrMapTool_h

namespace {
 
    // set the fields in the indices
    const int planeMult = 2;
    const int geoMult   = 2;
    const int gtccMult  = 10;

    const int nRc       = 9;
    const int nCc       = 8;

    // assign the ends to the 8 cables
    const int endArray[nCc] = {0, 1, 1, 0, 1, 0, 0, 1};
}

#include "GaudiKernel/IAlgTool.h"


static const InterfaceID IID_ITkrMapTool("ITkrMapTool", 1 , 0);

class ITkrMapTool : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrMapTool; }

    virtual int elecToGeo(int gtcc, int gtrc)  = 0;
    virtual int geoToElec(int plane, int end)  = 0;

    virtual void geoToElec(int plane, int end, int& gtcc, int& gtrc) = 0;
    virtual void elecToGeo(int gtcc, int gtrc, int& plane, int& end) = 0;

private:

};

#endif
