/**
 * @class ITkrGhostTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrGhostTool.h,v 1.7 2010/04/08 20:54:04 lsrea Exp $
 */
#ifndef ITkrGhostTool_h
#define ITkrGhostTool_h

//namespace {
//    // set the fields in the indices
//    const int viewMult  = 1;
//    const int layerMult = 2;
//    const int towerMult = 1000;
//}

#include "GaudiKernel/IAlgTool.h"

#include "TkrUtil/TkrTowerBits.h"

static const InterfaceID IID_ITkrGhostTool("ITkrGhostTool", 3 , 0);
class ITkrGhostTool : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrGhostTool; }

    virtual StatusCode getTkrVector(unsigned short& tkrVector) = 0;
    virtual StatusCode calculateTkrVector(
        Event::TkrClusterCol* pCol, unsigned short& towerBits) = 0;
    virtual StatusCode calculateTkrVector(
        Event::TkrDigiCol* pCol, unsigned short& towerBits) = 0;
    virtual StatusCode flagSingles()   = 0;
    virtual StatusCode flagEarlyHits(Event::TkrClusterCol* col=0) = 0;
    virtual StatusCode flagEarlyTracks() = 0;
    virtual StatusCode flagEarlyVertices() = 0;
    virtual StatusCode flagEarlyCalClusters() = 0;

};

#endif
