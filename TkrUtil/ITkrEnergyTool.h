/**
 * @class ITkrEnergyTool
 *
 * @brief Implements an interface for a Gaudi Tool for implementing code which will augment
 *        energy from the Calorimeter with that "seen" in the tracker
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrEnergyTool.h,v 1.1 2005/06/22 21:40:58 usher Exp $
 */
#ifndef ITkrEnergyTool_h
#define ITkrEnergyTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

static const InterfaceID IID_ITkrEnergyTool("ITkrEnergyTool", 1 , 0);

class ITkrEnergyTool : virtual public IAlgTool
{
public:

    /// @brief Defines the method to determine total event energy
    virtual double getTotalEnergy(const Event::TkrTrack* track, double calEnergy) = 0;

    virtual double getEvtEnergyEstimation(const Event::TkrTrack* track) = 0;

    virtual double GetEnergyUB2Correction(int method, int tkr1firstlayer, double tkr1zdir, double energy) = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrEnergyTool; }
};

#endif
