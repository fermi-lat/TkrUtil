/**
 * @class ITkrMcTracksTool
 *
 * @brief Interface to the tool for returning information from the tables relating Monte Carlo 
 *        McParticles, McPositionHits and TkrRecon TkrClusters, which give information about the
 *        Monte Carlo tracks in the Tracker. 
 *
 * @author Tracy Usher
 */

#ifndef ITkrMcTracksTOOL_H
#define ITkrMcTracksTOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McSiLayerHit.h"

static const InterfaceID IID_ITkrMcTracksTool("ITkrMcTracksTool", 2 , 0);

class ITkrMcTracksTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrMcTracksTool; }

    /// @brief Return the number of Monte Carlo tracks
    virtual int                         getNumMcTracks()=0;

    /// @brief Returns information about the event
    virtual const unsigned long         getClassificationBits()=0;

    /// @brief Returns primary McParticle
    virtual const Event::McParticleRef  getPrimaryParticle()=0;

    /// @brief Returns secondary particles
    virtual int                         getNumSecondaries()=0;
    virtual const Event::McParticleRef  getSecondary(int mcPartIdx)=0;

    /// @brief Returns associated particles
    virtual int                         getNumAssociated()=0;
    virtual const Event::McParticleRef  getAssociated(int mcPartIdx)=0;

    /// @brief Returns a vector of hits associated as one McParticle track
    virtual const Event::McPartToHitVec getMcPartTrack(const Event::McParticleRef mcPart)=0;

    /// @brief Returns the layer number of a given McPositionHit on the track
    virtual const int                   getTrackHitLayer(const Event::McParticleRef mcPart, int hitIdx)=0;

    /// @brief Returns number of Tracker (cluster) hits for a given track
    virtual const int                   getNumClusterHits(const Event::McParticleRef mcPart)=0;

    /// @brief Returns number of shared Tracker (cluster) hits for a given track
    virtual const int                   getNumSharedTrackHits(const Event::McParticleRef mcPart)=0;

    /// @brief Returns number of gaps and their sizes (in layers) for a given track
    virtual const int                   getNumGaps(const Event::McParticleRef mcPart)=0;
    virtual const int                   getGapSize(const Event::McParticleRef mcPart, int gapIdx)=0;
    virtual const int                   getGapStartHitNo(const Event::McParticleRef mcPart, int gapIdx)=0;

    /// @brief Returns the "straightness" of a given track
    virtual const double                getTrackStraightness(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40)=0;
    virtual const Hep3Vector            getTrackDirection(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40)=0;

    /// @brief Returns energy loss information within the tracker volume
    virtual const double                getTrackTotEneLoss(const Event::McParticleRef mcPart)= 0;
    virtual const double                getTrackELastHit(const Event::McParticleRef mcPart)= 0;
    virtual const double                getTrackBremELoss(const Event::McParticleRef mcPart, int& nTotRad)= 0;
    virtual const double                getTrackDeltaELoss(const Event::McParticleRef mcPart, int& nTotDlta, int& nHitDlta)=0;
    virtual const int                   getTrackDeltaRange(const Event::McParticleRef mcPart, double& aveRange, double& maxRange)=0;

    /// @brief Compares two tracks and returns information on shared hits (if any)
    virtual const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart)=0;
    virtual const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart1, const Event::McParticleRef mcPart2)=0;

};
#endif
