/**
 * @class TkrTrkParams
 *
 * @brief Implementation of a Vector for the generic Kalman Filter. This version based on CLHEP HepVector
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/TkrTrkParams.h,v 1.1 2004/09/18 18:38:42 usher Exp $
 */

#ifndef TkrTrkParams_h
#define TkrTrkParams_h

#include "CLHEP/Matrix/Vector.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"


class TkrTrkParams : public HepVector, virtual public Event::ITkrTrackParamsAccess
{
public:

    // Destructor
    virtual ~TkrTrkParams() {}

    // Constructors from HepMatrix
    TkrTrkParams() : HepVector() {}
    // Default constructor. Gives vector of length 0.
    // Another Vector can be assigned to it.

    TkrTrkParams(int p) : HepVector(p) {}
    TkrTrkParams(int p, int i) : HepVector(p, i) {}

    TkrTrkParams(const TkrTrkParams &m1) : HepVector(m1) {}
    // Constructor. Gives vector of length p.

#ifdef HEP_USE_RANDOM
    TkrTrkParams(int p, HepRandom &r) : HepVector(p,r) {}
#endif

    TkrTrkParams(const HepVector &v) : HepVector(v) {}
    TkrTrkParams(const HepMatrix &m) : HepVector(m) {}

    inline TkrTrkParams(Event::TkrTrackParams& m1);

    inline void setParams(Event::TkrTrackParams* params);
    inline void getParams(Event::TkrTrackParams* params);
};

TkrTrkParams::TkrTrkParams(Event::TkrTrackParams& m1) : HepVector(4)
{
    (*this)(1) = m1(1);
    (*this)(2) = m1(2);
    (*this)(3) = m1(3);
    (*this)(4) = m1(4);
}

inline void TkrTrkParams::setParams(Event::TkrTrackParams* params)
{
    params->setxPosition((*this)(1));
    params->setxSlope((*this)(2));
    params->setyPosition((*this)(3));
    params->setySlope((*this)(4));

    return;
}

inline void TkrTrkParams::getParams(Event::TkrTrackParams* params)
{
    (*this)(1) = params->getxPosition();
    (*this)(2) = params->getxSlope();
    (*this)(3) = params->getyPosition();
    (*this)(4) = params->getySlope();

    return;
}
#endif

