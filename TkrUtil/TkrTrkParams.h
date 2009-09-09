/**
 * @class TkrTrkParams
 *
 * @brief Implementation of a Vector for the generic Kalman Filter. This version based on CLHEP CLHEP::HepVector
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/TkrTrkParams.h,v 1.3 2006/03/21 01:15:47 usher Exp $
 */

#ifndef TkrTrkParams_h
#define TkrTrkParams_h

#include "CLHEP/Matrix/Vector.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"


class TkrTrkParams : public CLHEP::HepVector, virtual public Event::ITkrTrackParamsAccess
{
public:

    // Destructor
    virtual ~TkrTrkParams() {}

    // Constructors from CLHEP::HepMatrix
    TkrTrkParams() : CLHEP::HepVector() {}
    // Default constructor. Gives vector of length 0.
    // Another Vector can be assigned to it.

    TkrTrkParams(int p) : CLHEP::HepVector(p) {}
    TkrTrkParams(int p, int i) : CLHEP::HepVector(p, i) {}

    TkrTrkParams(const TkrTrkParams &m1) : ITkrTrackParamsAccess(), CLHEP::HepVector(m1) {}
    // Constructor. Gives vector of length p.

#ifdef HEP_USE_RANDOM
    TkrTrkParams(int p, CLHEP::HepRandom &r) : CLHEP::HepVector(p,r) {}
#endif

    TkrTrkParams(const CLHEP::HepVector &v) : CLHEP::HepVector(v) {}
    TkrTrkParams(const CLHEP::HepMatrix &m) : CLHEP::HepVector(m) {}

    inline TkrTrkParams(Event::TkrTrackParams& m1);

    inline void setParams(Event::TkrTrackParams* params);
    inline void getParams(Event::TkrTrackParams* params);
};

TkrTrkParams::TkrTrkParams(Event::TkrTrackParams& m1) : CLHEP::HepVector(4)
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

