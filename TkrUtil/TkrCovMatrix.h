/**
 * @class TkrCovMatrix
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/TkrCovMatrix.h,v 1.1 2004/09/18 18:38:42 usher Exp $
 */

#ifndef TkrCovMatrix_h
#define TkrCovMatrix_h

#include "CLHEP/Matrix/Matrix.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class TkrCovMatrix : public HepMatrix, virtual public Event::ITkrTrackParamsAccess
{
public:

    // Destructor
    virtual ~TkrCovMatrix() {}

    // Constructors from HepMatrix
    TkrCovMatrix() : HepMatrix() {}
    TkrCovMatrix(int p, int q): HepMatrix(p,q) {}
    TkrCovMatrix(int p, int q, int i) : HepMatrix(p,q,i) {}
    TkrCovMatrix(int p, int q, HepRandom &r) : HepMatrix(p,q,r) {}
    TkrCovMatrix(const HepMatrix &m1) : HepMatrix(m1) {}
    TkrCovMatrix(const TkrCovMatrix &m1) : HepMatrix(m1) {}
    // Copy constructor.

    TkrCovMatrix(const HepSymMatrix &m1) : HepMatrix(m1) {}
    TkrCovMatrix(const HepDiagMatrix &m1) : HepMatrix(m1) {}
    TkrCovMatrix(const HepVector &m1) : HepMatrix(m1) {}
    // Constructors from SymMatrix, DiagMatrix and Vector.

    inline TkrCovMatrix(Event::TkrTrackParams& m1);

    inline void setParams(Event::TkrTrackParams* params);
    inline void getParams(Event::TkrTrackParams* params);
};

TkrCovMatrix::TkrCovMatrix(Event::TkrTrackParams& m1) : HepMatrix(4,4) 
{
    (*this)(1,1) = m1(1,1);
    (*this)(1,2) = m1(1,2);
    (*this)(1,3) = m1(1,3);
    (*this)(1,4) = m1(1,4);
    (*this)(2,1) = m1(2,1);
    (*this)(2,2) = m1(2,2);
    (*this)(2,3) = m1(2,3);
    (*this)(2,4) = m1(2,4);
    (*this)(3,1) = m1(3,1);
    (*this)(3,2) = m1(3,2);
    (*this)(3,3) = m1(3,3);
    (*this)(3,4) = m1(3,4);
    (*this)(4,1) = m1(4,1);
    (*this)(4,2) = m1(4,2);
    (*this)(4,3) = m1(4,3);
    (*this)(4,4) = m1(4,4);
}

inline void TkrCovMatrix::setParams(Event::TkrTrackParams* params)
{
    params->setxPosxPos((*this)(1,1));
    params->setxPosxSlp((*this)(1,2));
    params->setxPosyPos((*this)(1,3));
    params->setxPosySlp((*this)(1,4));
    params->setxSlpxSlp((*this)(2,2));
    params->setxSlpyPos((*this)(2,3));
    params->setxSlpySlp((*this)(2,4));
    params->setyPosyPos((*this)(3,3));
    params->setyPosySlp((*this)(3,4));
    params->setySlpySlp((*this)(4,4));

    return;
}

inline void TkrCovMatrix::getParams(Event::TkrTrackParams* params)
{
    (*this)(1,1) = params->getxPosxPos();
    (*this)(1,2) = params->getxPosxSlp();
    (*this)(1,3) = params->getxPosyPos();
    (*this)(1,4) = params->getxPosySlp();
    (*this)(2,2) = params->getxSlpxSlp();
    (*this)(2,3) = params->getxSlpyPos();
    (*this)(2,4) = params->getxSlpySlp();
    (*this)(3,3) = params->getyPosyPos();
    (*this)(3,4) = params->getyPosySlp();
    (*this)(4,4) = params->getySlpySlp();

    (*this)(2,1) = (*this)(1,2);
    (*this)(3,1) = (*this)(1,3);
    (*this)(3,2) = (*this)(2,3);
    (*this)(4,1) = (*this)(1,4);
    (*this)(4,2) = (*this)(2,4);
    (*this)(4,3) = (*this)(3,4);

    return;
}

#endif

