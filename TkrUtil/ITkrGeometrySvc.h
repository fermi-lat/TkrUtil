/** @file ITkrGeometrySvc
@brief Supplies useful geometry information for TKR digitization and reconstruction

  $Header$
*/

#ifndef __ITKRGEOMETRYSVC_H
#define __ITKRGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"

#include "GlastSvc/Reco/IKalmanParticle.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"

#include "CLHEP/Geometry/Point3D.h"

#include <string>

/** 
 * @class ITkrGeometrySvc
 *
 * @brief Abstract interface to TkrGeometrySvc (q.v.)
 * 
 * @author Tracy Usher
 */

static const InterfaceID IID_ITkrGeometrySvc("ITkrGeometrySvc", 2 , 0); 

class ITkrGeometrySvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    static const InterfaceID& interfaceID() { return IID_ITkrGeometrySvc; }
    
    //Retrieve stored information

    virtual int    numXTowers()=0;
    virtual int    numYTowers()=0;
    virtual int    numViews()=0;
    virtual int    numLayers()=0;
    virtual int    numNoConverter()=0;
    virtual int    numSuperGlast()=0;

    virtual int    numPlanes()=0;

    virtual double towerPitch()=0;
    virtual double trayWidth()=0;
    virtual double trayHeight()=0;
    
    virtual double ladderGap()=0;
    virtual double ladderInnerGap()=0;
    virtual int    ladderNStrips()=0;
    
    virtual double siStripPitch()=0;
    virtual double siResolution()=0;
    virtual double siThickness()=0;
    virtual double siDeadDistance()=0;

    // Digi and Reco layers differ in the ordering
    virtual int ilayer(int layer)=0; // deprecated
    /// convert from digi<->recon layer number
    virtual int reverseLayerNumber(int layer)=0;

    /// Return the strip position (in local coordinates) given the stripId
    virtual HepPoint3D getStripPosition( int tower, int layer, int view, 
        double stripId) = 0;

    /// Return z position for layer and view
    virtual double getReconLayerZ(int layer, int view) = 0;
    /// Return average z position for a layer
    virtual double getReconLayerZ(int layer) = 0;

    /// Provide access to the propagator
    virtual IKalmanParticle* getPropagator() = 0;
    /// Provide access to the failure mode service
    virtual ITkrFailureModeSvc* getTkrFailureModeSvc() = 0;
    /// Provide acess to the alignment service
    virtual ITkrAlignmentSvc* getTkrAlignmentSvc() = 0;
    /// Provide access to the bad strips service
    virtual ITkrBadStripsSvc* getTkrBadStripsSvc() = 0;
};

#endif