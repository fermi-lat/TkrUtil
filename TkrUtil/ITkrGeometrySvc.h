/** @file ITkrGeometrySvc
@brief Supplies useful geometry information for TKR digitization and reconstruction

  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrGeometrySvc.h,v 1.3 2003/03/13 19:06:06 lsrea Exp $
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

static const InterfaceID IID_ITkrGeometrySvc("ITkrGeometrySvc", 3 , 0); 

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
    virtual int    nWaferAcross()=0;
    
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

    /// calculate the tray number, botTop from layer, view
    virtual void layerToTray (int layer, int view, int& tray, int& botTop) = 0;
    /// calculate layer, view from tray, botTop
    virtual void trayToLayer (int tray, int botTop, int& layer, int& view) = 0;    
    /// calculate layer (digi format) and view from plane
    virtual void planeToLayer (int plane, int& layer, int& view)=0;
    


};

#endif