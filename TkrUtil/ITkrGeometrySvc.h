/** @file ITkrGeometrySvc.h
 @brief Abstract interface to TkrGeometrySvc (q.v.)

  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrGeometrySvc.h,v 1.9 2003/07/18 22:27:19 lsrea Exp $
*/

#ifndef __ITKRGEOMETRYSVC_H
#define __ITKRGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"

#include "GlastSvc/Reco/IKalmanParticle.h"
#include "GlastSvc/Reco/IPropagator.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"

#include "CLHEP/Geometry/Point3D.h"

#include <string>

/** 
 * @class ITkrGeometrySvc
 *
 * @brief Abstract interface to TkrGeometrySvc (q.v.)
 * 
 * @author Tracy Usher
 */

static const InterfaceID IID_ITkrGeometrySvc("ITkrGeometrySvc", 7 , 0); 

namespace {
    enum convType { ABSENT = -1, NOCONV = 0, STANDARD, SUPER, ALL, NTYPES};
}

class ITkrSplitsSvc;

class ITkrGeometrySvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    static const InterfaceID& interfaceID() { return IID_ITkrGeometrySvc; }
    
    //Retrieve stored information

    virtual int    numXTowers()     const = 0;
    virtual int    numYTowers()     const = 0;
    virtual int    numViews()       const = 0;
    virtual int    numLayers()      const = 0;
    virtual int    numNoConverter() const = 0;
    virtual int    numSuperGlast()  const = 0;
    virtual int    numRegular()     const = 0;

    virtual int    numPlanes()      const = 0;

    virtual double towerPitch()     const = 0;
    virtual double trayWidth()      const = 0;
    virtual double trayHeight()     const = 0;
    
    virtual double ladderGap()      const = 0;
    virtual double ladderInnerGap() const = 0;
    virtual int    ladderNStrips()  const = 0;
    virtual int    nWaferAcross()   const = 0;
    virtual double siWaferSide()    const = 0;
    
    virtual double siStripPitch()   const = 0;
    virtual double siResolution()   const = 0;
    virtual double siThickness()    const = 0;
    virtual double siDeadDistance() const = 0;

    virtual double calZTop()        const = 0;
	virtual double calZBot()        const = 0;
	virtual double calXWidth()      const = 0;
	virtual double calYWidth()      const = 0;


    // Digi and Reco layers differ in the ordering
    /// convert from digi<->recon layer number
    virtual int reverseLayerNumber(int layer) const = 0;

    /// Return the strip position (in local coordinates) given the stripId
    virtual HepPoint3D getStripPosition( int tower, int layer, int view, 
        double stripId) const = 0;

    /// Return z position for layer and view
    virtual double getReconLayerZ(int layer, int view=2) const = 0;

    /// Return radlen for the converter in a layer
    virtual double getReconRadLenConv(int layer) const = 0;
    /// Return radlen for the rest of the layer
    ///    counting down from the bottom of the converter
    virtual double getReconRadLenRest(int layer) const = 0;
    /// Return converter type for a layer
    virtual convType getReconLayerType(int layer) const = 0;
    
    virtual int getNumType(convType type) const = 0;
    /// get average radlen of converter for each type
    virtual double getAveConv(convType type) const = 0;
    /// get average radlen of rest for each type)
    ///    counting down from the bottom of the converter
    virtual double getAveRest(convType type) const = 0;

    /// Provide access to the old propagator
    virtual IKalmanParticle*    getPropagator() const = 0;
    /// Provide access to the new propagator
    virtual IPropagator*        getG4PropagationTool() const = 0;
    /// Provide access to the failure mode service
    virtual ITkrFailureModeSvc* getTkrFailureModeSvc() const = 0;
    /// Provide acess to the alignment service
    virtual ITkrAlignmentSvc*   getTkrAlignmentSvc() const = 0;
    /// Provide access to the bad strips service
    virtual ITkrBadStripsSvc*   getTkrBadStripsSvc() const = 0;
    /// Provide access to the splitss service
    virtual ITkrSplitsSvc*      getTkrSplitsSvc() const = 0;


    /// calculate the tray number, botTop from layer, view
    virtual void layerToTray (int layer, int view, int& tray, int& botTop) const = 0;
    /// calculate layer, view from tray, botTop
    virtual void trayToLayer (int tray, int botTop, int& layer, int& view) const = 0;    
    /// calculate layer (digi format) and view from plane
    virtual void planeToLayer (int plane, int& layer, int& view) const = 0;
};

#endif
