/** @file ITkrGeometrySvc.h
 @brief Abstract interface to TkrGeometrySvc (q.v.)

  $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrGeometrySvc.h,v 1.27 2006/11/02 19:34:47 lsrea Exp $
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
#include "TkrUtil/ITkrToTSvc.h"

#include "CLHEP/Geometry/Point3D.h"

#include <string>

/** 
 * @class ITkrGeometrySvc
 *
 * @brief Abstract interface to TkrGeometrySvc (q.v.)
 * 
 * @author Tracy Usher
 */

static const InterfaceID IID_ITkrGeometrySvc("ITkrGeometrySvc", 14 , 0); 

namespace {
    enum convType { ABSENT = -1, NOCONV = 0, STANDARD, SUPER, ALL, NCONVTYPES};
    convType ANYCONV = ALL;
    enum limitType { LOW, HIGH };
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
    virtual int    numTrays()       const = 0;
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
    virtual double siActiveWaferSide() const = 0;
    virtual double ladderPitch()    const = 0;
    virtual double waferPitch()     const = 0;
    virtual int    chipsPerLadder() const = 0;
    virtual int    stripsPerChip () const = 0;

    
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
    virtual int reversePlaneNumber(int plane) const = 0;

    /// Return the strip position (in local coordinates) given the stripId
    virtual HepPoint3D getStripPosition( int tower, int layer, int view, 
        double stripId) const = 0;

    /// Return z position for layer and view
    virtual double getReconLayerZ(int layer, int view=2) const = 0;

        /// return the z position for a digiLayer and view
    virtual double getLayerZ     (int digiLayer,  int view=2) const = 0;
    /// returns Z of *Layer* (average of x and y plane)
    /// TkrId of either plane will work
    virtual double getLayerZ     (const idents::TkrId& tkrId) const = 0;

    /// new stuff, based on plane and TkrId;
    virtual int    getPlane (const idents::TkrId& tkrId) const = 0;
    virtual int getPlane(double z) const = 0;

    virtual double getPlaneZ(int plane) const = 0;
    virtual double getPlaneZ(const idents::TkrId& tkrId) const = 0;
    virtual int    getLayer (int plane) const = 0;
    virtual int    getLayer (const idents::TkrId& tkrId) const = 0;
    virtual int    getView  (int plane) const = 0;
    virtual int    getView  (const idents::TkrId& tkrId) const = 0;

    /// Return radlen for the converter in a layer
    virtual double getReconRadLenConv(int layer) const = 0;
    virtual double getRadLenConv(int layer) const = 0;
    /// Return radlen for the rest of the layer
    ///    counting down from the bottom of the converter
    virtual double getReconRadLenRest(int layer) const = 0;
    virtual double getRadLenRest(int layer) const = 0;
    /// Return converter type for a layer
    virtual convType getReconLayerType(int layer) const = 0;
    virtual convType getLayerType(int layer) const = 0;
    
    virtual int getNumType(convType type) const = 0;
    /// get average radlen of converter for each type
    virtual double getAveConv(convType type) const = 0;
    /// get average radlen of rest for each type)
    ///    counting down from the bottom of the converter
    virtual double getAveRest(convType type) const = 0;

    /// Provide access to the old propagator
    virtual IKalmanParticle*    getPropagator()        const = 0;
    /// Provide access to the new propagator
    virtual IPropagator*        getG4PropagationTool() const = 0;
    /// Provide access to the failure mode service
    virtual ITkrFailureModeSvc* getTkrFailureModeSvc() = 0;
    /// Provide acess to the alignment service
    virtual ITkrAlignmentSvc*   getTkrAlignmentSvc()   = 0;
    /// Provide access to the bad strips service
    virtual ITkrBadStripsSvc*   getTkrBadStripsSvc()   = 0;
    /// Provide access to the splits service
    virtual ITkrSplitsSvc*      getTkrSplitsSvc()      = 0;
    /// Provide access to the ToT service
    virtual ITkrToTSvc*         getTkrToTSvc()         = 0;

    /// calculate the tray number, botTop from layer, view
    virtual void layerToTray (int layer, int view, int& tray, int& botTop) const = 0;
    /// calculate layer, view from tray, botTop
    virtual void trayToLayer (int tray, int botTop, int& layer, int& view) const = 0;    
    /// calculate layer (digi format) and view from plane
    virtual void planeToLayer (int plane, int& layer, int& view) const = 0;

    /// does the tower exist?
    virtual bool isTower(int tower)      const = 0;
    /// get tower type
    virtual int  getTowerType(int tower) const = 0;
    /// get lowest and highest tower in either direction
    virtual int  getLimitingTower(int view, enum limitType) const = 0;
    /// get the actual limits in either direction
    virtual double getLATLimit   (int view, enum limitType) const = 0;
    /// are we in the "active" LAT?
    virtual bool   isInActiveLAT (Point pos) const = 0;

    // definitions of plane, layer
    virtual int trayToPlane(int tray, int botTop) const = 0; 
    virtual int trayToBiLayer(int tray, int botTop) const = 0;
    virtual int planeToTray(int plane) const = 0;
    virtual int planeToBotTop(int plane) const = 0;
    virtual int getBottomTrayFlag() const = 0;
    virtual int getTopTrayFlag()    const = 0;
    virtual unsigned int getDefaultClusterStatus() const = 0;

    // changes added at end for minimal disruption, should reorganize later
    virtual int getPlaneSeparation(
        const idents::TkrId& id1, const idents::TkrId& id2) const = 0;

    virtual double truncateCoord( double x, double pitch, 
        int numElements, int& elementNumber, bool reverse = false) const = 0;
    virtual double getConvZ(int layer) const = 0;
    virtual bool isTopPlaneInLayer(int plane) const = 0;
    virtual double gettkrZBot() const = 0;

    // new stuff here for minimal disruption
    // put this back... Johann uses it!
    virtual bool inTower(int view, const Point p, int& iXTower, int& iYTower,
        double& xActiveDist, double& yActiveDist, double& xGap, double &yGap) const = 0;

};

#endif
