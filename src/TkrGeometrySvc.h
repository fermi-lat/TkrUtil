
#ifndef TKRGEOMETRYSVC_H
#define TKRGEOMETRYSVC_H 1

/** 
 * @class TkrGeometrySvc
 *
 * @brief Supplies the geometry constants and calculations to the 
 * TkrRecon Package
 *
 * The constants all flow from GlastDetSvc.
 * 
 * @author Leon Rochester
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrGeometrySvc.h,v 1.11 2004/03/12 05:49:22 lsrea Exp $
 */

#include "GaudiKernel/Service.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"


class TkrGeometrySvc : public Service,
        virtual public ITkrGeometrySvc
{
public:
    enum { NVIEWS=2, NLAYERS=18, NTOWERS=16};
    
    TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    //Retrieve stored information

    int    numXTowers() const     {return m_numX;} 
    int    numYTowers() const     {return m_numY;} 
    int    numViews()   const     {return m_nviews;} 
    int    numLayers()  const     {return getNumType(ALL);}
    int    numNoConverter() const {return getNumType(NOCONV);}
    int    numSuperGlast()  const {return getNumType(SUPER);}
    int    numRegular()     const {return getNumType(STANDARD);}
    int    nWaferAcross()   const {return m_nWaferAcross;}
    double siWaferSide()    const {return m_siWaferSide;}

    int    numPlanes()      const {return numLayers();}

    double towerPitch()     const {return m_towerPitch;}
    double trayWidth()      const {return m_trayWidth;}
    double trayHeight()     const {return m_trayHeight;}
    
    double ladderGap()      const {return m_ladderGap;}
    double ladderInnerGap() const {return m_ladderInnerGap;}
    int    ladderNStrips()  const {return m_ladderNStrips;}
    
    double siStripPitch()   const {return m_siStripPitch;}
    double siResolution()   const {return m_siResolution;}
    double siThickness()    const {return m_siThickness;}
    double siDeadDistance() const {return m_siDeadDistance;}

	double calZTop()        const {return m_calZTop;}
	double calZBot()        const {return m_calZBot;}
	double calXWidth()      const {return m_calXWidth;}
	double calYWidth()      const {return m_calYWidth;}

    /// reverse the numbering of the bilayers (goes both ways)
    int reverseLayerNumber(int layer) const {return numLayers()-layer-1;}

    /// return the position of a strip, will accept int or double
    HepPoint3D getStripPosition(int tower, int layer, int view, double stripid) const;

    /// return the z position for a reconLayer and view
    double getReconLayerZ(int layer, int view=2) const;
    /// return the average z position for a reconLayer
    //double getReconLayerZ(int layer);
    /// return the rad length of the converter for each layer
    double getReconRadLenConv(int layer) const { 
        return m_radLenConv[reverseLayerNumber(layer)];
    }
    /// return the rad length of the rest of the layer, for each layer
    ///   counting down from the bottom of the converter
    double getReconRadLenRest(int layer) const { 
        return m_radLenRest[reverseLayerNumber(layer)];
    }

    /// return converter type for a layer
    convType getReconLayerType(int layer) const;

    /// return number of layers of each type
    int getNumType(convType type) const { return m_numLayers[(int)type];}
    /// get average radlen of converter for each type
    double getAveConv(convType type) const { return m_aveRadLenConv[(int)type];}
    /// get average radlen of rest for each type)
    ///    counting down from the bottom of the converter
    double getAveRest(convType type) const {return m_aveRadLenRest[(int) type];}


    /// Provide access to the old propagator
    IKalmanParticle* getPropagator()     const {return m_KalParticle;}
    /// provide access to the new propagator
    IPropagator* getG4PropagationTool()  const {return m_G4PropTool;}

    /// Provide access to the failure mode service
    ITkrFailureModeSvc* getTkrFailureModeSvc() const { return m_tkrFail;}
    /// Provide access to the bad strips service
    ITkrBadStripsSvc*   getTkrBadStripsSvc()   const { return m_badStrips;}
    /// Provide access to the alignment service
    ITkrAlignmentSvc*   getTkrAlignmentSvc()   const { return m_tkrAlign;}
    /// Provide access to splits service
    ITkrSplitsSvc*      getTkrSplitsSvc()      const { return m_tkrSplits;}
    /// Provide access to ToT service
    ITkrToTSvc*         getTkrToTSvc()         const { return m_tkrToT;}

    /// calculate the tray number, botTop from layer, view
    void layerToTray (int layer, int view, int& tray, int& botTop) const;
    /// calculate layer, view from tray, botTop
    void trayToLayer (int tray, int botTop, int& layer, int& view) const;
    /// calculate layer (digi format) and view from plane
    void planeToLayer (int plane, int& layer, int& view) const;
    
 
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrGeometrySvc::interfaceID(); 
    }
    /// return the service type
    const IID& type() const;
    
private:
    
    /// number of Towers in X
    int    m_numX; 
    /// number of Towers in Y
    int    m_numY;              
    /// two views, always!
    int    m_nviews;        
    int m_nWaferAcross;

    /// Distance between centers of adjacent towers
    double m_towerPitch;
    /// Width of the ladders and the gaps between them
    double m_trayWidth;
    /// from top of one tray to the next (actually pitch) 
    /// -- Not really a constant, uses the smallest
    double m_trayHeight;   
    
    /// gap between adjacent ladders   
    double m_ladderGap;     
    /// gap between SSDs on the same ladder
    double m_ladderInnerGap;
    /// number of strips in a ladder
    int    m_ladderNStrips; 
    
    /// distance between two strips
    double m_siStripPitch; 
    /// nominally, (strip pitch)/sqrt(12)
    double m_siResolution;
    /// thickness of the silicon
    double m_siThickness;
    /// width of the dead region around the edge of a wafer
    double m_siDeadDistance;
    /// size of Wafer
    double m_siWaferSide;
    /// z coordinate of top of the top CAL crystal
	double m_calZTop;
	/// z coordinate of bottom of the bottom CAL crystal
	double m_calZBot;
	/// (maximum) width in x of active CsI across the entire instrument
	double m_calXWidth;
	/// (maximum) width in y of active CsI across the entire instrument
	double m_calYWidth;
    /// z positions of all the layers (digi convention)
    double m_layerZ[NLAYERS][NVIEWS];

    /// radiation lengths of converter, by recon layer
    double m_radLenConv[NLAYERS];
    /// radiation lengths of remainder of layer
    /// we count from the bottom of the converter to the top of the next lower converter
    double m_radLenRest[NLAYERS];
    /// number of layers of each type
    int m_numLayers[NTYPES];
    /// average radlen of converter
    double m_aveRadLenConv[NTYPES];
    /// average radlen of the rest
    double m_aveRadLenRest[NTYPES];

    /// Returns minimum trayHeight... I hope we can stop using this soon
    StatusCode getMinTrayHeight(double trayHeight);

    /// returns z position of X, Y or average plane for each layer
    StatusCode fillLayerZ();

    /// fill rad lengths
    StatusCode fillPropagatorInfo();

	/// get the relevant CAL info();
	StatusCode getCalInfo();

    /// return converter type for a layer
    convType getDigiLayerType(int digiLayer) const;


    /// pointer to the detector service
    IGlastDetSvc * m_pDetSvc;

    /// array to hold the tower part of the volumeIds of the silicon planes
    idents::VolumeIdentifier m_volId_tower[NTOWERS];
    /// array to hold the tray part of the volumeIds of the silicon planes
    idents::VolumeIdentifier m_volId_layer[NLAYERS][NVIEWS];

    /// Pointer to the old-style propagator needed by the track fit
    IKalmanParticle* m_KalParticle;
    /// pointer to the new propagation tool
    IPropagator* m_G4PropTool;
    /// pointer to the failure mode service
    ITkrFailureModeSvc* m_tkrFail;
    /// pointer to alignment service
    ITkrAlignmentSvc*   m_tkrAlign;
    /// pointer to bad strips service
    ITkrBadStripsSvc*   m_badStrips;
    /// pointer to splits service
    ITkrSplitsSvc*      m_tkrSplits;
    /// number of the test tower (used to find zLayer, etc.)
    int m_testTower;
    /// pointer to the ToT service
    ITkrToTSvc*         m_tkrToT;

};

#endif // TKRGEOMETRYSVC_H
