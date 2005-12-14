
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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrGeometrySvc.h,v 1.25 2005/08/17 00:41:30 lsrea Exp $
 */

#include "GaudiKernel/Service.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"


class TkrGeometrySvc : public Service,
        virtual public ITkrGeometrySvc
{
public:
    enum { NVIEWS=2 /*, NLAYERS=18, NTOWERS=16, NPLANES=36*/};
    
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
    double siActiveWaferSide() const {return m_siWaferSide - 2.*m_siDeadDistance;}
    double ladderPitch()    const {return m_siWaferSide + m_ladderGap; }
    double waferPitch()     const {return m_siWaferSide + m_ladderInnerGap; }
    int    chipsPerLadder() const {return m_chipsPerLadder; }
    int    stripsPerChip () const {return m_ladderNStrips/m_chipsPerLadder; }

    int    numPlanes()      const {return m_numPlanes; }

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

    double gettkrZBot() const  {return m_tkrZBot; }
    
    double calZTop()        const {return m_calZTop;}
	double calZBot()        const {return m_calZBot;}
	double calXWidth()      const {return m_calXWidth;}
	double calYWidth()      const {return m_calYWidth;}

    /// reverse the numbering of the bilayers (goes both ways)
    int reverseLayerNumber(int layer) const {return numLayers()-layer-1;}
    /// same for planes
    int reversePlaneNumber(int plane) const {return numPlanes()-plane-1;}

    /// return the position of a strip, will accept int or double
    HepPoint3D getStripPosition(int tower, int layer, int view, double stripid) const;

    /// return the z position for a reconLayer and view
    double getReconLayerZ(int reconLayer, int view=2) const;
    /// return the z position for a digiLayer and view
    double getLayerZ     (int digiLayer,  int view=2) const;
    /// returns Z of *Layer* (average of x and y plane)
    /// TkrId of either plane will work
    double getLayerZ     (const idents::TkrId& tkrId) const {
        return getLayerZ(getLayer(tkrId));
    }
    /// z-position of converters
    double getConvZ      (int layer) const {
        if (getLayerType(layer)!=NOCONV) { 
            return m_convZ[layer]; 
        } else {
            // this is mainly so that a non-crazy value is returned
            //   in fact this will never be asked for if things are working correctly
            //   (except maybe for reverse tracking?)
            return std::max(getLayerZ(layer, 0), getLayerZ(layer,1));
        }
    }

    /// new stuff, based on plane and TkrId;
    int    getPlane (const idents::TkrId& tkrId) const {
        return 2*tkrId.getTray() + tkrId.getBotTop() - getBottomTrayFlag();
    }
	/// returns the nearest plane to the given z
    int getPlane(double z) const;

    /// returns number of planes between two objects specified by TkrId
    int getPlaneSeparation(const idents::TkrId& id1, const idents::TkrId& id2) const;

    double getPlaneZ(int plane) const {
        return m_planeZ[plane];
    }
    double getPlaneZ(const idents::TkrId& tkrId) const {
        return getPlaneZ(getPlane(tkrId));
    }
    int    getLayer (int plane) const {
        return m_planeToLayer[plane];
    }
    int    getLayer (const idents::TkrId& tkrId) const {
        return getLayer(getPlane(tkrId));
    }
    int    getView  (int plane) const {
        return m_planeToView [plane];
    }
    int    getView  (const idents::TkrId& tkrId) const {
        return getView(getPlane(tkrId));
    }
    bool   isTopPlaneInLayer(int plane) const {
        return m_isTopPlaneInLayer[plane];
    }

    /// return the rad length of the converter for each layer
    double getReconRadLenConv(int layer) const { 
        return getRadLenConv(reverseLayerNumber(layer));
    }
    double getRadLenConv(int layer) const { 
        return m_radLenConv[layer];
    }
    /// return the rad length of the rest of the layer, for each layer
    ///   counting down from the bottom of the converter
    double getReconRadLenRest(int layer) const { 
        return getRadLenRest(reverseLayerNumber(layer));
    }
    double getRadLenRest(int layer) const { 
        return m_radLenRest[layer];
    }

    /// return converter type for a layer
    convType getReconLayerType(int layer) const;
    convType getLayerType(int layer) const;

    /// return number of layers of each type
    int getNumType(convType type) const { return m_numLayers[(int)type];}
    /// get average radlen of converter for each type
    double getAveConv(convType type) const { return m_aveRadLenConv[(int)type];}
    /// get average radlen of rest for each type)
    ///    counting down from the bottom of the converter
    double getAveRest(convType type) const {return m_aveRadLenRest[(int) type];}

    ///does the tower exist?
    bool isTower(int tower)     const {return m_towerType[tower]>-1;}
    /// get the tower Type
    int getTowerType(int tower) const {return m_towerType[tower];}
    /// get limiting tower
    int getLimitingTower(int view, limitType type) const {
        if (view<0 || view> 1) return -1;
        return (view==0 ? m_xLim[type] : m_yLim[type]);
    }
    /// get the limits of the LAT as implemented
    double getLATLimit (int view, limitType type) const;
    /// are we in the "active" LAT?
    bool isInActiveLAT (Point pos) const;

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

        // definitions of plane, layer
    int trayToPlane(int tray, int botTop) const {
        return 2*tray + botTop - getBottomTrayFlag() ;
    }
    int trayToBiLayer(int tray, int botTop) const {
        return tray + botTop - getBottomTrayFlag() ;
    }
    /*
    void planeToTray(int plane, int& tray, int& face) const {
        tray = planeToTray(plane);
        face = planeToBotTop(plane);
    }
    */
    int planeToTray(int plane) const {
        return (plane+getBottomTrayFlag())/2;
    }
    int planeToBotTop(int plane) const {
        return (plane+getBottomTrayFlag())%2;
    }
    int getBottomTrayFlag() const { return (m_bottomTrayNumber>-1 ? 1 : 0); }
    int getTopTrayFlag()    const { return (m_topTrayNumber>-1    ? 1 : 0); }

    unsigned int getDefaultClusterStatus() const;

    double truncateCoord( double x, double pitch, 
        int numElements, int& elementNumber, bool reverse = false) const;

private:
    
    /// number of Towers in X
    int    m_numX; 
    /// number of Towers in Y
    int    m_numY;              
    /// two views, always!
    int    m_nviews;        
    int m_nWaferAcross;
    int m_numPlanes;
    int m_numTowers;
    int m_numTrays;

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
    ///
    int    m_chipsPerLadder;
    
    /// distance between two strips
    double m_siStripPitch; 
    /// nominally, (strip pitch)/sqrt(12)
    double m_siResolution;
    /// resolution factor, default = 1/sqrt(12)
    double m_siResolutionFactor;
    /// thickness of the silicon
    double m_siThickness;
    /// width of the dead region around the edge of a wafer
    double m_siDeadDistance;
    /// size of Wafer
    double m_siWaferSide;
    /// nominal bottom of the TKR envelope
    double m_tkrZBot;
    /// z coordinate of top of the top CAL crystal
	double m_calZTop;
	/// z coordinate of bottom of the bottom CAL crystal
	double m_calZBot;
	/// (maximum) width in x of active CsI across the entire instrument
	double m_calXWidth;
	/// (maximum) width in y of active CsI across the entire instrument
	double m_calYWidth;
    /// z positions of all the planes (digi convention)
    std::vector<double> m_planeZ;
    /// the two planes in a layer are separated by more than this:
    double m_layerSeparation;
    /// radiation lengths of converter, by recon layer *** really digi layer, I think!
    std::vector<double> m_radLenConv;
    /// radiation lengths of remainder of layer
    /// we count from the bottom of the converter to the top of the next lower converter
    std::vector<double> m_radLenRest;
    /// position of radiators
    std::vector<double> m_convZ;
    /// number of layers of each type
    int m_numLayers[NTYPES];
    /// average radlen of converter
    double m_aveRadLenConv[NTYPES];
    /// average radlen of the rest
    double m_aveRadLenRest[NTYPES];

    std::vector<int> m_planeToView;
    std::vector<int> m_planeToLayer;
    std::vector<bool> m_isTopPlaneInLayer;
    std::vector<int> m_layerToPlane[NVIEWS];

    int m_topTrayNumber;
    int m_bottomTrayNumber;

    /// Retrieves the basic constants from GlastDetSvc
    StatusCode getConsts();
    /// Initializes arrays of private data
    void initializeArrays();
    ///
    void makeTowerIds();
    ///
    void makeLayerIds();
    ///
    StatusCode getTowerLimits();
    ///
    void getTowerType();

    /// Returns minimum trayHeight... I hope we can stop using this soon
    StatusCode getMinTrayHeight(double& trayHeight);

    /// fill rad lengths
    StatusCode fillPropagatorInfo();

	/// get the relevant CAL info();
	StatusCode getCalInfo();

    /// find the test tower
    StatusCode getTestTower();
    /// find the legal volumes and store the info
    StatusCode getVolumeInfo();

    /// pointer to the detector service
    IGlastDetSvc * m_pDetSvc;

    /// array of tower types; type is number of exposed edges, -1 means no tower
    std::vector<int>  m_towerType;
    /// lowest and highest actual tower in x and y
    int  m_xLim[2];
    int  m_yLim[2];
    /// array to hold the tower part of the volumeIds of the silicon planes
    std::vector<idents::VolumeIdentifier> m_volId_tower;
    /// array to hold the tray part of the volumeIds of the silicon planes
    std::vector<idents::VolumeIdentifier> m_volId_layer[NVIEWS];

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
    /// root for testTower
    idents::VolumeIdentifier m_testTowerId;
    /// pointer to the ToT service
    ITkrToTSvc*         m_tkrToT;

};

#endif // TKRGEOMETRYSVC_H
