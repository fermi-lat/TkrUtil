/** 
 @file TkrAlignmentSvc.h
 @brief Maintains list of alignment constants

 First version 23-Jan-2003
 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrAlignmentSvc.h,v 1.24 2013/04/26 18:08:25 lsrea Exp $
*/

#ifndef TKRALIGNMENTSVC_H
#define TKRALIGNMENTSVC_H 


#include "GaudiKernel/Service.h"
#include "GaudiKernel/ContainedObject.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "TkrUtil/IndexedVector.h"

#include "GaudiKernel/ContainedObject.h"


#include <string>
#include <vector>
#include <iostream>


namespace {
    const int ntypes = 6;
    const std::string itemType[ntypes] = {"TOWER", "TRAY", "FACE", "LADDER", "WAFER", "NONE"};
    enum aType {TOWER, TRAY, FACE, LADDER, WAFER, NONE};
}

/// A small class to define an alignment const item

/** 
  @class AlignmentItem
  @brief holds an element type, an element number and an AlignmentConsts, read in from a file
 */
class AlignmentItem
{
public:

    AlignmentItem (int type, int number=0, AlignmentConsts consts=AlignmentConsts()) :
        m_number(number), m_consts(consts)
    {     
        switch (type) {
            case 0: m_type = TOWER; break;
            case 1: m_type = TRAY;  break;
            case 2: m_type = FACE;  break;
            case 3: m_type = LADDER; break;
            case 4: m_type = WAFER;  break;
            default: m_type = NONE;
        }
    }

    AlignmentItem (aType type, int number=0, AlignmentConsts consts=AlignmentConsts()) :
        m_type(type), m_number(number), m_consts(consts)
    {}
    AlignmentItem () :
        m_type(TOWER), m_number(0), m_consts(AlignmentConsts())
    {}
    
    ~AlignmentItem() {} 

    /// return the type
    aType getType()             { return m_type; }
    /// return the index
    int   getNumber()           { return m_number; }
    /// return the constants
    const AlignmentConsts getConsts() { return m_consts; }
private:
    /// type, TOWER, TRAY, FACE etc.
    aType m_type;
    /// index of the type (i.e., which tray)
    int   m_number;
    /// aligment consts with respect to containing element
    AlignmentConsts m_consts;
    
};

#include "GlastSvc/GlastDetSvc/IGeometry.h"

/** 
 @class TkrAlignmentGeomVisitor
 @brief Retrieves the necessary geometry info
*/

class TkrAlignmentGeomVisitor: public IGeometry
{
public:
    
    TkrAlignmentGeomVisitor() {}
    
    ~TkrAlignmentGeomVisitor() {}
    
    /// Standard interface to the detModel
    virtual IGeometry::VisitorRet pushShape(ShapeType s, const UintVector& id, 
        std::string name, std::string material, const DoubleVector& params, 
        VolumeType type, SenseType sense);
    
    /// called to signal end of nesting 
    virtual void popShape();
    
    /// Need a setMode in case someone wants other than default 
    virtual void setMode(std::string pmode) {
        m_mode = pmode;
    }
    /// Implements getMode for IGeometry interface
    virtual std::string getMode() {return m_mode;}
    
    double getTrayBotHeight() { return m_trayBotHeight;}
    double getTrayTopHeight() { return m_trayTopHeight;}
private:
    
    double m_trayBotHeight;
    double m_trayTopHeight;
    
    /// mode for traversing geometry
    std::string m_mode;
};


/** 
 @class TkrAlignmentSvc
 @brief Maintains list of alignment constants
*/

class TkrAlignmentSvc : public Service, virtual public ITkrAlignmentSvc
{
public:
    
    enum {SIM_SHIFT = 0, REC_SHIFT = 1};
    
    /// sets type to SIM or REC
    void SetCalibType(calibType type) const {m_calibType = type;}

    /// Constructor of this form must be provided
    TkrAlignmentSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrAlignmentSvc() {}
    
    /// reads in alignment data and creates in-memory constants
    StatusCode initialize();
    /// not much here
    StatusCode finalize();
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts
        (calibType type, int tower, int layer, int view, int ladder, int wafer) const;
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts
        (calibType type, idents::VolumeIdentifier id) const;
    
    /// moves the entry and exit of an MCPositionHit according to alignment consts
    void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit, HepVector3D dir) const;
    
    HepVector3D deltaReconPoint(const HepPoint3D& point, const HepVector3D& dir, 
        int layer, int view, 
        alignTask task, 
        const AlignmentConsts* consts
    ) const;
	
    /// Get the volId and the local coordinates for the point to be aligned
    idents::VolumeIdentifier getGeometryInfo(int layer, int view, 
        const HepPoint3D& globalPoint, HepPoint3D& alignmentPoint) const;

    HepPoint3D getTowerCoordinates(const HepPoint3D& globalPoint,
        int& nXTower, int& nYTower) const;
   
    /// true means make alignment corrections at digi time
    bool alignSim() const {return ((m_fileFlag&(1<<SIM_SHIFT))>0); }
    /// true means make alignment corrections at recon time
    bool alignRec() const {return ((m_fileFlag&(1<<REC_SHIFT))>0); }
    
    /// queryInterface - required for a service
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);
    /// required for a service
    static const InterfaceID& interfaceID() { 
        return ITkrAlignmentSvc::interfaceID(); 
    }
    /// returns the service type
    const InterfaceID& type() const; 

    /// update to latest pointer when calibration changes
    void update(CalibData::TkrTowerAlignCalib* pTowerAlign, 
        CalibData::TkrInternalAlignCalib* pInternalAlign);

    void moveReconPoint(HepPoint3D& point, const HepVector3D& dir, 
        int layer, int view, alignTask task, const AlignmentConsts* consts) const;
    
private:   
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts(calibType type, int index) const;
    
    /// to calculate delta of a point
    HepVector3D getDelta(int view,  const HepPoint3D& point, const HepVector3D& dir,
        const AlignmentConsts* alConsts) const;

    void applyDelta(
        double pointX, double pointY, double alphaX, double alphaY,
        const AlignmentConsts* alConsts, 
        double& deltaPointX, double& deltaPointY
        ) const;

    /// set the const mode
    void setMode(std::string pmode) {
        m_mode = pmode;
        if     (m_mode=="sim") {m_calibType = SIM;}
        else if(m_mode=="rec") {m_calibType = REC;}
        else                   {m_calibType = UNKNOWN_TYPE;}
    }
    
    /// gets relevant coordinates from the geometry
    StatusCode getGeometry();
    /// gets the alignment data, or prepares the test case
    StatusCode getData(std::string fileName);
    /// runs the tests
    StatusCode doTest();
    /// fill the tray constants from the tower constants
    StatusCode fillTrayConsts();
    /// do the transformation from tower to tray
    void calculateTrayConsts(AlignmentConsts& thisTray) const;
    /// fill the face constants
    StatusCode fillFaceConsts();
    /// do the transformation from tray to face
    void calculateFaceConsts(AlignmentConsts& thisPlane) const;
    /// fill the ladder constants
    StatusCode fillLadderConsts();
    /// do the transformation from face to ladder
   void calculateLadderConsts(AlignmentConsts& thisLadder) const;

    /// fill the wafer constants
    StatusCode fillWaferConsts();
    /// do the transformation from ladder to wafer
    void calculateWaferConsts(AlignmentConsts& thisWafer) const;
   
    /// reads the alignment items from a file
    StatusCode readFromFile();
    /// generate the per-wafer alignment constants
    StatusCode processConstants();
    /// does some logic on the next item in the read-in list
    bool getNextItem(aType type, AlignmentItem& item);

    AlignmentConsts makeConsts(CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot) {
        return AlignmentConsts(0.001*disp.x(), 0.001*disp.y(), 0.001*disp.z(), 
            0.001*rot.x(), 0.001*rot.y(), 0.001*rot.z());
    }
    
    /// File name for sim constants, old text file system
    std::string m_simFile;
    /// File name for rec constants, old text file system
    std::string m_recFile;
    /// test mode flag
    int m_testMode;
    /// file flag: bit 0 means sim file exists, bit 1 means rec file exists
    int m_fileFlag;
    /// maximum allowed delta (sqrt(deltaX**2 + sqrt(deltaY**2)
    double m_maxDelta;
    /// scale for sim file (for testing)
    double m_simScale;
    /// scale for rec file (for testing)
    double m_recScale;
    /// same for each component;
    std::vector<double> m_simScaleVec;
    std::vector<double> m_recScaleVec;
        
    /* Although the information about the alignment is developed in terms of trays and faces,
       it's stored in the arrays by layer and view, since this is the preferred set of
       variables for reconstruction.

       Yes, I know it's confusing... I just confused myself about all this!
    */
    
    typedef IndexedVector<AlignmentConsts> ConstsVec;

    // array size will be nTowers*nLayers*nViews*nLadders*nWafers (currently nLadders = nWafers)
    /// array to hold simulation constants   
    mutable ConstsVec m_simConsts; 
    /// array to hold reconstruction constants  
    mutable ConstsVec m_recConsts; 
    
    // the following consts are used to construct the wafer constants 
    //   from the tower, tray and ladder constants.
    
    /// z of tray center in tower
    std::vector<double> m_trayZ;
    /// z of plane in tray, vs. tray number and botTop
    mutable IndexedVector<double> m_faceZ;
    
    /// holds alignment consts for the towers during construction
    mutable AlignmentConsts m_towerConsts;
    /// holds alignment consts for the trays during construction
    mutable AlignmentConsts m_trayConsts;
    /// holds alignment consts for the planes during construction
    mutable AlignmentConsts m_faceConsts;
    /// holds alignment consts for the ladders during construction
    mutable AlignmentConsts m_ladderConsts;
    /// hold alignment consts for the wafers during construction
    mutable AlignmentConsts m_waferConsts;
    /// current element being constructed
    mutable int m_tower;
    mutable int m_tray;
    mutable int m_face;
    mutable int m_ladder;
    mutable int m_wafer;

    std::string m_mode;
    mutable calibType m_calibType;
    bool m_printConsts;
    std::ifstream* m_dataFile;

    ITkrGeometrySvc* m_tkrGeom;
    IGlastDetSvc* m_pDetSvc;

    typedef std::vector<AlignmentItem*> alignCol;
    /// list of items read in
    alignCol m_itemCol;
    /// item being processed
    alignCol::iterator m_pItem;

    /// some useful constants
    int m_nTowers;
    int m_nTrays;
    int m_nLayers;
    int m_nViews;
    int m_nFaces;
    int m_nLadders;
    int m_nWafers;

};



//! Serialize the object for writing

inline StreamBuffer& AlignmentConsts::serialize( StreamBuffer& s ) const {
    s << m_deltaX
        << m_deltaY
        << m_deltaZ
        << m_rotX
        << m_rotY
        << m_rotZ;   
    return s;
}

//! Serialize the object for reading
inline StreamBuffer& AlignmentConsts::serialize( StreamBuffer& s )       {
    s >> m_deltaX
        >> m_deltaY
        >> m_deltaZ
        >> m_rotX
        >> m_rotY
        >> m_rotZ;
    
    return s;
}

//! Fill the ASCII output stream
std::ostream& operator<<(std::ostream &s, AlignmentConsts consts) {
    s << "AlignmentConsts: "
            << "delta(" 
            << std::setprecision(6)
        << consts.getDeltaX() << ", "
        << consts.getDeltaY() << ", "
        << consts.getDeltaZ() << ") "
        << " rot("
        << consts.getRotX() << ", "
        << consts.getRotY() << ", "
        << consts.getRotZ() << ")" ; 
    return s;
}

#endif // TKRALIGNMENTSVC_H


