/** 
 @file TkrAlignmentSvc.h
 @brief Maintains list of alignment constants

 First version 23-Jan-2003
 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrAlignmentSvc.h,v 1.7 2004/06/12 05:03:47 lsrea Exp $
*/

#ifndef TKRALIGNMENTSVC_H
#define TKRALIGNMENTSVC_H 


#include "GaudiKernel/Service.h"
#include "GaudiKernel/ContainedObject.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

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
    AlignmentConsts getConsts() { return m_consts; }
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
        VolumeType type);
    
    /// called to signal end of nesting 
    virtual void popShape();
    
    /// Need a setMode in case someone wants other than default 
    virtual void setMode(std::string pmode) {m_mode = pmode;}
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
    
    /// Constructor of this form must be provided
    TkrAlignmentSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrAlignmentSvc() {}
    
    /// reads in alignment data and creates in-memory constants
    StatusCode initialize();
    /// not much here
    StatusCode finalize();
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts
        (constType type, int tower, int layer, int view, int ladder, int wafer) const;
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts
        (constType type, idents::VolumeIdentifier id) const;
    
    /// moves the entry and exit of an MCPositionHit according to alignment consts
    void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit) const;
    
    /// moves a TkrCluster according to deltaX and deltaY alignment consts
    void moveCluster(int tower, int layer, int view, int ladder, HepPoint3D& point) const;

    void moveReconHit(int tower, int layer, int view, int ladder,
        HepPoint3D& point, HepVector3D dir) const;

    /// Get the volId and the local coordinates for the point to be aligned
    virtual idents::VolumeIdentifier getGeometryInfo(int layer, int view, HepPoint3D globalPoint, 
        HepPoint3D& alignmentPoint) const;
   
    /// true means make alignment corrections at digi time
    bool alignSim() const {return ((m_fileFlag&(1<<SIM_SHIFT))>0); }
    /// true means make alignment corrections at recon time
    bool alignRec() const {return ((m_fileFlag&(1<<REC_SHIFT))>0); }
    
    /// queryInterface - required for a service
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);
    /// required for a service
    static const InterfaceID& interfaceID() { 
        return ITkrAlignmentSvc::interfaceID(); 
    }
    /// returns the service type
    const IID& type() const;    
    
private:
    
    
    /// returns the constants for a given element
    const AlignmentConsts* getConsts(constType type, int index) const;
    
    /// to calculate delta of a point
    HepVector3D getDelta(int view,  const HepPoint3D& point, const HepVector3D& dir,
        const AlignmentConsts* alConsts) const;
    
    /// returns index
    int getIndex(int tower, int layer, int view, int ladder, int wafer) const;
    
    /// gets relevant coordinates from the geometry
    StatusCode getGeometry();
    /// gets the alignment data, or prepares the test case
    StatusCode getData(std::string fileName);
    /// runs the tests
    StatusCode doTest();
    /// fill the tray constants from the tower constants
    StatusCode fillTrayConsts();
    /// do the transformation from tower to tray
    void calculateTrayConsts(AlignmentConsts& thisTray);
    /// fill the face constants
    StatusCode fillFaceConsts();
    /// do the transformation from tray to face
    void calculateFaceConsts(AlignmentConsts& thisPlane);
    /// fill the ladder constants
    StatusCode fillLadderConsts();
    /// do the transformation from face to ladder
   void calculateLadderConsts(AlignmentConsts& thisLadder);

    /// fill the wafer constants
    StatusCode fillWaferConsts();
    /// do the transformation from ladder to wafer
    void calculateWaferConsts(AlignmentConsts& thisWafer);
   
    /// reads the alignment items from a file
    StatusCode readFromFile();
    /// generate the per-wafer alignment constants
    StatusCode processFile();
    /// does some logic on the next item in the read-in list
    bool getNextItem(aType type, AlignmentItem& item);
    
    /// File name for sim constants
    std::string m_simFile;
    /// File name for rec constants
    std::string m_recFile;
    /// test mode flag
    int m_testMode;
    /// file flag: bit 0 means sim file exists, bit 1 means rec file exists
    int m_fileFlag;
    /// maximum allowed delta (sqrt(deltaX**2 + sqrt(deltaY**2)
    double m_maxDelta;
    
    /// dimension of arrays
    // the flight instrument only has (4 ladders)x(4 wafers) but to allow
    //   for possible conversion to BFEM/BTEM I've started with 5 each.
    enum {NLAYERS = 18, NVIEWS = 2, NTOWERS = 16, NLADDERS= 4, NWAFERS = 4,
        NELEMENTS = NLAYERS*NVIEWS*NTOWERS*NLADDERS*NWAFERS};
    
    /// array to hold simulation constants  [ max needed: 9216 = 16*18*2*4*4 ]   
    AlignmentConsts m_simConsts[NELEMENTS];
    /// array to hold reconstruction constants  [ max needed: 9216 = 16*18*2*4*4 ]   
    AlignmentConsts m_recConsts[NELEMENTS];
    
    // the following consts are used to construct the wafer constants 
    //   from the tower, tray and ladder constants.
    
    /// z of tray center in tower
    double m_trayZ[NLAYERS+1];
    /// z of plane in tray, vs. tray number and botTop
    double m_faceZ[NLAYERS+1] [NVIEWS];
    
    /// holds alignment consts for the towers
    AlignmentConsts m_towerConsts;
    /// holds alignment consts for the trays
    AlignmentConsts m_trayConsts;
    /// holds alignment consts for the planes
    AlignmentConsts m_faceConsts;
    /// holds alignment consts for the ladders
    AlignmentConsts m_ladderConsts;
    /// hold alignment consts for the wafers
    AlignmentConsts m_waferConsts;
    
    /// current element being constructed
    std::string m_mode;
    std::ifstream* m_dataFile;
    int m_tower;
    int m_tray;
    int m_face;
    int m_ladder;
    int m_wafer;
    
    ITkrGeometrySvc* m_pGeoSvc;
    IGlastDetSvc* m_pDetSvc;

    typedef std::vector<AlignmentItem*> alignCol;
    /// list of items read in
    alignCol m_itemCol;
    /// item being processed
    alignCol::iterator m_pItem;
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

inline std::ostream& AlignmentConsts::fillStream( std::ostream& s ) const {
    s << "class AlignmentConsts: "
        << "delta(X,Y,Z) (" 
        << m_deltaX << ", "
        << m_deltaY << ", "
        << m_deltaZ << ") "
        << " rot(X,Y,Z) ("
        << m_rotX << ", "
        << m_rotY << ", "
        << m_rotZ << ")"
        << std::endl;
    return s;
}

#endif // TKRALIGNMENTSVC_H


