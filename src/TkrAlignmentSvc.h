/** 
 @file TkrAlignmentSvc.h
 @brief Maintains list of alignment constants

 First version 23-Jan-2003
 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrAlignmentSvc.h,v 1.4 2003/04/11 20:51:47 lsrea Exp $
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
    ///
    StatusCode fillLadderConsts();
    ///
    StatusCode fillWaferConsts();
    
    /// reads bad strips from file file
    StatusCode readFromFile();
    
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
    enum {NLAYERS = 18, NVIEWS = 2, NTOWERS = 16, NLADDERS= 5, NWAFERS = 5,
        NELEMENTS = NLAYERS*NVIEWS*NTOWERS*NLADDERS*NWAFERS};
    
    /// array to hold simulation constants  [ max needed: 14400 = 16*18*2*5*5 ]   
    AlignmentConsts m_simConsts[NELEMENTS];
    /// array to hold reconstruction constants  [ max needed: 14400 = 16*18*2*5*5 ]   
    AlignmentConsts m_recConsts[NELEMENTS];
    
    // the following consts are used to construct the wafer constants 
    //   from the tower, tray and ladder constants.
    
    /// z of tray center in tower
    double m_trayZ[NLAYERS+1];
    /// z of plane in tray, vs. tray number and botTop
    double m_planeZ[NLAYERS+1] [NVIEWS];
    
    /// holds alignment consts for the towers
    AlignmentConsts m_towerConsts;
    /// holds alignment consts for the trays
    AlignmentConsts m_trayConsts;
    /// holds alignment consts for the ladders
    AlignmentConsts m_ladderConsts;
    /// hold alignment consts for the wafers
    AlignmentConsts m_waferConsts;
    
    /// current elememt being constructed
    std::string m_mode;
    std::ifstream* m_dataFile;
    int m_tower;
    int m_tray;
    int m_botTop;
    int m_ladder;
    int m_wafer;
    
    ITkrGeometrySvc* m_pGeoSvc;
    IGlastDetSvc* m_pDetSvc;
};


/*
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
*/

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


