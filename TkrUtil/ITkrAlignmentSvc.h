/** @file ITkrAlignmentSvc.h
@brief AlignmentConsts class & Abstract interface to TkrAlignmentSvc (q.v.) 
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrUtil/TkrUtil/ITkrAlignmentSvc.h,v 1.21 2008/05/20 01:14:20 lsrea Exp $
*/


#ifndef __ITKRALIGNMENTSVC_H
#define __ITKRALIGNMENTSVC_H 1

#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/ContainedObject.h"

#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "CalibData/Tkr/TkrTowerAlignCalib.h"
#include "CalibData/Tkr/TkrInternalAlignCalib.h"

#include <string>
#include <vector>
#include <iostream>

// TU: Hacks for CLHEP 1.9.2.2 and beyond
#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#ifndef HepVector3D
typedef HepGeom::Vector3D<double> HepVector3D;
#endif

static const InterfaceID IID_ITkrAlignmentSvc("ITkrAlignmentSvc", 8, 0); 

namespace {
    enum alignTask {NULLTASK= 0, APPLYCONSTS=1, FINDTOWERCONSTS=2, FINDWAFERCONSTS=3};

	// bits in the flags word: 
	// 0  0  X  X | X  X  X  X
	//
	//       r  r   r  d  d  d
	//       o  o   o  z  y  x
	//       t  t   t
	//       z  y   x

}


/// A small class to define alignment constants

/** 
  @class AlignmentConsts 
  @brief holds the six alignment constants  
 */

class AlignmentConsts 
{
public:
    
    AlignmentConsts (double deltaX=0., double deltaY=0., double deltaZ=0.,
        double rotX=0., double rotY=0., double rotZ=0.):
        m_deltaX(deltaX), m_deltaY(deltaY), m_deltaZ(deltaZ),
        m_rotX(rotX), m_rotY(rotY), m_rotZ(rotZ)
    {}
    
    /*
    AlignmentConsts (CLHEP::Hep3Vector disp = CLHEP::Hep3Vector(0., 0., 0.),
        CLHEP::Hep3Vector rot = CLHEP::Hep3Vector(0., 0., 0.)):
        m_deltaX(disp.x()), m_deltaY(disp.y()), m_deltaZ(disp.z()),
        m_rotX(rot.x()), m_rotY(rot.y()), m_rotZ(rot.z())
    {}
    */

   ~AlignmentConsts() {} 

    /// get the displacement in X
    double getDeltaX()  const {return m_deltaX;}
    /// get the displacement in Y
    double getDeltaY()  const {return m_deltaY;}
    /// get the displacement in Z
    double getDeltaZ()  const {return m_deltaZ;}
    /// get the rotation around the X axis
    double getRotX()    const {return m_rotX;}
    /// get the rotation around the Y axis
    double getRotY()    const {return m_rotY;}
    /// get the rotation around the Z axis
    double getRotZ()    const {return m_rotZ;}

    /// check for a zero AlignmentConsts
    bool isNull() { return (m_deltaX==0. && m_deltaY==0. && m_deltaZ==0.
        && m_rotX==0. &&m_rotY==0. && m_rotZ==0.);} 

    //! Fill the ASCII output stream
    friend std::ostream& operator<<( std::ostream& s , AlignmentConsts consts);  
        
    //! Serialize the object for reading
    StreamBuffer& serialize( StreamBuffer& s );
    //! Serialize the object for writing
    StreamBuffer& serialize( StreamBuffer& s ) const;
    
private:

    /// displacement in X
    double m_deltaX;
    /// displacement in Y
    double m_deltaY;
    /// displacement in Z
    double m_deltaZ;
    /// rotation around X axis
    double m_rotX;
    /// rotation around Y axis
    double m_rotY;
    /// rotation around Z axis
    double m_rotZ;   
};

/** 
  @class ITkrAlignmentSvc 
  @brief Abstract interface to TkrAlignmentSvc (q.v.)  
  @author Leon Rochester 23-Jan-2003
 */


class ITkrAlignmentSvc : public virtual IInterface
{
public:


    static const InterfaceID& interfaceID() { return IID_ITkrAlignmentSvc; }
   
    enum calibType {SIM=0, REC, NCALIBTYPES, UNKNOWN_TYPE};

	enum alignMode {XANDY=0x3, ROTZ=0x20, ANGLE=0x1c, USEALL=0x3f};

    /// retrieve the alignment consts for element tower, layer, view, ladder, wafer
    virtual const AlignmentConsts* getConsts(calibType type, int tower, 
        int layer, int view, int ladder=0, int wafer=0) const = 0;
    /// retrieve the alignment consts for a given volume Id
    virtual const AlignmentConsts* getConsts(calibType type, 
        idents::VolumeIdentifier id) const = 0;
    /// move the McHit by the alignment consts
    virtual void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit, 
        HepVector3D dir=HepVector3D(0.0, 0.0, 0.0)) const = 0;

    /// move the recon hit by the alignment consts
  //  virtual void moveReconPoint(HepPoint3D& point, const HepVector3D& dir, 
  //      int layer, int view, alignTask task = APPLYCONSTS, 
  //      const AlignmentConsts* consts = 0,
  //      const unsigned flags = USEALL) const = 0;

    /// move the recon hit by the alignment consts
    virtual HepVector3D deltaReconPoint(
        const HepPoint3D& point, const HepVector3D& dir, 
        int layer, int view, 
        unsigned flags = USEALL, 
        alignTask task = APPLYCONSTS, 
        const AlignmentConsts* consts = 0 
    ) const = 0;

    /// Get the volId and the local coordinates for the point to be aligned
    //virtual idents::VolumeIdentifier getGeometryInfo(int layer, int view, 
    //    const HepPoint3D& globalPoint, HepPoint3D& alignmentPoint) const = 0;
    virtual HepPoint3D getTowerCoordinates(const HepPoint3D& globalPoint,
        int& nXTower, int& nYTower) const = 0;
    /// true = perform alignment during simulation
    virtual bool alignSim() const = 0;
    /// true = perform alignment during reconstruction
    virtual bool alignRec() const = 0;

    /// update to latest pointer when calibration changes
    virtual void update(CalibData::TkrTowerAlignCalib* pTowerAlign, 
        CalibData::TkrInternalAlignCalib* pInternalAlign) = 0;
    /// set SIM/REC type
    virtual void SetCalibType(calibType type) const = 0;

};

typedef bool useFlags [6];

#endif
