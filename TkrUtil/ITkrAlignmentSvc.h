/** @file ITkrAlignmentSvc.h
@brief AlignmentConsts class & Abstract interface to TkrAlignmentSvc (q.v.) 
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrAlignmentSvc.h,v 1.16 2005/12/20 02:35:57 lsrea Exp $
*/


#ifndef __ITKRALIGNMENTSVC_H
#define __ITKRALIGNMENTSVC_H 1

#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/ContainedObject.h"

#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

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

static const InterfaceID IID_ITkrAlignmentSvc("ITkrAlignmentSvc", 6, 0); 

namespace {
    enum constType {SIM=0, REC=1, UNKNOWN_TYPE=2};
    enum alignTask {NULLTASK= 0, APPLYCONSTS=1, FINDTOWERCONSTS=2, FINDWAFERCONSTS=3};
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
   
    /// retrieve the alignment consts for element tower, layer, view, ladder, wafer
    virtual const AlignmentConsts* getConsts(constType type, int tower, 
        int layer, int view, int ladder=0, int wafer=0) const = 0;
    /// retrieve the alignment consts for a given volume Id
    virtual const AlignmentConsts* getConsts(constType type, 
        idents::VolumeIdentifier id) const = 0;
    /// move the McHit by the alignment consts
    virtual void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit) const = 0;
    /// move the recon hit by the alignment consts
    virtual void moveReconPoint(HepPoint3D& point, const HepVector3D& dir, 
        int layer, int view, alignTask task = APPLYCONSTS, 
        const AlignmentConsts* consts = 0) const = 0;
    /// move the recon hit by the alignment consts
    virtual HepVector3D deltaReconPoint(
        const HepPoint3D& point, const HepVector3D& dir, 
        int layer, int view, alignTask task= APPLYCONSTS, 
        const AlignmentConsts* consts = 0) const = 0;
    /// Get the volId and the local coordinates for the point to be aligned
    //virtual idents::VolumeIdentifier getGeometryInfo(int layer, int view, 
    //    const HepPoint3D& globalPoint, HepPoint3D& alignmentPoint) const = 0;
    virtual HepPoint3D getTowerCoordinates(const HepPoint3D& globalPoint,
        int& nXTower, int& nYTower) const = 0;
    /// true = perform alignment during simulation
    virtual bool alignSim() const = 0;
    /// true = perform alignment during reconstruction
    virtual bool alignRec() const = 0;
};

#endif
