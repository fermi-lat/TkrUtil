#ifndef __ITKRALIGNMENTSVC_H
#define __ITKRALIGNMENTSVC_H 1

#include "GaudiKernel/IInterface.h"

#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"

#include <string>
#include <vector>

//----------------------------------------------
//
//   TkrAlignmentSvc
//
//   Tracker Alignment Service. Supplies the alignment constants
// 
//----------------------------------------------
//             Leon Rochester, 23-Jan-2003
//----------------------------------------------

static const InterfaceID IID_ITkrAlignmentSvc("ITkrAlignmentSvc", 1 , 0); 

enum constType {SIM=0, REC=1};

/// A small class to define alignment constants

class AlignmentConsts 
{
public:
    
    AlignmentConsts (double deltaX=0., double deltaY=0., double deltaZ=0.,
        double rotX=0., double rotY=0., double rotZ=0.):
        m_deltaX(deltaX), m_deltaY(deltaY), m_deltaZ(deltaZ),
        m_rotX(rotX), m_rotY(rotY), m_rotZ(rotZ)
    {}
    
    ~AlignmentConsts() {} 

    double getDeltaX()  const {return m_deltaX;}
    double getDeltaY()  const {return m_deltaY;}
    double getDeltaZ()  const {return m_deltaZ;}
    double getRotX()    const {return m_rotX;}
    double getRotY()    const {return m_rotY;}
    double getRotZ()    const {return m_rotZ;}
    
private:

    double m_deltaX;
    double m_deltaY;
    double m_deltaZ;
    double m_rotX;
    double m_rotY;
    double m_rotZ;   
};

class ITkrAlignmentSvc : public virtual IInterface
{
public:

    static const InterfaceID& interfaceID() { return IID_ITkrAlignmentSvc; }
   
    virtual AlignmentConsts* getConsts(constType type, int tower, int layer, int view, 
        int ladder=0, int wafer=0) const = 0;

    virtual AlignmentConsts* getConsts(constType type, idents::VolumeIdentifier id) const = 0;

    virtual void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit) const = 0;

    virtual void moveCluster(int tower, int layer, int view, int ladder,
        HepPoint3D& point) const = 0;

    virtual bool alignSim() const = 0;
    virtual bool alignRec() const = 0;
};

#endif