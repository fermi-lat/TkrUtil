
#ifndef TKRALIGNMENTSVC_H
#define TKRALIGNMENTSVC_H 

/** 
* @class AlignmentSvc
*
* @brief Maintains list of alignment constants
*
* First version 23-Jan-2003
*
* @author Leon Rochester
*
* $Header$
*/

#include "GaudiKernel/Service.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

#include <string>
#include <vector>

class TkrAlignmentSvc : public Service, virtual public ITkrAlignmentSvc
{
public:

    enum {SIM_SHIFT = 0, REC_SHIFT = 1};
    
    /// Constructor of this form must be provided
    TkrAlignmentSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrAlignmentSvc() {}
    
    //. reads in bad strips and creates in-memory constants
    StatusCode initialize();
    StatusCode finalize();
    
    /// returns the constants for a given element
    AlignmentConsts* getConsts
        (constType type, int tower, int layer, int view, int ladder, int wafer) const;

    /// returns the constants for a given element
    AlignmentConsts* getConsts
        (constType type, idents::VolumeIdentifier id) const;

    /// moves the entry and exit of an MCPositionHit according to alignment consts
    void moveMCHit(idents::VolumeIdentifier id, 
        HepPoint3D& entry, HepPoint3D &exit) const;

    /// moves a TkrCluster according to deltaX and deltaY alignment consts
    void moveCluster(int tower, int layer, int view, int ladder, HepPoint3D& point) const;

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
    AlignmentConsts* getConsts(constType type, int index) const;

    /// to calculate delta of a point
    HepVector3D getDelta(int view,  const HepPoint3D& point, const HepVector3D& dir,
        const AlignmentConsts* alConsts) const;
    
    /// returns index
    int getIndex(int tower, int layer, int view, int ladder, int wafer) const;
    
    /// reads bad strips from file file
    void readFromFile(std::ifstream* simFile, AlignmentConsts* xxxConsts);
    
    /// File name for sim constants
    std::string m_simFile;
    /// File name for rec constants
    std::string m_recFile;
    /// test mode flag
    int m_testMode;
    /// file flag: bit 0 means sim file exists, bit 1 means rec file exists
    int m_fileFlag;
   
    /// dimension of arrays
    enum {NLAYERS = 18, NVIEWS = 2, NTOWERS = 16, NLADDERS= 5, NWAFERS = 5,
        NELEMENTS = NLAYERS*NVIEWS*NTOWERS*NLADDERS*NWAFERS};
    
    /// array to hold simulation constants  [ max needed: 14400 = 16*18*2*5*5 ]   
    AlignmentConsts m_simConsts[NELEMENTS];
    /// array to hold reconstruction constants  [ max needed: 14400 = 16*18*2*5*5 ]   
    AlignmentConsts m_recConsts[NELEMENTS];
};


#endif // TKRALIGNMENTSVC_H


