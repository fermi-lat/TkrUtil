
#ifndef TKRBADSTRIPSSVC_H
#define TKRBADSTRIPSSVC_H 

/** 
* @class BadStripsSvc
*
* @brief Maintains lists of bad strips, and provides access methods 
*
* First version 3-Jun-2001
*
* The bad strips are currently kept in ascii files, which are read in under
* the control of the jobOptions file. In the ascii files, strips are
* marked as hot or dead, but in memory, strips are only bad.
*
* This will be changing shortly when we interface to Joanne's calibration 
* database.
*
* The service creates an array of vectors.  The singly indexed array
* corresponds to a doubly indexed array by tower and layer.
* 
* The original design was a vector of vectors, but this was abandoned
* because some of the code failed on unix.
* 
* The use of the bad strips in the clustering algorithm depends
* on mixing good and bad strips and still being able to sort in ascending 
* strip order.
*
* To this end strips are tagged good or bad by adding a high-order bit. 
* Thus a tagged strip doesn't produce a legal strip number. The service has 
* methods for manipulating the tags, which are currently accessed by 
* TkrMakeClusters.
*
* The tagging arrangement can be extended, if necessary. For Example, 
* we may want to differentiate dead from hot strips in the reconstruction.
*
* @author Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Services/TkrBadStripsSvc.h,v 1.9 2002/09/18 23:41:00 lsrea Exp $
*/

#include "GaudiKernel/Service.h"

#include "TkrUtil/ITkrBadStripsSvc.h"

#include <string>
#include <vector>

class TkrBadStripsSvc : public Service, virtual public ITkrBadStripsSvc
{
public:
    
    /// Constructor of this form must be provided
    TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrBadStripsSvc() {}
    
    //. reads in bad strips and creates in-memory vectors
    StatusCode initialize();
    StatusCode finalize();
    
    /// converts from (tower, layer, view) to index into array
    int getIndex(int tower, int layer, 
        idents::GlastAxis::axis) const ;
    /// adds a strip to a badstrip vector
    void addStrip(stripCol* v, TaggedStrip taggedStrip);
    /// returns a pointer to a vector of bad strips for a given 
    /// (tower, layer and view)
    const stripCol* getBadStrips(int tower, int layer, 
        idents::GlastAxis::axis) const;
    /// returns a pointer to a vector of bad strips for a given array index
    const stripCol* getBadStrips(int index) const;
    /// returns true if the strip (tower, layer, view, strip) is bad
    bool isBadStrip(int tower, int layer, 
        idents::GlastAxis::axis, int strip) const;
    /// returns true if the given strip is found in the vector pointed 
    /// to by stripCol
    bool isBadStrip(const stripCol* v, int strip) const;
    
    /// queryInterface - required for a service
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);
    /// required for a service
    static const InterfaceID& interfaceID() { 
        return ITkrBadStripsSvc::interfaceID(); 
    }
    /// returns the service type
    const IID& type() const;   
    
private:
    
    //void makeCol(const int size); // left over from original attempt
    
    /// reads bad strips from file file
    void readFromFile(std::ifstream* file);

    /// File name for constants
    std::string m_badStripsFile;  
    
    /// implicit dimension of 256 array
    enum {NLAYERS = 18, NVIEWS = 2, NTOWERS = 16, NELEMENTS = NLAYERS*NVIEWS*NTOWERS};
    
    /// array to hold bad strips vectors  [ max needed: 576 = 16*18*2 ]   
    stripCol m_stripsCol[NELEMENTS];
};

    
#endif // TKRBADSTRIPSSVC_H
    
    
