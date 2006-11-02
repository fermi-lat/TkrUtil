/** 
 @file TkrBadStripsSvc.h

 @brief Maintains lists of bad strips, and provides access methods 

 First version 3-Jun-2001
  @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrBadStripsSvc.h,v 1.12 2006/03/21 01:15:48 usher Exp $
*/

#ifndef TKRBADSTRIPSSVC_H
#define TKRBADSTRIPSSVC_H 

/** 
 @class TkrBadStripsSvc

 @brief Maintains lists of bad strips, and provides access methods 

 First version 3-Jun-2001

 The bad strips are kept in ascii files or in the TkrBadStrips database.
 If the files are made known to the service, they are used instead of 
 the database, providing a flexible development tool. The files are not
 for production use, and maintenance is not guarranteed.

 The files are read in under
 the control of the jobOptions file. In the ascii files, strips are
 marked as hot or dead, but in memory, strips are only bad.

 The service creates an array of vectors.  The singly indexed array
 corresponds to a doubly indexed array by tower and layer.
 
 The original design was a vector of vectors, but this was abandoned
 because some of the code failed on unix.
 
 The use of the bad strips in the clustering algorithm depends
 on mixing good and bad strips and still being able to sort in ascending 
 strip order.

 To this end strips are tagged good or bad by adding a high-order bit. 
 Thus a tagged strip doesn't produce a legal strip number. The service has 
 methods for manipulating the tags, which are currently accessed by 
 TkrMakeClusters.

 The tagging arrangement can be extended, if necessary. For Example, 
 we may want to differentiate dead from hot strips in the reconstruction.

*/

#include "GaudiKernel/Service.h"

#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrBadStripsSvcCalib.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/IndexedVector.h"

#include <string>

namespace {
/**
@class BadVisitor

  Minimal class derived from CalibData::BadStripsVisitor to
  check out BadStrips visitor interface.
*/
    class BadVisitor : public CalibData::BadStripsVisitor {
    public:
        BadVisitor(MsgStream* log=0) : m_log(log), m_nStrips(0) {}
        
        void setLog(MsgStream* log) {m_log = log;}
        
        virtual CalibData::eVisitorRet badTower(unsigned int row, unsigned int col,
            int badness);
        
        virtual CalibData::eVisitorRet badPlane(unsigned int row, unsigned int col, 
            unsigned int tray, bool top,
            int badness, bool allBad,
            const CalibData::StripCol& strips);
        
        void setService(ITkrBadStripsSvcCalib* pBadStrips) {m_pBadStrips = pBadStrips;}
        void setService(ITkrGeometrySvc*       tkrGeom)    {m_tkrGeom    = tkrGeom;}

        bool isEmpty() { return m_nStrips==0; }
        // manipulate flags to trigger generation of bad clusters
        bool generateBadClusters();
 
    private:
        MsgStream* m_log;
        ITkrBadStripsSvcCalib* m_pBadStrips;
        ITkrGeometrySvc*       m_tkrGeom;
        int m_nStrips;
    };
}

class TkrBadStripsSvc : public Service, virtual public ITkrBadStripsSvcCalib,
virtual public ITkrBadStripsSvc
{
public:
    
    /// Constructor of this form must be provided
    TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrBadStripsSvc() {}
    
    //. reads in bad strips and creates in-memory vectors
    StatusCode initialize();
    StatusCode finalize();

    /// sets type to SIM or REC
    void SetCalibType(calibType type) const {m_calibType = type;}
    
    /// adds a strip to a badstrip vector
    void addStrip(stripCol* v, TaggedStrip taggedStrip);
    /// returns a pointer to a vector of bad strips for a given 
    /// (tower, layer and view)
    const stripCol* getBadStrips(int tower, int layer, 
        idents::GlastAxis::axis) const;
    /// returns true if the strip (tower, layer, view, strip) is bad
    bool isBadStrip(int tower, int layer, 
        idents::GlastAxis::axis, int strip) const;
    /// returns true if the given strip is found in the vector pointed 
    /// to by stripCol
    bool isBadStrip(const stripCol* v, int strip) const;
		
    std::ostream& fillStream( std::ostream& s ) const;        

    bool empty() const { return m_empty[m_calibType]; }

    // for bad clusters
    void  setBadClusterCol(Event::TkrClusterCol* pClus)    {
        if (m_pBadClus[m_calibType])  {delete m_pBadClus[m_calibType];}
        m_pBadClus[m_calibType] = pClus; 
    }
    void  setBadIdClusterMap(Event::TkrIdClusterMap* pMap) {
        if(m_pBadMap[m_calibType]) { delete m_pBadMap[m_calibType];}
        m_pBadMap[m_calibType] = pMap; 
    }

    Event::TkrClusterCol*  getBadClusterCol() const    { return m_pBadClus[m_calibType]; }
    Event::TkrIdClusterMap* getBadIdClusterMap() const { return m_pBadMap[m_calibType]; }
    //Event::TkrDigiCol* getBadDigiCol() const           { return m_pBadDigi; }
    StatusCode makeBadDigiCol(Event::TkrDigiCol* pDigis);

    /// called by TrkCalibAlg to cause an update of the strip lists
    StatusCode update(CalibData::BadStrips* pDead, CalibData::BadStrips* pHot);

    
    /// queryInterface - required for a service
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);
    /// required for a service
    static const InterfaceID& interfaceID() { 
        return ITkrBadStripsSvc::interfaceID(); 
    }
    /// returns the service type
    const InterfaceID& type() const;   
    
private:
	  
    //void makeCol(const int size); // left over from original attempt
    
    /// reads bad strips from file file
    void readFromFile(std::ifstream* file);
	
    StatusCode doInit();

    StatusCode generateBadClusters() ;    	
	
    /// File name for constants
    std::string m_commonBadStripsFile;
    std::string m_badStripsFile[2];  
    
    /// array to hold bad strips vectors  [ max needed: 576 = 16*18*2 ]   
    mutable IndexedVector<stripCol> m_stripsCol[2];
	
	BadVisitor* m_visitor;

    ITkrGeometrySvc* m_tkrGeom;
	
    bool m_empty[2];

    Event::TkrDigiCol*      m_pBadDigi[2];
    Event::TkrClusterCol*   m_pBadClus[2];
    Event::TkrIdClusterMap* m_pBadMap[2];
    mutable bool m_generateBadClusters;
    mutable calibType m_calibType;
};

//! Fill the ASCII output stream

std::ostream& operator<<(std::ostream &s, stripCol* v) {
    int size = v->size();
    if (size) {
        s << " size " << size << std::endl << " strips " ;
        for (int j=0;j<size;j++) {
            int strip = (*v)[j].getStripNumber();
            s << strip<< " " ;
        }
        s << std::endl;
    }
    return s;
}

inline std::ostream& TkrBadStripsSvc::fillStream( std::ostream& s ) const 
{
    int nTowers = m_stripsCol[0].getDim(0);
    int nLayers = m_stripsCol[0].getDim(1);
    int nViews  = m_stripsCol[0].getDim(2);
    int tower, layer, view;

    s << "class TkrBadStripsSvc bad strip lists from fillStream: " << std::endl;
    for(tower=0; tower<nTowers;++tower) {
        for(layer=0; layer<nLayers; ++layer) {
            for(view=0; view<nViews; ++view) {
                s << "Twr/lyr/view (" << tower << ", " << layer << ", " << view  << ")";
                const stripCol* v = &m_stripsCol[m_calibType](tower, layer, view);
                s << v ;
            }
        }
    }
    return s;
}

#endif // TKRBADSTRIPSSVC_H


