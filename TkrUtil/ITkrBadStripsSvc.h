#ifndef __ITKRBADSTRIPSSVC_H
#define __ITKRBADSTRIPSSVC_H 1

#include "GaudiKernel/IInterface.h"
#include "idents/GlastAxis.h"

#include <string>
#include <vector>
#include <fstream>

//----------------------------------------------
//
//   TkrBadStripsSvc
//
//   Tracker BadStrips Service. Supplies the bad strips 
//   for use by SiClustersAlg, for example
//----------------------------------------------
//             Leon Rochester, 3-June-2001
//----------------------------------------------

static const InterfaceID IID_ITkrBadStripsSvc(907, 1 , 0); 

/// A small class to define tagged strips 
class TaggedStrip 
{
public: 
    enum {tagShift = 12, stripMask = 0xFFF, BIG = stripMask};
    
    TaggedStrip(int stripNumber = 0, int tag = 0)
    {
        m_stripNumber = (stripNumber & stripMask);
        m_tag = (tag & stripMask);
    }
    
    ~TaggedStrip() {}
    
    int getStripNumber() const { return m_stripNumber;}
    int getTag()         const { return m_tag;}
    bool isTaggedBad()   const { return getTag()>0;}
    static TaggedStrip makeLastStrip() { return TaggedStrip(BIG, BIG);}
    static TaggedStrip makeTaggedStrip(const int &strip) 
    {
        return TaggedStrip(strip);
    }
    
    operator int() const
    {
        return ((m_stripNumber & stripMask) << tagShift) 
            | ((m_tag & stripMask));
    }
    
private:
    int m_stripNumber;
    int m_tag;
};


typedef std::vector<TaggedStrip> stripCol;
typedef stripCol::iterator       stripCol_it;
typedef stripCol::const_iterator stripCon_it;


class ITkrBadStripsSvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    static const InterfaceID& interfaceID() { return IID_ITkrBadStripsSvc; }
   
    virtual int getIndex(int tower, int layer, 
        idents::GlastAxis::axis axis) const = 0;
    virtual const stripCol* getBadStrips(int tower, int layer, 
        idents::GlastAxis::axis axis) const = 0;
    virtual const stripCol* getBadStrips(int index) const = 0;
    virtual bool isBadStrip(int tower, int layer, 
        idents::GlastAxis::axis axis, int strip) const = 0;
    virtual bool isBadStrip(const stripCol* v, int strip) const = 0;
 };

#endif