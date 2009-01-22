/**
 * @class ITkrGhostTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrGhostTool.h,v 1.2 2008/12/01 19:44:34 lsrea Exp $
 */
#ifndef ITkrGhostTool_h
#define ITkrGhostTool_h

#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ITkrGhostTool("ITkrGhostTool", 2 , 0);

class TkrTowerBits
{
public:
    TkrTowerBits () :
      m_xBits(0), m_yBits(0){}
      virtual ~TkrTowerBits() {}

      void setBit(Event::TkrCluster* clus) {
          //get layer, view, set appropriate bit
          int layer = clus->getLayer();
          idents::TkrId id  = clus->getTkrId();
          int view = id.getView();
          (view==0 ? setXBit(layer) : setYBit(layer));
      }
      void setBit(Event::TkrDigi* digi) {
          //get layer, view, set appropriate bit
          int layer = digi->getBilayer();
          int view = digi->getView();
          (view==0 ? setXBit(layer) : setYBit(layer));
      }

      void setXBit(int layer)  {m_xBits |= (1<<layer);}

      void setYBit(int layer)  {m_yBits |= (1<<layer);}

      bool isTriggered() {
          bool ret = true;
          // make up layer bits, loop, return true if true
          unsigned int layerBits = m_xBits&m_yBits;
          int i;
          int mask0 = 7;
          for(i=0;i<16;++i) {
              int mask = (mask0<<i);
              //std::cout << i << std::hex << (mask<<i) << std::dec << std::endl;
              if((layerBits&mask)==mask) return true;
          }
          return false;
      }

      unsigned int getXBits() {return m_xBits;}
      unsigned int getYBits() {return m_yBits;}
      unsigned int getBits()  {return m_xBits&m_yBits;}

      // these are the layers that participated in the software trigger
      unsigned int getTriggeredBits() {
          unsigned int layerBits = m_xBits&m_yBits;
          unsigned int trigBits  = 0;
          int i;
          int mask0 = 7;
          for(i=0;i<16;++i) {
              int mask = (mask0<<i);
              //std::cout << i << std::hex << (mask<<i) << std::dec << std::endl;
              if((layerBits&mask)==mask) trigBits |=mask;
          }
          return trigBits;
      }     
   
private:
      int    m_xBits;
      int    m_yBits;
};

typedef std::vector<TkrTowerBits*> towerVec;
typedef towerVec::iterator towerIter;


class ITkrGhostTool : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrGhostTool; }

    virtual StatusCode getTkrVector(unsigned short& tkrVector) = 0;
    virtual StatusCode calculateTkrVector(
        Event::TkrClusterCol* pCol, unsigned short& towerBits) = 0;
    virtual StatusCode calculateTkrVector(
        Event::TkrDigiCol* pCol, unsigned short& towerBits) = 0;
    virtual StatusCode flagSingles()   = 0;
    virtual StatusCode flagEarlyHits(Event::TkrClusterCol* col=0) = 0;
    virtual StatusCode flagEarlyTracks() = 0;
    virtual StatusCode flagEarlyVertices() = 0;

};

#endif
