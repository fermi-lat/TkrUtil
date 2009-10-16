/// @file TkrTowerBits.h

/** @class TkrTowerBits
@brief handles the calculation of the TKR trigger

@author Leon Rochester
$Header$
*/

#ifndef TkrTowerBits_h
#define TkrTowerBits_h

class TkrTowerBits
{
public:
    TkrTowerBits () :
      m_xBits(0), m_yBits(0), m_debug(false) {}
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
      void setBit(int layer, int view) {
          (view==0 ? setXBit(layer) : setYBit(layer));
      }

      void setXBit(int layer)  {m_xBits |= (1<<layer);}

      void setYBit(int layer)  {m_yBits |= (1<<layer);}

      bool isTriggered() {
          bool ret = true;
          // make up layer bits, loop, return true if true
          unsigned int layerBits = m_xBits&m_yBits;
          int i;
          if(m_debug) std::cout << "tower trigger/x/y/l " <<std::oct << m_xBits << " " 
              << m_yBits << " " << layerBits << std::dec << std::endl;
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

      void clear() { m_xBits = m_yBits = 0; }
      void setDebug(bool debug) { m_debug = debug;}
   
private:
      int    m_xBits;
      int    m_yBits;
      bool   m_debug;
};

typedef std::vector<TkrTowerBits*> towerVec;
typedef towerVec::iterator towerIter;

#endif

