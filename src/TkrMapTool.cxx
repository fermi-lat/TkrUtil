/**
* @class TkrMapTool
*
* @brief Implements a Gaudi Tool for setting the candidate track energies before 
*        the track fit
*
* @author The Tracking Software Group
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrMapTool.cxx,v 1.2.56.1 2012/01/20 01:55:28 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 

#include "TkrUtil/ITkrMapTool.h"

#include <iomanip>
#include <map>

class TkrMapTool : public AlgTool, virtual public ITkrMapTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrMapTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrMapTool() {}

    /// @brief calculations associated with the TEM diagnostic info

    StatusCode initialize();

private:
    int elecToGeo(int gtcc, int gtrc);
    int geoToElec(int plane, int end);
    
    void geoToElec(int plane, int end, int& gtcc, int& gtrc);    
    void elecToGeo(int gtcc, int gtrc, int& plane, int& end);
        
    int elecIndex(int gtcc, int gtrc) { return gtccMult*gtcc+gtrc; };

    std::map< int, int>  m_tkrMap;
    std::map< int, int>  m_tkrInverseMap;
};

static ToolFactory<TkrMapTool> s_factory;
const IToolFactory& TkrMapToolFactory = s_factory;

TkrMapTool::TkrMapTool(const std::string& type, 
                       const std::string& name, 
                       const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrMapTool>(this);

    return;
}

StatusCode TkrMapTool::initialize()
{
    // Purpose and Method: 
    // Inputs:  
    // Outputs: 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    setProperties();


    //This maps the electronic space into the geometrical one
    //Two cables map on to each plane, one at each end
    //
    // m_tkrMap[gtccMult*CC+RC] = PP; (gtccMult==10)
    //   CC = gtcc
    //   RC = gtrc
    //   PP = plane

    /*
    // this is the explicit mapping -- 
    m_tkrMap[ 0] = m_tkrMap[10] =  0;    m_tkrMap[20] = m_tkrMap[30] =  1;
    m_tkrMap[60] = m_tkrMap[70] =  2;    m_tkrMap[40] = m_tkrMap[50] =  3;
    m_tkrMap[ 1] = m_tkrMap[11] =  4;    m_tkrMap[21] = m_tkrMap[31] =  5;
    m_tkrMap[61] = m_tkrMap[71] =  6;    m_tkrMap[41] = m_tkrMap[51] =  7;
    m_tkrMap[ 2] = m_tkrMap[12] =  8;    m_tkrMap[22] = m_tkrMap[32] =  9;
    m_tkrMap[62] = m_tkrMap[72] = 10;    m_tkrMap[42] = m_tkrMap[52] = 11;
    m_tkrMap[13] = m_tkrMap[ 3] = 12;    m_tkrMap[23] = m_tkrMap[33] = 13;
    m_tkrMap[63] = m_tkrMap[73] = 14;    m_tkrMap[43] = m_tkrMap[53] = 15;
    m_tkrMap[14] = m_tkrMap[ 4] = 16;    m_tkrMap[24] = m_tkrMap[34] = 17;
    m_tkrMap[74] = m_tkrMap[64] = 18;    m_tkrMap[44] = m_tkrMap[54] = 19;
    m_tkrMap[15] = m_tkrMap[ 5] = 20;    m_tkrMap[25] = m_tkrMap[35] = 21;
    m_tkrMap[75] = m_tkrMap[65] = 22;    m_tkrMap[45] = m_tkrMap[55] = 23;
    m_tkrMap[ 6] = m_tkrMap[16] = 24;    m_tkrM[ap[26] = m_tkrMap[36] = 25;
    m_tkrMap[76] = m_tkrMap[66] = 26;    m_tkrMap[46] = m_tkrMap[56] = 27;
    m_tkrMap[17] = m_tkrMap[ 7] = 28;    m_tkrMap[37] = m_tkrMap[27] = 29;
    m_tkrMap[77] = m_tkrMap[67] = 30;    m_tkrMap[57] = m_tkrMap[47] = 31;
    m_tkrMap[ 8] = m_tkrMap[18] = 32;    m_tkrMap[28] = m_tkrMap[38] = 33;
    m_tkrMap[68] = m_tkrMap[78] = 34;    m_tkrMap[58] = m_tkrMap[48] = 35;
    */

    // this is the same as the above

    int i;

    for(i=0;i<nRc;++i) {
        m_tkrMap[i]    = m_tkrMap[10+i] = 4*i;
        m_tkrMap[20+i] = m_tkrMap[30+i] = 4*i+1;
        m_tkrMap[60+i] = m_tkrMap[70+i] = 4*i+2;
        m_tkrMap[40+i] = m_tkrMap[50+i] = 4*i+3;
    }

    int j;

    for (i=0; i<nRc; ++i) {
        for(j=0; j<nCc; ++j) {
            int ind = j*gtccMult + i;
            int plane = m_tkrMap[ind];
            int end   = endArray[j];
            m_tkrInverseMap[planeMult*plane+end] = ind;
            //std::cout << plane << " " << end << " " << ind << std::endl;
        }
    }

    return sc;
}

int TkrMapTool::elecToGeo(int gtcc, int gtrc)
{
    return m_tkrMap[elecIndex(gtcc, gtrc)];
}

int TkrMapTool::geoToElec(int plane, int end)
{
    return m_tkrInverseMap[planeMult*plane + end];
}

void TkrMapTool::geoToElec(int plane, int end, int& gtcc, int& gtrc)
{
    int res = geoToElec(plane, end);
    gtcc = res/gtccMult;
    gtrc = res%gtccMult;
}

void TkrMapTool::elecToGeo(int gtcc, int gtrc, int& plane, int& end)
{
    int res = elecToGeo(gtcc, gtrc);
    plane = res;
    end   = endArray[gtcc];
}
