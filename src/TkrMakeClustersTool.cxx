/** @file TkrMakeClustersTool.cxx
* @brief Tool to make TkrClusters, both good and bad

 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrMakeClustersTool.cxx,v 1.5 2006/11/02 19:34:48 lsrea Exp $
*/

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrMakeClustersTool.h"

#include <vector>
#include <map>
#include "geometry/Point.h"  


class TkrMakeClustersTool : public AlgTool, virtual public ITkrMakeClustersTool 
{
public:

    enum clusterType { STANDARDCLUSTERS, BADCLUSTERS };

    TkrMakeClustersTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~TkrMakeClustersTool() { }

    StatusCode initialize();

    StatusCode makeClusters(
        Event::TkrClusterCol* pClus, Event::TkrIdClusterMap* clusMap,
        Event::TkrDigiCol* pTkrDigiCol,
        ITkrBadStripsSvc::clusterType clType=ITkrBadStripsSvc::STANDARDCLUSTERS);

private:

    /// gets the position of a cluster 
    Point position(int tower, int ilayer, int v, int strip0, int stripf) const;
    /// returns true if the two hits have a gap between them
    bool isGapBetween(const TaggedStrip &lowHit, const TaggedStrip &highHit) const;
    /// returns true if the cluster is "good"
    bool isGoodCluster( const TaggedStrip &lowHit, 
        const TaggedStrip &highHit, int nBad) const;

    /// get the list of bad strips
    const stripCol* getBadStrips(int tower, int digiLayer, 
        int view) const;
    /// makes a first guess at the corrected ToT for a cluster
    float calculateMips(Event::TkrDigi* pDigi, int strip0, int stripf, 
        int& rawToT, int& end) const;

    /// Keep pointer to the geometry service
    ITkrGeometrySvc*  m_tkrGeom;  
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* m_pBadStrips;
    /// Keep pointer to the ToT service
    ITkrToTSvc* m_pToT;
    /// if STANDARDCLUSTERS, usual clustering; if BADCLUSTERS, construct bad-cluster list
    ITkrBadStripsSvc::clusterType m_type;
    TaggedStrip m_lastStrip;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrMakeClustersTool> s_factory;
const IToolFactory& TkrMakeClustersToolFactory = s_factory;

// Standard Constructor
TkrMakeClustersTool::TkrMakeClustersTool(const std::string& type, 
                                         const std::string& name, 
                                         const IInterface* parent)
                                         : AlgTool( type, name, parent )
{    
    // Declare additional interface
    declareInterface<ITkrMakeClustersTool>(this); 
}

StatusCode TkrMakeClustersTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    m_tkrGeom = 0;
    if( serviceLocator() ) {
        sc = serviceLocator()->service( "TkrGeometrySvc", m_tkrGeom, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return sc;
        }
    }

    log << MSG::INFO << "TkrMakeClustersTool successfully initialized" << endreq;
    return sc;
}

StatusCode TkrMakeClustersTool::makeClusters(
    Event::TkrClusterCol* pClus, Event::TkrIdClusterMap* clusMap,
    Event::TkrDigiCol* pTkrDigiCol, ITkrBadStripsSvc::clusterType type)
{
    // Purpose: Makes Clusters from TkrDigis
    // Method:  Digis are scaned and grouped into contiguous groups
    // Inputs:  Digis, pointers to geometry and badstrips services
    // Outputs:  Clusters
    // Dependencies: None
    // Restrictions and Caveats:  None

    m_pBadStrips = m_tkrGeom->getTkrBadStripsSvc();
    m_pToT       = m_tkrGeom->getTkrToTSvc();
    m_type       = type;
    m_lastStrip  = TaggedStrip::makeLastStrip();

    unsigned int defaultStatus = m_tkrGeom->getDefaultClusterStatus();



    //Initialize the cluster lists...
    pClus->clear();

    Event::TkrDigiCol::const_iterator ppDigi = pTkrDigiCol->begin();
    int nclusters = 0;  // for debugging

    for (; ppDigi!= pTkrDigiCol->end(); ppDigi++) {
        // each digi contains the digitized hits from one layer of one tower
        Event::TkrDigi* pDigi = *ppDigi;

        int layer      =  pDigi->getBilayer();
        int view       =  pDigi->getView();
        int tower      = (pDigi->getTower()).id();

        int  towerX    = pDigi->getTower().ix();
        int  towerY    = pDigi->getTower().iy();
        int  tray      = 0;
        int  botTop    = 0;
        int  measure   = view == idents::GlastAxis::X 
            ? idents::TkrId::eMeasureX : idents::TkrId::eMeasureY;

        m_tkrGeom->layerToTray(layer, view, tray, botTop);

        idents::TkrId hitId(towerX, towerY, tray, (botTop==1), measure);

        // debug: std::cout << "digi t/l/v " << tower << " " << layer << " " << view << std::endl;
        // copy the hits, and make them into TaggedStrips

        int nHits = pDigi->getNumHits();

        stripCol stripHits(nHits);

        std::transform(pDigi->begin(), pDigi->end(), stripHits.begin(), 
            TaggedStrip::makeTaggedStrip);
        std::sort(stripHits.begin(), stripHits.end()); //paranoia

        // get the list of bad strips; pointer is zero if no list or empty list
        const stripCol* badStrips = getBadStrips(tower, layer, view);
        int badStripsSize = 0;
        if (badStrips) badStripsSize = badStrips->size();

        // now make the combined list
        nHits += badStripsSize;
        stripCol mergedHits(nHits+1);  // leave room for sentinel
        if (badStrips) {
            std::merge(stripHits.begin(), stripHits.end(), badStrips->begin(), badStrips->end(),
                mergedHits.begin());
        } else {
            std::copy(stripHits.begin(), stripHits.end(), mergedHits.begin());
        }
        mergedHits[nHits] = m_lastStrip; // big and bad; end of loop

        // the first strip of the current potential cluster
        TaggedStrip lowStrip  = *mergedHits.begin();  
        // the last strip of the current cluster
        TaggedStrip highStrip = lowStrip;       
        // the next strip
        TaggedStrip nextStrip = lowStrip;       
        int nBad = 0;
        bool kept;  // for debugging

        // Loop over the rest of the strips building clusters enroute.
        // Keep track of bad strips.
        // Loop over all hits, except the sentinel, which is there 
        // to provide a gap

        // don't let the loop go to the end... the code looks one ahead!
        stripCol_it itLast = mergedHits.end();
        itLast--;
        for( stripCol_it it = mergedHits.begin(); it!=itLast; ) {
            if(nextStrip.isTaggedBad()) nBad++;
            // here we get the next hit, and increment the iterator
            // at the same time!
            // debug: std::cout << "this pointer " << it <<" next " << it+1 << " end " << mergedHits.end() << std::endl;
            nextStrip = *(++it);
            // debug: std::cout << " got past next!" << std::endl;

            //If we have a gap, then make a cluster
            if (isGapBetween(highStrip, nextStrip)) {
                // there's a gap... see if the current cluster is good...
                // debug: std::cout << std::endl << "Test Cluster: " << lowStrip.getStripNumber() << " "
                //       << highStrip.getStripNumber() << " " << nBad ;
                if (kept = (isGoodCluster(lowStrip, highStrip, nBad))) {
                    // debug: std::cout << "good!" << std::endl;
                    // it's good... make a new cluster
                    int strip0 = lowStrip.getStripNumber();
                    int stripf = highStrip.getStripNumber();
                    Point pos = position(tower, layer, view, strip0, stripf);

                    // code to generate 1st order corrected ToT
                    int end;
                    int rawToT = 0;
                    float ToT = 0.0;
                    if(m_type!=ITkrBadStripsSvc::BADCLUSTERS) {
                        ToT = calculateMips(pDigi, strip0, stripf, rawToT, end);
                    }
                    unsigned int status = defaultStatus | 
                        ((end<<Event::TkrCluster::shiftEND)&Event::TkrCluster::maskEND);

                    Event::TkrCluster* cl = new Event::TkrCluster(hitId, strip0, stripf, 
                        pos, rawToT, ToT, status, nBad);


                    // for tests
                    //if(m_type == ITkrBadStripsSvc::BADCLUSTERS) {
                    //    std::cout << tower << " " << nHits << " " <<
                    //        tray << " " << botTop << " " << strip0 << " " << stripf-strip0+1 << std::endl;
                    //}
                    pClus->push_back(cl);
                    nclusters++;
                    if(clusMap!=0) (*clusMap)[hitId].push_back(cl);
                    // for tests
                    //std::cout << rawToT << " " << ToT << std::endl;
                } 
                lowStrip = nextStrip;  // start a new cluster with this strip
                nBad = 0;
            }
            // debug: std::cout << "on to next strip..." << std::endl;
            highStrip = nextStrip; // add strip to this cluster
        }
    }
    return StatusCode::SUCCESS;
}

float TkrMakeClustersTool::calculateMips(Event::TkrDigi* pDigi, 
                                     int strip0, int stripf, int& rawToT, int& end) const
{
    int layer      =  pDigi->getBilayer();
    int view       =  pDigi->getView();
    int tower      = (pDigi->getTower()).id();

    int lastStrip = pDigi->getLastController0Strip();
    if(strip0<=lastStrip) {
        if (stripf<=lastStrip) end = 0; else end = 2;
    } else end = 1;

    if (end<2) { 
        rawToT = pDigi->getToT(end);
    } else {
        // bit of a kludge for when the cluster overlaps the splitPoint
        rawToT = (pDigi->getToT(0)*(lastStrip-strip0+1)
            + pDigi->getToT(1)*(stripf-lastStrip))/(stripf-strip0+1);
    }

    typedef std::vector<double> mipsVec;
    typedef mipsVec::const_iterator mIter;
    unsigned int size = stripf-strip0+1;
    mipsVec totVec(size);                   
    double mips;
    int strip, i=0;                   
    // what is the flag for an invalid ToT?
    // could use bad strip, but ToTSvc should know too

    for (strip=strip0; strip<=stripf; ++strip, ++i) {
        int rawToT = pDigi->getToTForStrip(strip);
        mips = m_pToT->getMipsFromToT(rawToT, 
            tower, layer, view, strip);              
        totVec[i] = mips;
    }

    // now pick the answer
    // 1 strip is easy, for 2, pick the lowest, otherwise pick lowest of 
    // the interior strips
    // We'll do better when we do this for real

    double ToT;
    if (stripf==strip0) {
        ToT = totVec[0];
    } else if (stripf-strip0==1) {                       
        ToT = std::min(totVec[0], totVec[1] );
    } else {
        ToT = *std::min_element(++totVec.begin(), --totVec.end() );
    }
    return ToT;
}

Point TkrMakeClustersTool::position(int tower, int layer, int v,
                                int strip0, int stripf) const

{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  layer, view, first and last strip and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None

    // this converts from recon numbering to physical numbering of layers.
    // since we allow clusters to cross gaps (for "bad" clusters) 
    // we need to take average of position of first and last strip

    double firstStrip = strip0;
    double lastStrip  = stripf;
    HepPoint3D p = m_tkrGeom->getStripPosition(tower, layer, 
        v, firstStrip);
    p += m_tkrGeom->getStripPosition(tower, layer, 
        v, lastStrip);
    p *= 0.5;
    Point p1(p.x(), p.y(), p.z());
    return p1;

}

bool TkrMakeClustersTool::isGapBetween(const TaggedStrip &lowStrip, 
                                   const TaggedStrip &highStrip) const
{
    // Purpose and Method: decides whether there is a gap between two strips
    // Inputs:  strip numbers
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None

    //Get the actual hit strip number from the tagged strips
    int lowHit  = lowStrip.getStripNumber();
    int highHit = highStrip.getStripNumber();

    // gap between hits
    if (highHit > (lowHit + 1)) { return true; }

    // edge of chip
    int nStrips = m_tkrGeom->ladderNStrips();
    if(m_type==ITkrBadStripsSvc::STANDARDCLUSTERS 
        && (lowHit/nStrips) < (highHit/nStrips)) {return true; }

        return false;
}


bool TkrMakeClustersTool::isGoodCluster(const TaggedStrip &lowStrip, 
                                    const TaggedStrip &highStrip, 
                                    int nBad) const
{
    // Purpose and Method: Finds out if a cluster is "good"
    // Inputs: first and last strip, and number of bad strips
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None

    // always keep the cluster if we're constructing the bad clusters
    if (m_type==ITkrBadStripsSvc::BADCLUSTERS && highStrip!=m_lastStrip) {
        return true;
    }
    //Get the actual hit strip number from the tagged strips
    int lowHit  = lowStrip.getStripNumber();
    int highHit = highStrip.getStripNumber();

    // Require at least 1 good hit in the cluster    
    if ((highHit-lowHit+1)<=nBad) return false;

    // Require 10 or fewer bad hits in the cluster
    if (nBad>10)    return false;
    return true;
}

const stripCol* TkrMakeClustersTool::getBadStrips(int tower, int digiLayer, 
                                              int view) const
{
    // Purpose and Method: get the list of bad strips, if it exists
    // Inputs: tower, digiLayer, view
    // Outputs:  pointer to strip list, zero if list not there, or empty
    // Dependencies: None
    // Restrictions and Caveats:  None

    const stripCol* badStrips = 0;
    if (m_pBadStrips) {
        badStrips = 
            m_pBadStrips->getBadStrips(tower, digiLayer, 
            static_cast<idents::GlastAxis::axis>(view) );
        if (badStrips && badStrips->size()<=0) {
            badStrips = 0;
        }
    }
    return badStrips;
}

