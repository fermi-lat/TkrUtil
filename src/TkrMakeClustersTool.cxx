/** @file TkrMakeClustersTool.cxx
* @brief Tool to make TkrClusters, both good and bad

 @author Leon Rochester

 $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrMakeClustersTool.cxx,v 1.9 2009/09/09 00:25:54 lsrea Exp $
*/

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataSvc.h"
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
    StatusCode calculateToT();

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
        int nBad, int& rawToT, int& end) const;

    /// Keep pointer to the geometry service
    ITkrGeometrySvc*  m_tkrGeom;  
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* m_pBadStrips;
    /// Keep pointer to the ToT service
    ITkrToTSvc* m_pToT;
    /// Data service
    DataSvc* m_dataSvc;
    /// if STANDARDCLUSTERS, usual clustering; if BADCLUSTERS, construct bad-cluster list
    ITkrBadStripsSvc::clusterType m_type;
    TaggedStrip m_lastStrip;
    ///
    bool m_test250;
};

// Static factory for instantiation of algtool objects
//static ToolFactory<TkrMakeClustersTool> s_factory;
//const IToolFactory& TkrMakeClustersToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrMakeClustersTool);

// Standard Constructor
TkrMakeClustersTool::TkrMakeClustersTool(const std::string& type, 
                                         const std::string& name, 
                                         const IInterface* parent)
                                         : AlgTool( type, name, parent )
{    
    // Declare additional interface
    declareInterface<ITkrMakeClustersTool>(this); 
    declareProperty("test250", m_test250 = false);
}

StatusCode TkrMakeClustersTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    setProperties();

    m_tkrGeom = 0;
    if( serviceLocator() ) {
        sc = serviceLocator()->service( "TkrGeometrySvc", m_tkrGeom, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return sc;
        }
    }
    m_pBadStrips = m_tkrGeom->getTkrBadStripsSvc();
    m_pToT       = m_tkrGeom->getTkrToTSvc();
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
        return sc;
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

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

    m_type       = type;
    m_lastStrip  = TaggedStrip::makeLastStrip();

    unsigned int defaultStatus = m_tkrGeom->getDefaultClusterStatus();

    // keep this handy to do ToT plots
    //bool debugToT = false;
    //if(debugToT) {
    //    int i,j,k,l;
    //    int rawToT = 250;
    //    for(i=0;i<16;++i) {
    //        for(j=0;j<18;++j) {
    //            for(k=0;k<2;++k) {
    //                for(l=0;l<1536;++l) {
    //                    bool isBad = m_pBadStrips->isBadStrip(i,j,(idents::GlastAxis::axis)k,l);
    //                    double eDep = 0.5;
    //                    //int rawToT = m_pToT->getGain(eDep, i,j,k,l);
    //                    float x = m_pToT->getMuonScale(i,j,k,l);
    //                    if(!isBad) {
    //                        std::cout << x << " " ;
    //                        if((l+1)%64==0) std::cout << std::endl;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

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
                        ToT = calculateMips(pDigi, strip0, stripf, nBad, rawToT, end);
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
                                     int strip0, int stripf, int nBad, int& rawToT, int& end) const
{
    int layer      =  pDigi->getBilayer();
    idents::GlastAxis::axis view       =  pDigi->getView();
    int tower      = (pDigi->getTower()).id();

    int lastStrip = pDigi->getLastController0Strip();
    if(strip0<=lastStrip) {
        if (stripf<=lastStrip) end = 0; else end = 2;
    } else end = 1;

    if (end<2) { 
        rawToT = pDigi->getToT(end);
    } else {
        // bit of a kludge for when the cluster overlaps the splitPoint
        int tot0 = pDigi->getToT(0);
        int tot1 = pDigi->getToT(1);
        if(tot0<255&&tot1<255) {
            rawToT = (pDigi->getToT(0)*(lastStrip-strip0+1)
                + pDigi->getToT(1)*(stripf-lastStrip))/(stripf-strip0+1);
        } else {
            rawToT = 255;
        }
    }

    // 255 is not a real ToT; it just says that no ToT was measured
    double ToT;
    if(rawToT==255) { ToT = -1.;}
    else {
        typedef std::vector<double> mipsVec;
        typedef mipsVec::const_iterator mIter;
        unsigned int size = stripf-strip0+1;
        mipsVec totVec(size);                   
        double mips;
        int strip, i=0;                   
        // what is the flag for an invalid ToT?
        // Set the ToT to -1.
        // could use bad strip, but ToTSvc should know too?

        // Possibility:
        // if the strip is bad, set ToT to 1000.,
        // then it will be effectively skipped
        // if there are no stips left, ToT = -1;

        for (strip=strip0; strip<=stripf; ++strip, ++i) {
            int localRawToT = pDigi->getToTForStrip(strip);
            // for "heavy ion" test
            if(m_test250) localRawToT = 250; 
            if(nBad>0 && m_pBadStrips->isBadStrip(tower, layer, view, strip)) {
                mips = 1000.0; 
            } else {
                mips = m_pToT->getMipsFromToT(localRawToT, 
                    tower, layer, view, strip);
            }

            totVec[i] = mips;
        }

        // now pick the answer
        // 1 strip is easy, for 2, pick the lowest, otherwise pick lowest of 
        // the interior strips
        // We'll do better when we do this for real
        //
        // Why the lowest?
        // Because: The same amount of energy (ToT) deposited in two strips
        // with different gain will give a lower rawTot in the strip
        // with the lower gain, so it will be swamped by the rawTot from 
        // the other strip. The strip with the higher gain generates
        // the lower ToT.
        //

        // if the strip is bad, set ToT to 1000.,
        // then it will be effectively skipped
        // if there are no strips left, ToT = -1;
        if (stripf==strip0) {
            ToT = totVec[0];
        } else if (stripf-strip0==1) {                       
            ToT = std::min(totVec[0], totVec[1] );
        } else {
            ToT = *std::min_element(++totVec.begin(), --totVec.end() );
        }
    }

    if(ToT==1000.) ToT = -1.;

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
StatusCode TkrMakeClustersTool::calculateToT()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    //get the clusters
    SmartDataPtr<Event::TkrClusterCol> 
        clusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);
    // Recover a pointer to the raw digi objects
    SmartDataPtr<Event::TkrDigiCol> digiCol(m_dataSvc,
        EventModel::Digi::TkrDigiCol);
    if(digiCol==0) return StatusCode::SUCCESS;
    // map the clusters
    std::map<int, Event::TkrDigi*> digiMap;
    int index;
    Event::TkrDigi* digi = 0;
    int tower, layer, view;
    Event::TkrDigiCol::const_iterator digiIter;
    for(digiIter=digiCol->begin();digiIter!=digiCol->end();++digiIter) {
        digi = *(digiIter);
        idents::TowerId towerId = digi->getTower();
        tower = towerId.id();
        layer = digi->getBilayer();
        view  = digi->getView();      
        index = 1000*tower + 10*layer + view;
        digiMap[index] = digi;
    }
    if(clusterCol==0||digiCol==0) return sc;
    int clusSize = clusterCol->size();
    int i;
    int end;
    int lastIndex = -1;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        int strip0 = clus->firstStrip();
        int stripf = clus->lastStrip();
        int rawToT = (int)clus->getRawToT();
        int nBad   = clus->getNBad();
        tower = clus->tower();
        layer  = clus->getLayer();
        int plane  = clus->getPlane();
        view   = m_tkrGeom->getView(plane);
        index = 1000*tower + 10*layer + view;
        if(index!=lastIndex) {
            digi = digiMap[index];
            lastIndex = index;
        }
        double ToT = calculateMips(digi, strip0, stripf, nBad, rawToT, end);
        clus->setMips(ToT);
    }
    return sc;
}
