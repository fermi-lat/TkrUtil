// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/TkrMeritTool.cxx,v 1.1 2003/01/16 23:46:11 lsrea Exp $

// Include files

//#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "GlastSvc/Reco/IReconTool.h"

#include <algorithm>
//#include <numeric>
//include <cmath>
//#include <string>

class TkrMeritTool : public AlgTool, virtual public IReconTool 
{
public:
    
    TkrMeritTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrMeritTool() { }

    StatusCode initialize();
    
    /// gets called to return variables
    StatusCode get(const std::string& name, double&);
    
private:
    StatusCode calculate();
    
    /// Finds the two best tracks associated with this vertex
    void doTracks(Event::TkrVertex& vertex);
    /// Calculates the suplus hit ratio
    StatusCode doExtraHits(
        const Event::TkrVertex& vertex,
        const Event::TkrClusterCol& clusters, double energy);
    /// Calculates the x and y intercepts with the plane of the skirt
    void doSkirt();
    /// Returns the radiation lengths before each plane
    double radLen ( int layer) const;
    
    /// run number at last call
    int m_lastRun;
    /// event number at last call
    int m_lastEvent;
    /// Status code returned at last call
    StatusCode m_lastSc;
    
    // from TDS
    // from track analysis
    
    /// number of vertices in this event
    double m_nVertices; 
    /// best and next-best track in the chosen vertex
    const Event::TkrFitTrackBase* m_thePair[2];
    
    /// point of intercept with plane of skirt
    Point m_skirtIntercept;     
    /// surplus hit ratio
    double m_surplus_hit_ratio;

    // some pointers to services
    
    /// pointer to tracker geometry
    ITkrGeometrySvc*  m_pGeom;
    /// pointer to event data service
    IDataProviderSvc* m_pEventSvc;
    /// pointer to tool service
    IToolSvc* m_pToolSvc;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrMeritTool> s_factory;
const IToolFactory& TkrMeritToolFactory = s_factory;


// Standard Constructor
TkrMeritTool::TkrMeritTool(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent)
                           : AlgTool( type, name, parent ), 
                           m_lastRun(-1), m_lastEvent(-1),
                           m_lastSc(StatusCode::FAILURE)
{    
    // Declare additional interface
    declareInterface<IReconTool>(this); 
}

StatusCode TkrMeritTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    if( serviceLocator() ) {
        sc = serviceLocator()->service( "TkrGeometrySvc", m_pGeom, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find EventSvc" << endreq;
            return sc;
        }
        sc = serviceLocator()->service( "EventDataSvc", m_pEventSvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find EventSvc" << endreq;
            return sc;
        }
        sc = serviceLocator()->service( "ToolSvc", m_pToolSvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find ToolSvc" << endreq;
            return sc;
        }
    }
    log << MSG::INFO << "TkrMeritTool successfully initialized" << endreq;
    return sc;
}


StatusCode TkrMeritTool::get(const std::string& name, double& x)
{
    // Returns requested variable

    // It only does the calculation once for all the calls in the same event

    // This is set up so that if the calculation fails, the input variable
    // is unchanged... that is, it can be set to a default value by the
    // caller before the call, and will remain at that value if the 
    // name is unrecognized, or the calculation fails

    // Get out fast if we don't recognize the name
    if (name!="test" 
        && name!="REC_Surplus_Hit_Ratio" 
        && name!="REC_Tkr_SkirtX"
        && name!="REC_Tkr_SkirtY") {
        return StatusCode::FAILURE;
    }
    
    StatusCode sc = StatusCode::SUCCESS;
    
    // simple check to see if anything is working
    if(name=="test") {
        x = 99;
        return sc;
    }

    sc = calculate();

    // if calculation fails, return without touching input variable
    
    if (sc==StatusCode::SUCCESS) {
        if (name=="REC_Surplus_Hit_Ratio") 
        {
            x = m_surplus_hit_ratio;
        } 
        else if (name=="REC_Tkr_SkirtX") 
        {
            x = m_skirtIntercept.x();
        } 
        else if (name=="REC_Tkr_SkirtY") 
        {
            x = m_skirtIntercept.y();
        } 
    }    
    return sc;
} 

StatusCode TkrMeritTool::calculate() 
{   
    double _minEnergy = 30.0;
    
    //MsgStream log(msgSvc(), name());
    
    
    // We do all the calculations the first time we're called for each event
    // and just return for subsequent calls

    SmartDataPtr<Event::EventHeader>   header      (m_pEventSvc, 
        EventModel::EventHeader);
    
    if (header->run()==m_lastRun && header->event()==m_lastEvent) {
        return m_lastSc;
    }

    m_lastEvent = header->event();
    m_lastRun   = header->run();
    m_lastSc = StatusCode::FAILURE;

    //m_surplus_hit_ratio = 0.0;
    //m_skirtIntercept = Point(0., 0., 0.);
    
    SmartDataPtr<Event::TkrVertexCol>  vertices    (m_pEventSvc, 
        EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::CalClusterCol> clusters    (m_pEventSvc, 
        EventModel::CalRecon::CalClusterCol);
    SmartDataPtr<Event::TkrClusterCol> tkrClusters (m_pEventSvc, 
        EventModel::TkrRecon::TkrClusterCol);
    
    /*
    SmartDataPtr<Event::McParticleCol> particles   (m_pEventSvc, 
        EventModel::MC::McParticleCol);
    SmartDataPtr<Event::CalXtalRecCol> xtalrecs    (m_pEventSvc, 
        EventModel::CalRecon::CalXtalRecCol);
    SmartDataPtr<Event::AcdRecon>      acdrec      (m_pEventSvc, 
        EventModel::AcdRecon::Event);
        */

    if (vertices==0 || tkrClusters==0) {
        return m_lastSc;
    }

    m_nVertices = vertices->size();   
    if (m_nVertices==0) return m_lastSc;

    double energy = 0.0;
    if (clusters) {
        if( clusters->num()>0 ){
            energy = (clusters->getCluster(0))->getEnergySum();
        }
    }
    if (energy < _minEnergy) energy = _minEnergy;

    
    Event::TkrVertex& vertex = *(vertices->front());
    
    doTracks(vertex);
    if (m_thePair[0]==0) return m_lastSc;
    
    if (doExtraHits(vertex, tkrClusters, energy).isFailure()) return m_lastSc;
    
    doSkirt();
    
    m_lastSc = StatusCode::SUCCESS;
    return m_lastSc;
}

void TkrMeritTool::doTracks(Event::TkrVertex& vertex) 
{
    // params of the "best" and "pair" track
    //  best is [0], pair is [1]
    
    SmartRefVector<Event::TkrFitTrackBase>::const_iterator tkBegin 
        = vertex.getTrackIterBegin();
    SmartRefVector<Event::TkrFitTrackBase>::const_iterator tkEnd   
        = vertex.getTrackIterEnd();
    SmartRefVector<Event::TkrFitTrackBase>::const_iterator tkIter;

    int numTracks = vertex.getNumTracks();
    
    Event::TkrFitTrackBase* skip = 0;
    for (int pass = 0; pass<std::min(2,numTracks); pass++) {
        m_thePair[pass] = 0;
        SmartRefVector<Event::TkrFitTrackBase>::const_iterator trackIter 
            = vertex.getTrackIterBegin();
        double tkrQual = -1;
        for(tkIter=tkBegin; tkIter!=tkEnd; tkIter++)
        {
            const Event::TkrFitTrackBase* track = *trackIter;
            if (pass==1 && skip==track) continue;// skip the best track 2nd time;
            if (track->getQuality() > tkrQual) {
                tkrQual = track->getQuality();
            m_thePair[pass] = track;
            skip = const_cast<Event::TkrFitTrackBase*> (track);
            }
        }
    }
}

// replace this with an array in TkrGeometrySvc, maybe one for 
//  just the converter, and another for the total

double TkrMeritTool::radLen(int layer) const {
    double radLenEff = 0.03; 
    if (layer>11) radLenEff = 0.18;
    if (layer>15) radLenEff = 0.0;
    radLenEff += 2*0.4/90.36 + 0.007;  // silicon + stuff
    return radLenEff;
}

StatusCode TkrMeritTool::doExtraHits(const Event::TkrVertex& vertex,
                               const Event::TkrClusterCol& /* clusters */,
                               double energy)
{    
    // Hardwired constants for this method
    const double _maxHRF            = 3.0;
    const double _coneSigmas        = 5.0;
    const double _maxSprd           = 500.;     // mm
    const double _sigSprd           = 2.021;    // = 3.5 sigma/sqrt(3)
    // for multiple scattering formula
    const double _msEnergyCoeff     = 14.0; //MeV
    const double _msLogCoeff        = 0.038;
    // for last-layer calculation
    const double _llEnergyCoeff = 10.0;  // MeV
    const double _llConstCoeff  = 4.0;
    const double _llLogCoeff    = 2.0;
    
    ITkrQueryClustersTool* pQuery;
    StatusCode sc = m_pToolSvc->retrieveTool("TkrQueryClustersTool", pQuery);
    if( sc.isFailure() ) {
        return StatusCode::FAILURE;
    }
    
    double CsICorrEnergy = energy;
    
    int norma = 1;
    
    // params of the vertex
    Point x0  = vertex.getPosition();
    Vector t0 = vertex.getDirection();
    
    // get the best track
    
    const Event::TkrFitTrackBase* track  = m_thePair[0];
    Vector t1 = track->getDirection();
    double tZ = fabs(t1.z());
    
    /*
    // a possible refinement: take larger angle of two tracks
    const Event::TkrFitTrack* track2 = m_thePair[1];   
    Vector t2(0,0,1);
    if (track2) {
    t2 = track2->getDirection();
    double tZ = std::min(fabs(t1.z()),fabs(t2.z()));
    }
    */
        
    // find the number of hits around the first hit
    
    int firstLayer = vertex.getLayer();
    Point firstHit  = vertex.getPosition();
    
    double zProj = fabs(t0.z());
    // no idea what this is meant to be...
    double hitRadFact = std::min(zProj + sqrt(1.-zProj*zProj)/zProj, _maxHRF);
    
    double radLenEff = radLen(firstLayer)/tZ;    
    double thetaCone = _msEnergyCoeff/(CsICorrEnergy/2.)*sqrt(radLenEff)
        *(1. + _msLogCoeff*log(radLenEff));
    
    double minSprd = _coneSigmas* thetaCone * m_pGeom->trayHeight();
    double dltSprd = minSprd;
    double norm    = sqrt(1-zProj*zProj);
    double xFact = 1, yFact = 1;
    if (norm>0.0001) {
        xFact = 1. + (hitRadFact - 1.)*fabs(t0.x())/norm;
        yFact = 1. + (hitRadFact - 1.)*fabs(t0.y())/norm;
    }
    double dltX = dltSprd*xFact;
    double dltY = dltSprd*yFact;
    
    /*
    std::cout << "rl " << radLenEff << " thcone " << thetaCone 
    << " dsprd " << dltSprd << std::endl;
    std::cout << "nrm " << norm << " x/yfact " << xFact << " " << yFact
    << " dlx/y " << dltX << " " <<  dltY << std::endl;
    */
    
    
    // get the hits int the first layer
    // dltX and dltY are too large  The hits in the first layer
    // have very small lever arm for separating!
    double Rec_Sum_Hits = pQuery->numberOfHitsNear(firstLayer, 
        .5*dltX, .5*dltY, firstHit);
    //double RSH1= Rec_Sum_Hits;
    
    double Rec_showerHits1 = Rec_Sum_Hits;
    double Rec_showerHits2 = 0;
    
    if(firstLayer > 11) {
        Rec_showerHits1 = 0; 
        Rec_showerHits2 = Rec_Sum_Hits;
    }    
    if(Rec_Sum_Hits < 2.) Rec_Sum_Hits = 2.;
    
    int    lastLayer 
        = (int)(_llConstCoeff 
        + _llLogCoeff*log(CsICorrEnergy/_llEnergyCoeff)  
        + firstLayer);
    lastLayer = std::min(m_pGeom->numLayers()-1, lastLayer);
    if(lastLayer - firstLayer < 5 
        && lastLayer == m_pGeom->numLayers()-1) Rec_Sum_Hits +=.5;
    
    //  Find the number of hits around the rest of the gamma trajectory
    double outHits = 0;
    double deltaS    = m_pGeom->trayHeight()/zProj;
    double disp      =  deltaS;
    
    xFact *= sqrt(t0.z() *t0.z() + t0.x()*t0.x());
    yFact *= sqrt(t0.z() *t0.z() + t0.y()*t0.y());
    
    //double RSH2= Rec_Sum_Hits;
    
    // now do the rest
    
    // only the layers between first and last get included in the count
    // for Rec_Sum_Hits. The rest are also included in the Shower_Hits
    for(int layer=firstLayer+1; layer<m_pGeom->numLayers(); layer++) {
        double thetaMS = _msEnergyCoeff/(CsICorrEnergy/2.)*sqrt(radLenEff)
            *(1.+ _msLogCoeff*log(radLenEff));
        double xSprd = std::min(thetaMS * disp * xFact * _sigSprd, _maxSprd); 
        double ySprd = std::min(thetaMS * disp * yFact * _sigSprd, _maxSprd);
        Point trkHit((Point)(disp*t0 + x0));
        double nearHits = pQuery->numberOfHitsNear(layer, xSprd, 
            ySprd, trkHit);
            /*
            std::cout << "lyr " << layer << " nhits " << nearHits << std::endl
            << "          x/ySprd " << xSprd <<" " <<ySprd << std::endl
            << "          rl, disp " << radLenEff <<" "<< disp << std::endl;
            */
        if(layer < 12) Rec_showerHits1 +=  nearHits;
        else      Rec_showerHits2 +=  nearHits;
        
        if( layer <= lastLayer) {
            Rec_Sum_Hits +=     nearHits;
            outHits +=  pQuery->numberOfHitsNear(layer, _maxSprd, _maxSprd,   
                trkHit) - nearHits;
        }
        disp += deltaS;
        radLenEff += radLen(layer)/tZ;
    }
    
    norma = lastLayer - firstLayer;
    
    m_surplus_hit_ratio = Rec_Sum_Hits / norma;
    //double Rec_Outside_Hit_Ratio = outHits / norma;
    
    /*
    std::cout << "en " << CsICorrEnergy 
    << "  dltX,Y " << dltX << " " << dltY << std::endl;
    std::cout << "first/last layer " << firstLayer << " " << lastLayer 
    << " norm " << norma
    << " RSH1,2,x " << RSH1 << " " << RSH2 << " " << Rec_Sum_Hits
    << std::endl;   
    std::cout << "results: surplus hit ratio " << m_surplus_hit_ratio 
    << " outside hit ratio " << Rec_Outside_Hit_Ratio
    << std::endl;
    */
    
    return StatusCode::SUCCESS;  
}

void TkrMeritTool::doSkirt() 
{
    // top of CsI, as determined from Event Display
    const double _skirtZ = -26.475; // skirt z in mm
    
    const Event::TkrFitTrackBase* track = m_thePair[0];
    Point endPoint = track->getPosition(Event::TkrFitTrackBase::End);
    Vector endDir  = track->getDirection(Event::TkrFitTrackBase::End);
    double deltaS  = (endPoint.z() - _skirtZ)/endDir.z();
    m_skirtIntercept = endPoint - deltaS*endDir;
    return;
}

