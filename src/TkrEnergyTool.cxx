/**
 * @class TkrEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrEnergyTool.cxx,v 1.30 2005/06/21 23:29:34 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "TkrUtil/ITkrEnergyTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrEnergyTool.h"

#include "GlastSvc/Reco/IPropagator.h"

class TkrEnergyTool : public AlgTool, virtual public ITkrEnergyTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrEnergyTool() {}

    /// @brief Initialization of the tool
    StatusCode initialize();

    /// @brief Method to determine the total event energy given that from the Calorimeter
    ///        and the "best" track from the Tracker.
    double getTotalEnergy(const Event::TkrTrack* track, double calEnergy);

private:
    /// Internal methods
    /// Return Point given track parameters and deltaZ from track start
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Pointer to the cluster tool
    ITkrQueryClustersTool* m_clusTool;

    /// Pointer to the propagator
    IPropagator*           m_propagator;
};

static ToolFactory<TkrEnergyTool> s_factory;
const IToolFactory& TkrEnergyToolFactory = s_factory;

// constants defined at file scope
namespace {

    // Some constants collected from the file:
    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;
}

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrEnergyTool::TkrEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrEnergyTool>(this);

    return;
}

StatusCode TkrEnergyTool::initialize()
{
    // Purpose and Method: Initialization of the tool
    // Inputs:  none
    // Outputs: none, sets data members
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;

    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }

    //Locate a pointer to the G4Propagator
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagator)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }
    
    return sc;
}

double TkrEnergyTool::getTotalEnergy(const Event::TkrTrack* track, double CalEnergy)
{   
    // Purpose and Method: Augments calorimeter energy by that found in the tracker
    // Inputs:  The "best" track from Tracker, total Calorimeter energy
    // Outputs: The best estimate of the total energy of the event 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    const Event::TkrTrackHit* hit = track->front();

    int    topPlane = m_tkrGeom->getPlane((hit->getTkrId()));
    int    topLayer = m_tkrGeom->getLayer(topPlane);
    bool   isTop    = m_tkrGeom->isTopPlaneInLayer(topPlane);
    double arc_len  = 0.; 
    double convZ    = m_tkrGeom->getConvZ(m_tkrGeom->getLayer(hit->getTkrId()));
    Vector dir_ini  = hit->getDirection(Event::TkrTrackHit::SMOOTHED); 
    double secTheta = fabs(1./dir_ini.z());

    //Back up the start to the middle of the preceding converter
    // if the first hit is just under that converter
    Point  x_ini    = hit->getPoint(Event::TkrTrackHit::SMOOTHED);

    //Back up the start to the middle of the preceding converter
    // if the first hit is just under that converter
    // otherwise, just up enuf to get a run on the first plane
    double deltaZ = ( isTop ? hit->getZPlane() - convZ : 1.);
    x_ini   += dir_ini * deltaZ;
    arc_len += deltaZ * secTheta;

    double arc_tot = (x_ini.z() - m_tkrGeom->calZTop()) * secTheta; // to z at top of cal
    double addRad  = 0.;

    if(isTop) addRad = 0.5*m_tkrGeom->getRadLenConv(topLayer);
    
    // Do the swim
    m_propagator->setStepStart(x_ini, dir_ini);
    m_propagator->step(arc_tot);

    double radKal = m_propagator->getRadLength();                 

    // Set up summed var's and loop over all layers between track start and cal
    int    numHits[NTYPES];
    int    numLayers[NTYPES];

    for(int i=0; i < NTYPES; ++i) 
    {
        numHits[i] = 0;
        numLayers[i] = 0;
    }

    // need to loop over layers, because track can end before the end of the tracker
    int layer = topLayer+1; // so the "while" works
    double sprdMax = m_tkrGeom->trayWidth()/2.;

    while(layer--) 

    {
        //HepMatrix Q = m_propagator->getMscatCov(arc_len, CalEnergy/2.);
        //double xms = Q(1,1);
        //double yms = Q(3,3);

        // 4.0 sigma and not smaller then 2 mm (was 2.5 sigma)& less than a tower
        double xSprd = 80.*secTheta; //std::min(sqrt(4.+xms*16.), sprdMax); 
        double ySprd = 80.*secTheta; //std::min(sqrt(4.+yms*16.), sprdMax);

        // Assume location of shower center in given by 1st track
        Point x_hit = getPosAtZ(track, arc_len); 
 
        convType type = m_tkrGeom->getLayerType(layer);
        
        // count the layers and close hits by type
        numHits[type] += m_clusTool->numberOfHitsNear(layer, xSprd, ySprd, x_hit, dir_ini);
        numLayers[type]++;
        numLayers[ALL]++;

        // Increment arc-length
        if (layer > 0) 
        {
            int nextlayer = layer-1;
            deltaZ = m_tkrGeom->getLayerZ(layer)-m_tkrGeom->getLayerZ(nextlayer);
            arc_len += fabs(deltaZ*secTheta); 
        }
    }
    
    // Energy from nearby hit counting
    double ene_trks = _thinCoeff*numHits[STANDARD] 
                    + _thickCoeff*numHits[SUPER]
                    + _noradCoeff*numHits[NOCONV]; // Coefs are MeV/hit
 
    //Just the radiators
    // addRad removes half of the first converter, if present
    double thinConvRadLen  = m_tkrGeom->getAveConv(STANDARD);
    double thickConvRadLen = m_tkrGeom->getAveConv(SUPER);
    double rad_nom         = (thinConvRadLen*numLayers[STANDARD] 
                           + thickConvRadLen*numLayers[SUPER] - addRad) * secTheta;

    // The non-radiator stuff
    double trayRadLen = m_tkrGeom->getAveRest(ALL);
    double rad_min    = numLayers[ALL] * trayRadLen * secTheta;
    double radHits    = rad_nom + rad_min;

    double rad_swim = std::max(radKal, radHits);

    double ene_total = CalEnergy;
    // why are we normalizing by rad_nom instead of rad_nom+rad_min???
    if (rad_nom > 0.) ene_total += ene_trks * rad_swim / rad_nom;

    return ene_total;
}
