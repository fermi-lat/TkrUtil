// Mainpage for doxygen

/** @mainpage package TkrUtil
 
  @author Leon Rochester
 
  @section description Description
 
  This package provides TKR utilities for common use, both inside and outside of TkrRecon.

  @section TkrCalibAlg TkrCalibAlg
  TkrCalibAlg can be called once or twice for each event. If called once, it is called 
  before any tracker code is executed. It 
  checks the calibration database for a calibration in the validity interval 
  for each calibration type, and
  if a new one is required (or it's the first time through) it calls the appropriate
  service(s) (see below), to cause them to update their lists.

  It can also be called twice, once at the beginning of the simulation sequence
  and a second time at the beginning of the reconstruction sequence. In this case, the two calls
  are distinguished by having different names in the Gaudi sequence:
  TkrCalibAlg/SimCalib and TkrCalibAlg/RecCalib. Each instance of the algorithm can be
  initialized separately, and controls an independent set of calibration constants.
 
  Those services that store local information about the calibrations (TkrBadStripsSvc, 
  TkrFailureModeSvc, and eventually, TkrAlignmentSvc), have been modified to be able to
  store two sets of internal information. Which set is active is controlled by a call to 
  SetCalibType() for that service, which call is provided automatically by TkrCalibAlg.

  The algorithm's access methods for TkrBadStripsSvc and TkrFailureModeSvc
  are distinct from the usual ones and
  modify the class members, so for safety, are implemented
  through a separate abstract interface. This requires the services to check against
  two interface ID's in their queryInterface() methods.
 
  The default flavor for TkrCalibAlg is set to "ideal". In this mode the algorithm does
  nothing. So the flavor property of the algorithm must be set in the jobOptions for 
  anything to happen. In addition, there are separate properties for each calibration
  type, and setting these overrides the general flavor for that type. 

  @section TkrBadStripsSvc TkrBadStripsSvc
  TkrBadStripsSvc creates a list of individual strips failures in the Tkr, and utilities
  to search the lists to allow digi and recon algorithms to deal with hits based
  on those lists.
 
  As is the case for TkrFailureMode Svc (q.v.), it takes its input from two sources: 
  an ASCII list 
  of (tower, plane) elements, each with its own list of badstrips (dead or hot); 
  and/or from the TkrBadStrips calibration data in the TCDS. 
 
  @section TkrFailureModeSvc TkrFailureModeSvc
  TkrFailureModeSvc creates a list of large-scale failures in the Tkr, 
  and utilities to search the lists to allow digi and recon algorithms 
  to deal with hits based on those lists.
 
  It takes input from two sources: a list of towers and/or a list of (tower, layer, view);
  or, from the TkrBadStrips calibration data in the TCDS, through a visitor. The first is
  for informal or one-off tests, and the second is ultimately for actual calibration data 
  from the detector. In case both sources are present, the input from each is merged into 
  a common list, with duplicates removed.
 
  It provides a method to see if a TKR plane is contained in the lists.
 
  @section TkrAlignmentSvc
  TkrAlignmentSvc provides alignment constants for simulation and reconstruction. The constants
  are read in when filenames are supplied. (See below.) This will eventually become a part of
  the calibration system.
 
  @section TkrGeometrySvc TkrGeometrySvc
  TkrGeometrySvc assembles the methods required by various TKR algorithms that deal
  with TKR geometry.
 
  In addition, it stores pointers for TkrFailureModeSvc, TkrAlignmentSvc, and
  TkrBadStripsSvc, TkrSplitsSvc, TkrToTSvc, and GlastPropagatorSvc, which simplifies 
  a lot of code, since many modules that use the geometry use the other services as well.


  @section TkrMakeClustersTool TkrMakeClustersTool

  TkrMakeClustersTool either makes clusters out of the digis in an event, taking into account
  dead strips, or (at initialization time) produces the list of bad clusters, 
  for use in the pattern recognition.
 
  @section TkrQueryClustersTool TkrQueryClustersTool

  TkrQueryClustersTool allows algorithms to access certain quantities based on TkrClusterCol,
  but not immediately available in the TDS. 

  @section TkrSplitsSvc TkrSplitsSvc
  TkrSplitsSvc maintains a list of splits for the low and high controllers for each
  plane in the detector. The default split is after chip 11, or strip 767. 
  Otherwise, the splits come in from the calibration database.

  There's also a back-door to set the splits through an input xml file. This is for 
  quick development use and is not guarranteed to be maintained. If a file is
  specified, it overrides the database calibration.

  This service also handles ReadController and CableController truncation. See parameters, below.
 
  @section TkrToTSvc TkrToTSvc
  TkrToTSvc maintains the ToT thresholds and gains for each strip in the detector. The default
  settings are: all strips the same, and set to produce the same result as before this
  service was introduced. "randomized" mode internally generates a random set of thresholds and gains 
  similar to the those measured in the Engineering module. In this case, the default values
  are ignored. Will eventually be part of the calibration database
 
 @section jobOptions jobOptions
  *
  @param TkrAlignmentSvc.simFile
  The name of the file containing the alignment constants to be used
  during digitization. Alignments of all elements down to the wafer may be specified. (See below.)
 
  @param TkrAlignmentSvc.recFile
  The name of the file containing the alignment constants to be used
  during reconstruction.  Format is the same as that of the simFile.
  @param TkrAlignmentSvc.testMode
  sets a fixed offset in x and y into every silicon wafer as a simple test
  @param TkrAlignment.maximumDelta
  Limits the maximum correction (default = 5 mm).  Possibly needed to protect against bad things
  happening during digitization, not to mention towers crashing together!
  
  *
  @param TkrBadStripsSvc.badStripsFile
  The name of the file containing
  a list of bad (dead and hot) strips (default: null). These will be merged with bad strips
  coming from the calibration database.

  The above property can be set individually for the simulation and reconstruction, if the 
  sequence is set up appropriately. in this case, the jobOptions looks like:

  TkrBadStripsSvc.simBadStripsFile = ...;

  TkrBadStripsSvc.recBadStripsFile = ...;
  *
  @param TkrGeometrySvc.siResolutionFactor
  This number multiplies the strip with to generate the intrinsic strip resolution. Default is 
  1.0/sqrt(12.0)
  @param TkrGeometrySvc.layerSeparation
  Used to determine whether adjacent planes belong to the same layer. If they are
  separated by less than this distance, they are in the same layer, otherwise not. Default is 8   mm.
  @param radLenCut
  Used to separate standard from superglast layers. Default is 0.10 radiation lengths.
  @param activeOffset
  The test propagator is launched with this offset from the edge of a tower. Default is 40 mm.
  *
  @param TkrCalibAlg.calibFlavor
  Sets the overall flavor of calibration requested. Defaults to "ideal", which does nothing.
  @param TkrCalibAlg.deadStripsCalibFlavor
  Overrides the overall flavor for dead strips
  @param TkrCalibAlg.hotStripsCalibFlavor
  Overrides the overall flavor for hot strips
  @param TkrCalibAlg.splitsCalibFlavor
  Overrides the overall flavor for splits
  @param TkrCalibAlg.chargeInjectionCalibFlavor
  Overrides the overall flavor for ToT charge injection calibration
  @param TkrCalibAlg.muonCalibFlavor
  Overrides the overall flavor for ToT muon calibration

  As discussed above, if TkrCalibAlg is being called twice, the correct form is:

  TkrSimCalib.... and TkrRecCalib....
 *
  @param TkrFailureModeSvc.towerList
  Provide a list of strings of the form "tower" indicating towers that will be made dead.
  @param TkrFailureModeSvc.layerList
  Provide a list of strings of the form "tower_layer_view" indicating 
  planes that will be made dead.
  *
  The above property can be set individually for the simulation and reconstruction, if the 
  sequence is set up appropriately. in this case, the jobOptions looks like:

  TkrFailureModeSvc.simTowerList = ...;

  TkrFailureModeSvc.recLayerList = ...;

  etc.

  @param TkrQueryClustersTool.towerFactor
  factor to multiply tower pitch to get test distance for the unmeasured direction
  *
  @param TkrSplitsSvc.splitsFile
  The name of the xml file containing the splits specification -- not guaranteed to be maintained
  @param TkrSplitsSvc.cableBufferSize
  The size of the cable controller buffer (Default is 128.)
  @param TkrSplitsSvc.defaultMaxStrips
  The size of the read controller buffer, if all are the same (Default is 64.)
  @param TkrSplitsSvc.maxStripsFile
  The name of a file containing the sizes of each individual read controller buffer. This information
  will end up either in the calibrations database or the conditions database. An example of the current
  implementation is src/test/maxStrips.xml.
 *
  @param TkrToTSvc.ToTFile
  The name of the file containing the thresholds and gains. Currently
  has no effect.
  @param TkrToTSvc.defaultThreshold
  For ideal mode, all thresholds are set to this value (default = 1.177 fC).
 
  @param TkrToTSvc.defaultGain
  For ideal mode, all gains are set to this value (default = 0.589 fC/microsecond)
  @param TkrToTSvc.defaultQuad
  For ideal mode, all quad terms are set to this value (default = 0.00490 (fC/microsecond)**2)
  @param TkrToTSvc.defaultMuonScale
  Scaling derived from calibration with cosmic muons (default = 1.0)
  @param TkrToTSvc.mevPerMip
  Most probable energy deposit for a vertical Mip (default = 0.113 MeV)
  @param TkrToTSvc.fCPerMip
  charge measured for most probable energy deposit (somewhat arbitrary, set to 5.0 fC)
  @param TkrToTSvc.countsPerMicrosecond
  internal clock ticks per microsecond (default = 5.0)
  @param TkrToTSvc.maxToT
  Maximum ToT value allowed
  @param TkrToTSvc.useSingleTowerConstants
  use a single set of ToT constants for all the towers (For tests, default = false)
  @param TkrToTSvc.baseTower
  If useSingleTowerConstants is true, the tower of the constants to be used
  @param TkrToTSvc.useDefaultIfMissing (default = false)
  Specifies treatment of missing values when using calibration data. If false,
  missing constants are set to zero; if true, standard default values are used.
  
  <hr>
  @section files Structure of Input Files
  The structure of an alignment file is as follows:
  <br>
 @verbatim
     //  comment  (all such lines are ignored)
                  (as are all blank lines)
     // one line for each element with non-zero constants
     // order is deltaX deltaY deltaZ rotX rotY rotZ
     // deltas in microns, rots in mrads.
     
     // elements are TOWER, TRAY, FACE, LADDER, WAFER
     // numbering follows the hardware convention
     
     // an element is specified by specifying the elements above it
     // any element can be entered as ELEMENT NUMBER or
     //                               ELEMENT NUMBER const1 const2 ... const6
     // In each case the transformation is with respect to the element
     // immediately above it in the hierarchy

     // Examples:
     
     // rotate tray 3 in tower 6 by 0.5 mrad, everything else nominal
     TOWER 6 
     TRAY  3     0.    0.  0.  0.     0.    0.5

     // rotate tower 11, and translate several ladders
     // here all the elements in the tower will be transformed
     // because the tower itself is transformed
     TOWER 11   50. -100. 25.  0.15  -0.05  0.7
     TRAY   5
     FACE   0
     LADDER 1   25.   3.   0.    0.   0.    0.
     LADDER 2   12.  -15.  0.    0.   0.    0.
     TRAY   9
     LADDER 3    0.   5.   0.    0.   0.    0.

 @endverbatim


  @section notes release.notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements
  <hr>
 
 */

