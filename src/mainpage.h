// Mainpage for doxygen

/** @mainpage package TkrUtil
 
  @author Leon Rochester
 
  @section description Description
 
  This package provides TKR utilities for common use, both inside and outside of TkrRecon.  
 
  @section TkrCalibAlg TkrCalibAlg
  TkrCalibAlg is called once per event, before any tracker code is executed. It 
  checks the calibration database for a calibration in the validity interval 
  for each calibration type, and
  if a new one is required (or it's the first time through) it calls the appropriate
  service(s) (see below), to cause them to update their lists.
 
  The algorithm's access methods for TkrBadStripsSvc and TkrFailureModeSvc
  are distinct from the usual ones and
  modify the class members, so for safety, are implemented
  through a separate abstract interface. This requires the services to check against
  two interface ID's in their queryInterface() methods.
 
  The default flavor for TkrCalibAlg is set to "ideal". In this mode the algorithm does
  nothing. So the flavor property of the algorithm must be set in the jobOptions for 
  anything to happen. In addition, there are separate properties for each calibration
  type, and setting these overrides the general flavor for that type.
  
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
 
  @section TkrBadStripsSvc TkrBadStripsSvc
  TkrBadStripsSvc creates a list of individual strips failures in the Tkr, and utilities
  to search the lists to allow digi and recon algorithms to deal with hits based
  on those lists.
 
  As is the case for TkrFailureMode Svc (q.v.), it takes its input from two sources: 
  an ASCII list 
  of (tower, plane) elements, each with its own list of badstrips (dead or hot); 
  and/or from the TkrBadStrips calibration data in the TCDS. 
 
  @section TkrSplitsSvc TkrSplitsSvc
  TkrSplitsSvc maintains a list of splits for the low and high controllers for each
  plane in the detector. The default split is after chip 11, or strip 767. 
  Otherwise, the splits come in from the calibration database.

  There's also a back-door to set the splits through an input xml file. This is for 
  quick development use and is not guarranteed to be maintained. If a file is
  specified, it overrides the database calibration.
 
  @section TkrToTSvc TkrToTSvc
  TkrToTSvc maintains the ToT thresholds and gains for each strip in the detector. The default
  settings are: all strips the same, and set to produce the same result as before this
  service was introduced. "randomized" mode internally generates a random set of thresholds and gains 
  similar to the those measured in the Engineering module. In this case, the default values
  are ignored. Will eventually be part of the calibration database
 
  @section TkrGeometrySvc TkrGeometrySvc
  TkrGeometrySvc assembles the methods required by various TKR algorithms that deal
  with TKR geometry.
 
  In addition, it stores pointers for TkrFailureModeSvc, TkrAlignmentSvc, and
  TkrBadStripsSvc, TkrSplitsSvc, TkrToTSvc, and GlastPropagatorSvc, which simplifies a lot of code, 
  since many modules that use
  the geometry use the other services as well.
 
  @section TkrQueryClustersTool TkrQueryClustersTool
  TkrQueryClustersTool allows algorithms to access certain quantities based on,
  but not immediately available in the TDS.
 
 @section jobOptions jobOptions
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
 *
  @param TkrFailureModeSvc.towerList
  Provide a list of strings of the form "tower" indicating towers that will be made dead.
  @param TkrFailureModeSvc.layerList
  Provide a list of strings of the form "tower_layer_view" indicating 
  planes that will be made dead.
 *
  @param TkrBadStripsSvc.badStripsFile
  The name of the file containing
  a list of bad (dead and hot) strips (default: null). These will be merged with bad strips
  coming from the calibration database.
 *
  @param TkrSplitsSvc.splitsFile
  The name of the xml file containing the splits specification -- not guaranteed to be maintained
 *
  @param TkrToTSvc.ToTFile
  The name of the file containing the thresholds and gains. Currently
  has no effect.
  @param TkrToTSvc.mode
  Current values are "EM" (to generate EM-like constants) and "ideal" (gives
  the default uniform constants below.
  @param TkrToTSvc.defaultThreshold
  For ideal mode, all thresholds are set to this value 
  (default = -2.92 microsec/fC).
 
  @param TkrToTSvc.defaultGain
  For ideal mode, all gains are set to this value (default = 2.50267833 microsec/fC)
  *
  @param TkrAlignmentSvc.simFile
  The name of the file containing the alignment constants to be used
  during digitization. Alignments of all elements down to the wafer may be specified.

  <br>
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
 
  @param TkrAlignmentSvc.recFile
  Currently no alignment is done during reconstruction. (But it will be soon!)
  The name of the file containing the alignment constants to be used
  during reconstruction. May be over-ridden by testMode, below. Format is the same
  as that of the simFile, except that we expect to generate constants only for the towers.
  @param TkrAlignment.maximumDelta
  Limits the maximum correction (default = 5 mm).  Possibly needed to protect against bad things
  happening during digitization, not to mention towers crashing together!
  
 
  <hr>
  @section notes release.notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements
  <hr>
 
 */

