// Mainpage for doxygen

/** @mainpage package TkrUtil
 *
 * @author Leon Rochester
 *
 * @section description Description
 *
 * This package provides utilities for common use in the TKR.    
 *
 * @section TkrFailureModeSvc TkrFailureModeSvc
 * TkrFailureModeSvc creates a list of large-scale failures in the Tkr, 
 * and utilities
 * to search the lists to allow digi and recon algorithms to deal with hits based
 * on those lists.
 *
 * It can take either a list of towers or a list of (tower, layer, view)
 * to create the lists of dead objects. It provides a method to see if a 
 * TKR plane is contained in the lists.
 *
 * @section TkrBadStripsSvc TkrBadStripsSvc
 * TkrBadStripsSvc creates a list of individual strips failures in the Tkr, 
 * and utilities
 * to search the lists to allow digi and recon algorithms to deal with hits based
 * on those lists.
 *
 * It currently takes an ASCII list of (tower, plane) elements, each with its
 * own list of badstrips (dead or hot). 
 *
 * @section TkrGeometrySvc TkrGeometrySvc
 * TkrGeometrySvc assembles the methods required by various TKR algorithms that deal
 * with TKR geometry.
 *
 * @section TkrQueryClustersTool TkrQueryClustersTool
 * TkrQueryClustersTool allows algorithms to access certain quantities based on,
 * but not immediately available in the TDS.
 *
 * @section TkrMeritTool TkrMeritTool
 * TkrMeritTool consolidates various calculations to generate various quantities
 * required by merit. It is accessed through the abstract interface IReconTool,
 * which is found in the Recon package.
 *
 * @section jobOptions jobOptions
 *
 * @param TkrFailureModeSvc.towerList
 * Provide a list of strings of the form "tower" indicating towers that will be made dead.
 * @param TkrFailureModeSvc.layerList
 * Provide a list of strings of the form "tower_layer_view" indicating 
 * planes that will be made dead.
 *
 * @param TkrBadStripsSvc.badStripsFile
 * The name of the file containing
 * a list of bad (dead and hot) strips. Will eventually be supplanted by
 * the calibration database.
 *
 * @param TkrAlignmentSvc.simFile
 * The name of the file containing the alignment constants to be used
 * during digitization. May be over-ridden by testMode, below.
 * @param TkrAlignmentSvc.recFile
 * The name of the file containing the alignment constants to be used
 * during reconstruction. May be over-ridden by testMode, below.
 * @param TkrAlignmentSvc.testMode
 * A non-zero testMode means that the alignment files will be generated
 * internally, rather than through external files.  Each bit in the word represents
 * a separate testMode, and multiple modes can be applied at the same time.
 * Currently:<br>  <br>
 * testMode bit 0:   shift all wafers in simulation by a fixed (dX, dY)
 * <br>
 * testMode bit 1:   shift all wafers in reconstruction by a fixed (dX, dY)
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 *
 */

