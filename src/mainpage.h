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
 * @section jobOptions jobOptions
 *
 * @param TkrFailureModeSvc.towerList
 * Provide a list of strings of the form "tower" indicating towers that will be made dead.
 * @param TkrFailureModeSvc.layerList
 * Provide a list of strings of the form "tower_layer_view" indicating 
 * planes that will be made dead.
 *
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

