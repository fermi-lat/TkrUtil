/** @file ITkrReasonsTool.h
@brief Abs. int. for TkrReasonsTool
@author Leon Rochester
$Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/TkrUtil/ITkrReasonsTool.h,v 1.1 2011/02/24 22:04:50 lsrea Exp $
*/

#ifndef _H_ITkrReasonsTool
#define _H_ITkrReasonsTool

#include "GaudiKernel/IAlgTool.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrReasonsTool("ITkrReasonsTool", 0 , 0); 

/** @class ITkrReasonsTool
* @brief Abstract interface: evaluates the reasons for missing clusters
*/

class   ITkrReasonsTool : virtual public IAlgTool {
public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrReasonsTool; }
    virtual void setParams(const Point& end_pos, const int next_plane=-1) 
        const = 0;
    virtual bool isFailed() const = 0;
    virtual Vector getEdgeDistance() const = 0;
    virtual Vector getGapDistance() const = 0;
    virtual double getTruncDistance() const = 0;
    virtual double getBadClusterDistance() const = 0;
    virtual double getMinimumDistance(const Point& end_pos, 
        const int next_plane=-1) const = 0;

};

#endif  // _H_ITkrReasonsTool
