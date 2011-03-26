/** @file ITkrHitTruncationTool.h
*/

/**
* @class ITkrHitTruncationTool
*
* @brief Interface to the tool that checks for truncation at fit time
*
* @author Leon Rochester
*/

#ifndef ITkrHitTruncATIONTOOL_H
#define ITkrHitTruncATIONTOOL_H

#include "GaudiKernel/IAlgTool.h"

#include "idents/TkrId.h"
#include "geometry/Vector.h"

static const InterfaceID IID_ITkrHitTruncationTool("ITkrHitTruncationTool", 1 , 0);

class ITkrHitTruncationTool : virtual public IAlgTool 
{
public:

    /// 
    virtual StatusCode analyzeDigis() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrHitTruncationTool; }

    virtual double getDistanceToTruncation(int tower, int tray, int face, int view,
                                           double localX) = 0;

    virtual double getDistanceToTruncation(int tower, int plane, Vector towerPos) = 0;

    virtual double getDistanceToTruncation(idents::TkrId tkrId, Vector towerPos) = 0;

};
#endif
