
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_propagationVolume.h
 *  \brief  Definition of class \c FM_propagationVolume
 *  \author Gabriel Peyré
 *  \date   4-12-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_propagationVolume_H_
#define _FM_PROPAGATIONVOLUME_H_

#include "../../gw/gw_core/GW_Config.h"
#include "../fm_core/FM_volume.h"

namespace FM {

	using namespace GW;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  FM_propagationVolume
 *  \brief  A volume to make fast marching computations.
 *  \author Birgit Stender
 *  \date   01-03-2011
 *
 *  Should contain propagation vertexes.
 */ 
/*------------------------------------------------------------------------------*/

class FM_propagationVolume:	public FM_volume
{

public:

	FM_propagationVolume();
	virtual ~FM_propagationVolume();

private:

};

} // End namespace FM

#endif // _FM_PROPAGATIONVOLUME_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
