
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volumeVerexIterator.h
 *  \brief  Definition of class \c FM_volumeVertexIterator
 *  \author Birgit Stender
 *  \date   27-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_VOLUMEVERTEXITERATOR_H_
#define _FM_VOLUMEVERTEXITERATOR_H_

#include "../../gw/gw_core/GW_Config.h"

namespace FM {

using namespace GW;

class FM_volume; 
class FM_volumeVertex;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  FM_volumeVertexIterator
 *  \brief  An iterator on the vertex around a given vertex.
 *  \author Birgit Stender
 *  \date   27-02-2011
 *
 *  To iterate on a volume vertex.
 */ 
/*------------------------------------------------------------------------------*/

class FM_volumeVertexIterator
{

public:

	FM_volumeVertexIterator(FM_volumeVertex* pStartVert, FM_volumeVertex* pVert, 
												 FM_volumeVertex* pVertNeighbor=NULL, GW_Bool pOnVolume = GW_True,	
												 GW_U32 nNbrIncrement = 0 );

	/* assignement */
	FM_volumeVertexIterator& operator=( const FM_volumeVertexIterator& it);

	/* evaluation */
	GW_Bool operator==( const FM_volumeVertexIterator& it);
	GW_Bool operator!=( const FM_volumeVertexIterator& it);

	/* indirection */
	FM_volumeVertex* operator*(  );

	/* progression */
	void operator++();

private:

	FM_volumeVertex* pStartVert_;
	FM_volumeVertex* pVert_;
	FM_volume* pVolume_;
	FM_volume* pStartVolume_;
	FM_volumeVertex* pVertNeighbor_;
	GW_Bool pOnVolume_;

	/** just for debug purpose */
	GW_U32 nNbrIncrement_;

};

} // End namespace FM

#endif // _FM_VOLUMEVERTEXITERATOR_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
