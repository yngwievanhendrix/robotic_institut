
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volumeIterator.h
 *  \brief  Definition of class \c FM_volumeIterator
 *  \author Birgit Stender
 *  \date   27-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_VOLUMEITERATOR_H_
#define _FM_VOLUMEITERATOR_H_

#include "../../gw/gw_core/GW_Config.h"


namespace FM {
	
	using namespace GW;

	class FM_volume; 
	class FM_volumeVertex;

	/*------------------------------------------------------------------------------*/
	/** 
	*  \class  FM_volumeIterator
	*  \brief  Iterator on a group of volumes surounding a vertex in a volume mesh.
	*  \author Birgit Stender
	*  \date   27-02-2011
	*
	*  Usefull iterator.
	*/ 
	/*------------------------------------------------------------------------------*/

	class FM_volumeIterator
	{
	public:

		FM_volumeIterator(  FM_volumeVertex* pStartVert, FM_volumeVertex* pVert, FM_volume* pVolume=NULL, GW_U32 nNbrIncrement = 0 );

		/* assignement */
		FM_volumeIterator& operator=( const FM_volumeIterator& it);

		/* evaluation */
		GW_Bool operator==( const FM_volumeIterator& it);
		GW_Bool operator!=( const FM_volumeIterator& it);

		/* indirection */
		FM_volume* operator*(  );

		/* progression : \todo take in acount NULL pointer */
		void operator++();

	private:

		FM_volumeVertex* pStartVert_;
		FM_volumeVertex* pVert_;
		FM_volume* pStartVolume_;
		FM_volume* pVolume_;

		/** just for debug purpose */
		GW_U32 nNbrIncrement_;
	};

} // End namespace FM

#endif // _FM_VOLUMEITERATOR_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////