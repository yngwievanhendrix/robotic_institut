/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volumeVertex.cpp
 *  \brief  Definition of class \c FM_volumeVertex
 *  \author Birgit Stender
 *  \date   23-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_volumeVertex.h"

using namespace FM;

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex constructor
/**
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Constructor.
 */
/*------------------------------------------------------------------------------*/
 FM_volumeVertex::FM_volumeVertex():
	GW_Vertex   (),
	pVolume_		( NULL )
{
	/* NOTHING */
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex destructor
/**
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Destructor.
 */
/*------------------------------------------------------------------------------*/
 FM_volumeVertex::~FM_volumeVertex()
{
	/* NOTHING */
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::SetVolume
/**
 *  \param  volume [FM_volume&] The parent volume we are managed by.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Set the parent volume.
 */
/*------------------------------------------------------------------------------*/
void FM_volumeVertex::SetVolume( FM_volume& volume )
{
	pVolume_ = &volume;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::GetVolume
/**
 *  \return [FM_volume*] The volume.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Set the volume we are managed by.
 */
/*------------------------------------------------------------------------------*/
FM_volume* FM_volumeVertex::GetVolume()
{
	return pVolume_;
}
	
/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVolume
/**
*  \return [GW_Face*] The volume.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Set the volume we are managed by.
*/
/*------------------------------------------------------------------------------*/
const FM_volume* FM_volumeVertex::GetVolume() const
{
	return pVolume_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::GetNumberNeighbor
/**
 *  \return [GW_U32] Answer.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Get the number of neighbors around this vertex.
 */
/*------------------------------------------------------------------------------*/
GW_U32 FM_volumeVertex::GetNumberNeighbor()
{
	GW_U32 nNum = 0;
	for( FM_volumeVertexIterator it = this->BeginVolumeVertexIterator(); it!=this->EndVolumeVertexIterator(); ++it )
		nNum++;
	return nNum;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::IsBoundaryVertex
/**
 *  \return [GW_U32] Answer.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Test if the vertex is a boundary one.
 */
/*------------------------------------------------------------------------------*/
/*
GW_Bool FM_volumeVertex::IsBoundaryVertex()
{
	for( FM_volumeVertexIterator it = this->BeginVolumeVertexIterator(); it!=this->EndVolumeVertexIterator(); ++it )
	{
		if( it.GetLeftVolume()==NULL || it.GetRightVolume()==NULL )
			return GW_True;
	}
	return GW_False;
}
*/

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::BeginVolumeIterator
/**
 *  \return [FM_volumeIterator] The iterator.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  Begin iterator on the surrounding of the vertex.
 */
/*------------------------------------------------------------------------------*/
FM_volumeIterator FM_volumeVertex::BeginVolumeIterator()
{
	if( this->GetVolume()==NULL )
	{
		return this->EndVolumeIterator();
	}
	else
	{
		if (this->GetVolume()->GetNextVertex(*this)!=NULL)
			return FM_volumeIterator( this, this);
		else
			return this->EndVolumeIterator();

	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::EndVolumeIterator
/**
 *  \return [FM_volumeIterator] The iterator.
 *  \author Birgit Stender
 *  \date   27-02-2011
 * 
 *  End iterator for the surrounding of the vertex.
 */
/*------------------------------------------------------------------------------*/
FM_volumeIterator FM_volumeVertex::EndVolumeIterator()
{
	return FM_volumeIterator(NULL,NULL);
}


/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::BeginVolumeVertexIterator
/**
*  \return [FM_volumeVertexIterator] The iterator.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Begin iterator on the surrounding of the vertex.
*/
/*------------------------------------------------------------------------------*/
FM_volumeVertexIterator FM_volumeVertex::BeginVolumeVertexIterator()
{
	if( this->GetVolume()==NULL )
	{
		return this->EndVolumeVertexIterator();
	}
	else
	{
		return FM_volumeVertexIterator( this, this->GetVolume()->GetNextVertex(*this));
	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeVertex::EndVolumeVertexIterator
/**
*  \return [FM_volumeVertexIterator] The iterator.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  End iterator for the surrounding of the vertex.
*/
/*------------------------------------------------------------------------------*/
FM_volumeVertexIterator FM_volumeVertex::EndVolumeVertexIterator()
{
	return FM_volumeVertexIterator(NULL,NULL);
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////