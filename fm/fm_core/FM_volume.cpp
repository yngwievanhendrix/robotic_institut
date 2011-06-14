/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volume.cpp
 *  \brief  Definition of class \c FM_volume
 *  \author Birgit Stender
 *  \date   22-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_volume.h"

namespace FM {

	using namespace GW;

/*------------------------------------------------------------------------------*/
// Name : FM_volume constructor
/**
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Constructor.
 */
/*------------------------------------------------------------------------------*/
 FM_volume::FM_volume()
{
	for( GW_U32 i=0; i<4; ++i )
	{
		vertex_[i] = NULL;
		volumeNeighbors_[i] = NULL;
	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::virtual ~FM_volume
/**
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Destructor.
 */
/*------------------------------------------------------------------------------*/
 FM_volume::~FM_volume()
{
	for( GW_U32 i=0; i<4; ++i )
		GW_SmartCounter::CheckAndDelete( vertex_[i] );
}

 /*------------------------------------------------------------------------------*/
// Name : FM_volume::operator=
/**
 *  \param  volume [FM_volume&] volume
 *  \return [FM_volume&] *this
 *  \author Birgit Stender
 *  \date   03-03-2011
 * 
 *  Copy operator
 */
/*------------------------------------------------------------------------------*/
FM_volume& FM_volume::operator=(const FM_volume& volume)
{
	this->nID_ = volume.nID_;
	return *this;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::SetVolumeNeighbor
/**
 *  \param  pVolume [FM_Volume*] The volume. Can be NULL for surface volumes.
 *  \param  nFaceNum [GW_U32] The number of the face we are sharing with this volume.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Set one of the neigbhoring volumes.
 */
/*------------------------------------------------------------------------------*/

 void FM_volume::SetVolumeNeighbor(FM_volume* pVolume, GW_U32 nVertexNum)
{
	GW_ASSERT( nVertexNum<4 );
	volumeNeighbors_[nVertexNum] = pVolume;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::SetVolumeNeighbor
/**
 *  \param  pVolume0 [FM_volume*] volume 0
 *  \param  pVolume1 [FM_volume*] volume 1
 *  \param  pVolume2 [FM_volume*] volume 2
 *  \param  pVolume3 [FM_volume*] volume 3
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Assign the 4 neighboring volumes at once
 */
/*------------------------------------------------------------------------------*/
 void FM_volume::SetVolumeNeighbor(FM_volume* pVolume0, FM_volume* pVolume1, FM_volume* pVolume2, FM_volume* pVolume3)
{
	this->SetVolumeNeighbor( pVolume0, 0 );
	this->SetVolumeNeighbor( pVolume1, 1 );
	this->SetVolumeNeighbor( pVolume2, 2 );
	this->SetVolumeNeighbor( pVolume3, 3 );
}

/*------------------------------------------------------------------------------*/
// Name : FM_Volume::GetVolumeNeighbor
/**
 *  \param  nEdgeNum [GW_U32] number of the volume.
 *  \return [FM_volume*] The volume.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Get one of the volume neighbors
 */
/*------------------------------------------------------------------------------*/
 FM_volume* FM_volume::GetVolumeNeighbor(GW_U32 nFaceNum)
{
	GW_ASSERT( nFaceNum<4 );
	return volumeNeighbors_[nFaceNum];
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVolumeNeighbor
/**
*  \param  nEdgeNum [GW_U32] number of the volume.
*  \return [GW_Face*] The volume.
*  \author Birgit Stender
*  \date   22-02-2011
* 
*  Get one of the volume neighbor
*/
/*------------------------------------------------------------------------------*/
 const FM_volume* FM_volume::GetVolumeNeighbor(GW_U32 nFaceNum) const
{
	GW_ASSERT( nFaceNum<4 );
	return volumeNeighbors_[nFaceNum];
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVolumeNeighbor
/**
 *  \param  Vert [FM_volumeVertex&] The vertex.
 *  \return [GW_volume*] The volume. Can be NULL if correct vertex not found.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Return the neighboring volume opposite to the given vertex.
 */
/*------------------------------------------------------------------------------*/
 FM_volume* FM_volume::GetVolumeNeighbor(const FM_volumeVertex& vert)
{
	for( GW_U32 i=0; i<4; ++i )
	{
		if( vertex_[i]==&vert )
			return volumeNeighbors_[i];
	}
	return NULL;
}


/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVolumeNeighbor
/**
*  \param  vert0 [FM_volumeVertex&] first vertex of the face.
*  \param  vert1 [FM_volumeVertex&] second one.
*  \param  vert2 [FM_volumeVertex&] third one.
*  \return [FM_volume*] The  volume. Can be NULL if correct opposite vertex not found.
*  \author Birgit Stender
*  \date   22-02-2011
* 
*  Return the neighboring volume opposite to the given vertex.
*/
/*------------------------------------------------------------------------------*/
 FM_volume* FM_volume::GetVolumeNeighbor( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 )
{
	GW_I32 nVertex = this->GetVertexNumber(vert0, vert1, vert2);
	if( nVertex<0 )
		return NULL;
	return this->GetVolumeNeighbor( nVertex );
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::SetVertex
/**
 *  \param  Vert [FM_volumeVertex&] The new vertex.
 *  \param  nNum [GW_U32] Its number.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Set the vertex.
 */
/*------------------------------------------------------------------------------*/
 void FM_volume::SetVertex(FM_volumeVertex& vert, GW_U32 nNum)
{
	GW_ASSERT( nNum<4 );
	/* check the previous one, delete it if needed */
	GW_SmartCounter::CheckAndDelete( vertex_[nNum] );

	if( vert.GetVolume()==NULL )
		vert.SetVolume( *this );

	vertex_[nNum] = &vert;
	vert.UseIt();
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::SetVertex
/**
 *  \param  vert0 [FM_volumeVertex&] vertex 0
 *  \param  vert1 [FM_volumeVertex&] vertex 1
 *  \param  vert2 [FM_volumeVertex&] vertex 2
 *  \param  vert3 [FM_volumeVertex&] vertex 3
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Set all vertexes at once.
 */
/*------------------------------------------------------------------------------*/
void FM_volume::SetVertex(FM_volumeVertex& vert0, FM_volumeVertex& vert1, FM_volumeVertex& vert2, FM_volumeVertex& vert3)
{
	this->SetVertex( vert0, 0 );
	this->SetVertex( vert1, 1 );
	this->SetVertex( vert2, 2 );
	this->SetVertex( vert3, 3 );
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertex
/**
 *  \param  nNum [GW_U32] The number of the vertex.
 *  \return [FM_volumeVertex*] the vertex.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Get one of the vertex.
 */
/*------------------------------------------------------------------------------*/
FM_volumeVertex* FM_volume::GetVertex(GW_U32 nNum)
{
	GW_ASSERT( nNum<4 );
	return vertex_[nNum];
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertex
/**
*  \param  nNum [GW_U32] The number of the vertex.
*  \return [FM_volumeVertex*] the vertex.
*  \author Birgit Stender
*  \date   22-02-2011
* 
*  Get one of the vertexes.
*/
/*------------------------------------------------------------------------------*/
const FM_volumeVertex* FM_volume::GetVertex(GW_U32 nNum) const
{
	GW_ASSERT( nNum<4 );
	return vertex_[nNum];
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertex
/**
 *  \param  vert0 [FM_volumeVertex&] first vertex of the face.
 *  \param  vert1 [FM_volumeVertex&] second one.
 *  \param  vert2 [FM_volumeVertex&] third one.
 *  \return [FM_volumeVertex*] the vertex. can be NULL if face doesn't exist.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Get the vertex opposite to a given face
 */
/*------------------------------------------------------------------------------*/
FM_volumeVertex* FM_volume::GetVertex( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 )
{
	GW_I32 nVertex = this->GetVertexNumber( vert0, vert1, vert2 );
	if( nVertex<0 )
		return NULL;
	return this->GetVertex( nVertex );
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetNextVertex
/**
 *  \param  Vert [FM_volumeVertex&] The vertex.
 *  \return [FM_volumeVertex*] The next one.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Get the vertex following the one given by the user.
 */
/*------------------------------------------------------------------------------*/
FM_volumeVertex* FM_volume::GetNextVertex( const FM_volumeVertex& vert, GW_Bool direction )
{
	for( GW_U32 i=0; i<4; ++i )
	{
		if( vertex_[i]==&vert )
			if (direction)
				return vertex_[(i+1)%4];
			else
				return vertex_[(i-1)%4];

	}
	return NULL;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertexNumber
/**
 *  \param  volume [FM_volume&] The neighboring volume.
 *  \return [GW_I32] The number. Return -1 if not found.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Get the number of the vertex adjacent to the neighboring volume given.
 */
/*------------------------------------------------------------------------------*/
GW_I32 FM_volume::GetVertexNumber( const FM_volume& volume )
{
	for( GW_U32 i=0; i<4; ++i )
	{
		if( volumeNeighbors_[i]==&volume )
			return i;
	}
	return -1;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertexNumber
/**
 *  \param  vert [FM_volumeVertex&] The vertex.
 *  \return [GW_I32] The number. Return -1 if not found.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Return the number of the vertex.
 */
/*------------------------------------------------------------------------------*/
GW_I32 FM_volume::GetVertexNumber( const FM_volumeVertex& vert )
{
	for( GW_U32 i=0; i<4; ++i )
	{
		if( vertex_[i]==&vert )
			return i;
	}
	return -1;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetVertexNumber
/**
 *  \param  vert0 [FM_volumeVertex&] vertex 1
 *  \param  vert1 [FM_volumeVertex&] vertex 2
 *  \param  vert2 [FM_volumeVertex&] vertex 3
 *  \return [GW_I32] The number of the fourth vertex. Return -1 if not found.
 *  \author Birgit Stender
 *  \date   22-02-2011
 * 
 *  Return the number of the face corresponding to the three given 
 *  vertexes.
 */
/*------------------------------------------------------------------------------*/
GW_I32 FM_volume::GetVertexNumber( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 )
{
	for( GW_U32 i=0; i<4; ++i )
	{
		if( vertex_[i]==&vert0 )
		{
			if( vertex_[(i+1)%3]==&vert1 )
				if ( vertex_[(i+2)%3]==&vert2 )
					return (i+3)%3;
				if ( vertex_[(i+3)%3]==&vert2 )
					return (i+2)%3;
			if( vertex_[(i+2)%3]==&vert1 )
				if ( vertex_[(i+3)%3]==&vert2 )
					return (i+1)%3;
				if ( vertex_[(i+1)%3]==&vert2 )
					return (i+3)%3;
			if( vertex_[(i+3)%3]==&vert1 )
				if ( vertex_[(i+1)%3]==&vert2 )
					return (i+2)%3;
				if ( vertex_[(i+2)%3]==&vert2 )
					return (i+1)%3;
		}
	}
	return -1;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::IsResponsibleFor
/**
 *  \param  vert [FM_volumeVertex&] The vertex.
 *  \return [GW_Bool] yes/no ?
 *  \author Birgit Stender
 *  \date   02-03-2011
 * 
 *  Return true if we are the parent of a given vertex.
 */
/*------------------------------------------------------------------------------*/
 GW_Bool FM_volume::IsResponsibleFor( const FM_volumeVertex& vert )
{
	for (GW_U32 i=0; i<4; ++i)
	{
		if (vertex_[i]!=NULL)
			if (vertex_[i]==&vert)
				return vertex_[i]->GetVolume()==this;
	}
	return GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::IsIncluded
/**
 *  \param  vert [FM_volumeVertex&] The vertex.
 *  \return [GW_Bool] yes/no ?
 *  \author Birgit Stender
 *  \date   07-03-2011
 * 
 *  Return true if the given vertex is included in the volume.
 */
/*------------------------------------------------------------------------------*/
 GW_Bool FM_volume::IsIncluded( const FM_volumeVertex& vert )
{
	for (GW_U32 i=0; i<4; ++i)
	{
		if (vertex_[i]!=NULL)
			if (vertex_[i]==&vert)
				return GW_True;
	}
	return GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volume::SetID
/**
*  \param  nID [GW_U32] New ID
*  \author Birgit Stender
*  \date   22-02-2011
* 
*  Set the number of the volume in the mesh.
*/
/*------------------------------------------------------------------------------*/
void FM_volume::SetID(GW_U32 nID)
{
	nID_ = nID;
}


/*------------------------------------------------------------------------------*/
// Name : FM_volume::GetID
/**
*  \return [GW_U32] The ID
*  \author Birgit Stender
*  \date   22-02-2011
* 
*  Get the ID of the volume.
*/
/*------------------------------------------------------------------------------*/
GW_U32 FM_volume::GetID() const
{
	return nID_;
}

} // End namespace FM


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
