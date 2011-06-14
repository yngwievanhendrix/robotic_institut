
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volume.h
 *  \brief  Definition of class \c FM_volume
 *  \author Birgit Stender
 *  \date   
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_VOLUME_H_
#define _FM_VOLUME_H_

#include "../../gw/gw_core/GW_Config.h"
#include "../../gw/gw_core/GW_SmartCounter.h"
#include "FM_volumeVertex.h"

namespace FM {

	using namespace GW;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  FM_volume
 *  \brief  A tetrahedron with its four vertexes and its four neighbors.
 *  \author Birgit Stender
 *  \date   2-12-2003
 *
 *  We use smart counter to check wether the face should descallocate a
 *	vertex. 
 *
 *	For vertex processing, a face is responsible for the vertex that caries
 *	a pointer on it. In fact, each vertex has a pointer on a single face
 *	that is responsible for that vertex.
 *
 *	Vertex, face and edge are labeled in a consistant way.
 *	This means :
 *		- edge are labeled by the number of the opposite vertex.
 *		- neighbor faces are labeled by the edge number, i.e. by the opposite vertex number.
 *
 *	We refer to these number as "edge number".
 */ 
/*------------------------------------------------------------------------------*/

class FM_volume:	public GW_SmartCounter
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    FM_volume();
	virtual ~FM_volume();
	virtual FM_volume& operator=(const FM_volume& volume);
    //@}

    //-------------------------------------------------------------------------
    /** \name Accessors */
    //-------------------------------------------------------------------------
    //@{
	void SetVolumeNeighbor(FM_volume* pVolume, GW_U32 nFaceNum);
	void SetVolumeNeighbor(FM_volume* pVolume0, FM_volume* pVolume1, FM_volume* pVolume2, FM_volume* pVolume3);
	FM_volume* GetVolumeNeighbor( GW_U32 nEdgeNum );
	const FM_volume* GetVolumeNeighbor( GW_U32 nEdgeNum ) const;
	FM_volume* GetVolumeNeighbor( const FM_volumeVertex& vert);
	FM_volume* GetVolumeNeighbor( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 );

	void SetVertex(FM_volumeVertex& vert, GW_U32 nNum);
	void SetVertex(FM_volumeVertex& vert0, FM_volumeVertex& vert1, FM_volumeVertex& vert2, FM_volumeVertex& vert3);
	FM_volumeVertex* GetVertex( GW_U32 nNum );
	const FM_volumeVertex* GetVertex( GW_U32 nNum ) const;
	FM_volumeVertex* GetVertex( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 );
	FM_volumeVertex* GetNextVertex( const FM_volumeVertex& vert, GW_Bool direction = GW_True );

	GW_I32 GetVertexNumber( const FM_volume& volume );
	GW_I32 GetVertexNumber( const FM_volumeVertex& vert );
	GW_I32 GetVertexNumber( const FM_volumeVertex& vert0, const FM_volumeVertex& vert1, const FM_volumeVertex& vert2 );

	GW_Bool IsResponsibleFor( const FM_volumeVertex& vert );
	GW_Bool IsIncluded( const FM_volumeVertex& vert );

	void	SetID(GW_U32 nID);
	GW_U32	GetID() const;

    //@}

private:

	/** our defining vertex */
	FM_volumeVertex*		vertex_[4];
	/** the 4 volumes arround us */
	FM_volume*		volumeNeighbors_[4];
	/** The ID of the face, given by the Mesh. Should be in the range [0,...,NbrVertex] */
	GW_U32 nID_;

};

/*------------------------------------------------------------------------------*/
/** \name a list of FM_volume */
/*------------------------------------------------------------------------------*/
//@{
typedef std::list<class FM_volume*> T_volumeList;
typedef T_volumeList::iterator IT_volumeList;
typedef T_volumeList::reverse_iterator RIT_volumeList;
typedef T_volumeList::const_iterator CIT_volumeList;
typedef T_volumeList::const_reverse_iterator CRIT_volumeList;
//@}


/*------------------------------------------------------------------------------*/
/** \name a vector of FM_volume */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class FM_volume*> T_volumeVector;
typedef T_volumeVector::iterator IT_volumeVector;
typedef T_volumeVector::reverse_iterator RIT_volumeVector;
typedef T_volumeVector::const_iterator CIT_volumeVector;
typedef T_volumeVector::const_reverse_iterator CRIT_volumeVector;
//@}

/*------------------------------------------------------------------------------*/
/** \name a map of FM_volume */
/*------------------------------------------------------------------------------*/
//@{
typedef std::map<GW_U32, class FM_volume*> T_volumeMap;
typedef T_volumeMap::iterator IT_volumeMap;
typedef T_volumeMap::reverse_iterator RIT_volumeMap;
typedef T_volumeMap::const_iterator CIT_volumeMap;
typedef T_volumeMap::const_reverse_iterator CRIT_volumeMap;
//@}



} // End namespace FM

#endif // _FM_VOLUME_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c)  Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
