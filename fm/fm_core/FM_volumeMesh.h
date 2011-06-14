
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volumeMesh.h
 *  \brief  Definition of class \c FM_volumeMesh
 *  \author Birgit Stender
 *  \date   22-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_VOLUMEMESH_H_
#define _FM_VOLUMEMESH_H_

#include "..\..\gw\gw_core\GW_Config.h"
#include "FM_Config.h"
#include "FM_volume.h"
#include "FM_volumeVertex.h"

namespace FM {

	using namespace GW;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  fm_volumeMesh
 *  \brief  A tetrahydral unstructured mesh 
 *  \author Birgit Stender
 *  \date   23-02-2011
 *
 *  Set up the connectivity using the following rules
 *		- vertex has a pointer to one of the faces it belongs to.
 *		- face has a pointer to one of the volumes it belongs to.
 *		- the 4 neighboring volume elements needs to be set up.
 *
 *	This class use smart pointer management. So it deletes pointers not used anymore
 *
 *	Explicit management of the volume and face numbers is required.
 */ 
/*------------------------------------------------------------------------------*/

class FM_volumeMesh
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    FM_volumeMesh();
    virtual ~FM_volumeMesh();
	virtual FM_volumeMesh& operator=(const FM_volumeMesh& v);
    //@}

    //-------------------------------------------------------------------------
    /** \name Resize manager. */
    //-------------------------------------------------------------------------
    //@{
	void SetNbrVolume( GW_U32 nNum );
	void SetNbrVertex( GW_U32 nNum );

	GW_U32 GetNbrVolume() const;
	GW_U32 GetNbrVertex() const;

	void Reset();
    //@}

	//-------------------------------------------------------------------------
    /** \name Vertex/Face management */
    //-------------------------------------------------------------------------
    //@{

	void SetVertex( GW_U32 nNum, FM_volumeVertex* pVert );
	FM_volumeVertex* GetVertex( GW_U32 nNum );
	const FM_volumeVertex* GetVertex( GW_U32 nNum ) const;

	/*
	void AddFace( FM_volumeFace& pFace );
	void SetFace( GW_U32 nNum, FM_volumeFace* pFace );
	GW_Face* GetFace( GW_U32 nNum );
	const GW_Face* GetFace( GW_U32 nNum ) const;
	*/

	void AddVolume( FM_volume& pVolume );
	void SetVolume( GW_U32 nNum, FM_volume* pVolume );
	FM_volume* GetVolume( GW_U32 nNum );
	const FM_volume* GetVolume( GW_U32 nNum ) const;
    //@}

	void BuildConnectivity();

	//-------------------------------------------------------------------------
    /** \name Class factory methods. */
    //-------------------------------------------------------------------------
    //@{
	virtual FM_volume& CreateNewVolume();
	virtual FM_volumeVertex& CreateNewVertex();
    //@}

	//typedef void (*volumeVertexIterate_Callback)( FM_volumeVertex& vert );
	//static void IterateConnectedComponent_volumeVertex( FM_volumeVertex& start_vert, volumeVertexIterate_Callback pCallback );
	//typedef void (*volumeIterate_Callback)( FM_volume& volume );
	//static void IterateConnectedCompontent_volume( FM_volume& start_volume, volumeIterate_Callback pCallback );

	/* helper*/
	void CheckIntegrity();

protected:

	/** contains all volumes of the mesh */
	T_volumeVector volumeVector_;
	/** contains all faces of the mesh */
	/*T_volumeFaceVector faceVector_;*/
	/** contains all faces of the mesh */
	T_volumeVertexVector vertexVector_;
	/** contains all vertex neighbors for each vertex */
	T_volumeVertexList* vertexNeighbors_;
	/** contains all volume neighbors of each vertex */
	T_volumeList* volumeNeighbors_;
};

/*------------------------------------------------------------------------------*/
/** \name a vector of FM_volumeMesh */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class FM_volumeMesh*> T_volumeMeshVector;
typedef T_volumeMeshVector::iterator IT_MeshVector;
typedef T_volumeMeshVector::reverse_iterator RIT_MeshVector;
typedef T_volumeMeshVector::const_iterator CIT_MeshVector;
typedef T_volumeMeshVector::const_reverse_iterator CRIT_MeshVector;
//@}


} // End namespace FM

#endif // _FM_VOLUMEMESH_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
