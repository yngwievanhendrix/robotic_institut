/*------------------------------------------------------------------------------*/
/** 
*  \file   FM_volumeMesh.cpp
*  \brief  Definition of class \c FM_volumeMesh
*  \author Birgit Stender
*  \date   23-02-2011
*/ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_volumeVertexIterator.h"
#include "FM_volumeMesh.h"

using namespace FM;

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh constructor
/**
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Constructor.
*/
/*------------------------------------------------------------------------------*/
FM_volumeMesh::FM_volumeMesh()
{
	vertexVector_.resize(0);
	volumeVector_.resize(0);
	vertexNeighbors_ = NULL;
	volumeNeighbors_ = NULL;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh destructor
/**
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  destructor.
*/
/*------------------------------------------------------------------------------*/
FM_volumeMesh::~FM_volumeMesh()
{
	for( GW_U32 i=0; i<this->GetNbrVertex(); i++ )
		GW_SmartCounter::CheckAndDelete( this->GetVertex(i) );
	for( GW_U32 i=0; i<this->GetNbrVolume(); i++ )
		GW_SmartCounter::CheckAndDelete( this->GetVolume(i) );
	GW_DELETE(vertexNeighbors_);
	GW_DELETE(volumeNeighbors_);
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::operator=
/**
*  \param  v [FM_volumeMesh&] Mesh.
*  \return [FM_volumeMesh&] *this
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Copy operator.
*/
/*------------------------------------------------------------------------------*/
FM_volumeMesh& FM_volumeMesh::operator=(const FM_volumeMesh& Mesh)
{
	this->SetNbrVertex( Mesh.GetNbrVertex() );
	this->SetNbrVolume( Mesh.GetNbrVolume() );

	for( GW_U32 i=0; i<this->GetNbrVertex(); ++i )
	{
		if( this->GetVertex(i)==NULL )
		{
			FM_volumeVertex& NewVert = this->CreateNewVertex();
			this->SetVertex(i, &NewVert);
		}
		FM_volumeVertex& NewVert = *this->GetVertex(i);
		const FM_volumeVertex& OriginalVert = *Mesh.GetVertex(i);
		NewVert = OriginalVert;

		/* resolve face attachement */
		/*const GW_Face* pFace = OriginalVert.GetFace();
		if( pFace!=NULL )
		NewVert.SetFace( *this->GetFace( pFace->GetID() ) );*/

		/* resolve volume attachement */
		const FM_volume* pVolume = OriginalVert.GetVolume();
		if( pVolume!=NULL )
			NewVert.SetVolume( *this->GetVolume( pVolume->GetID() ) );
	}
	for( GW_U32 i=0; i<this->GetNbrVolume(); ++i )
	{
		if( this->GetVolume(i)==NULL )
		{
			FM_volume& NewVolume = this->CreateNewVolume();
			this->SetVolume(i, &NewVolume);
		}
		FM_volume& NewVolume = *this->GetVolume(i);
		const FM_volume& OriginalVolume = *Mesh.GetVolume(i);
		NewVolume = OriginalVolume;
		/* resolve Vertex and neighbor */
		for( GW_U32 j=0; j<4; ++j )
		{
			GW_U32 VertID = OriginalVolume.GetVertex(j)->GetID();
			NewVolume.SetVertex( *this->GetVertex(VertID), j );
			const FM_volume* pNeigh = OriginalVolume.GetVolumeNeighbor(j);
			if( pNeigh==NULL )
				NewVolume.SetVolumeNeighbor(NULL, j);
			else
				NewVolume.SetVolumeNeighbor( this->GetVolume( pNeigh->GetID() ), j);
		}
	}

	return *this;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::GetNbrVertex
/**
*  \return [GW_U32] The number of vertexes
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Get the number of vertexes *allocated* in the mesh.
*/
/*------------------------------------------------------------------------------*/
GW_U32 FM_volumeMesh::GetNbrVertex() const
{
	return (GW_U32) vertexVector_.size();
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::GetNbrVolumes
/**
*  \return [GW_U32] the number.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Get the number of volumes *allocated* in the mesh.
*/
/*------------------------------------------------------------------------------*/
GW_U32 FM_volumeMesh::GetNbrVolume() const
{
	return (GW_U32) volumeVector_.size();
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::SetFace
/**
*  \param  nNum [GW_U32] The number of the new face.
*  \param  pFace [FM_volumeFace*] The new face.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Set one of the face.
*/
/*------------------------------------------------------------------------------*/
/*
void FM_volumeMesh::SetFace( GW_U32 nNum, FM_volumeFace* pFace )
{
GW_ASSERT( nNum<this->GetNbrFace() );
if( this->GetFace(nNum)!=NULL )
GW_SmartCounter::CheckAndDelete( this->GetFace(nNum) );
faceVector_[nNum] = pFace;
if( pFace!=NULL )
{
pFace->UseIt();
pFace->SetID(nNum);
}
}
*/

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::AddFace
/**
*  \param  Face [FM_volumeFace&] The new face.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Add a new face at the end of the list.
*/
/*------------------------------------------------------------------------------*/
/*
void FM_volumeMesh::AddFace( FM_volumeFace& Face )
{
this->SetNbrFace( this->GetNbrFace()+1 );
this->SetFace( this->GetNbrFace()-1, &Face );
}
*/

/*------------------------------------------------------------------------------*/
// Name : GW_Mesh::GetFace
/**
*  \param  nNum [GW_U32] The number
*  \return [GW_Face*] The face.
*  \author Gabriel Peyré
*  \date   2-15-2003
* 
*  Get one of the faces.
*/
/*------------------------------------------------------------------------------*/
/*
FM_volumeFace* FM_volumeMesh::GetFace( GW_U32 nNum )
{
GW_ASSERT( nNum<this->GetNbrFace() );
return faceVector_[nNum];
}
*/

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::SetVolume
/**
*  \param  nNum [GW_U32] The number of the new volume.
*  \param  pVolume [FM_volume*] The new volume.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Set one of the volumes.
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::SetVolume( GW_U32 nNum, FM_volume* pVolume )
{
	GW_ASSERT( nNum<this->GetNbrVolume() );
	if( this->GetVolume(nNum)!=NULL )
		GW_SmartCounter::CheckAndDelete( this->GetVolume(nNum) );
	volumeVector_[nNum] = pVolume;
	if( pVolume!=NULL )
	{
		pVolume->UseIt();
		pVolume->SetID(nNum);
	}
}


/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::AddVolume
/**
*  \param  volume [FM_volume&] The new volume.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Add a new volume at the end of the list.
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::AddVolume( FM_volume& volume )
{
	this->SetNbrVolume( this->GetNbrVolume()+1 );
	this->SetVolume( this->GetNbrVolume()-1, &volume );
}


/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::SetVertex
/**
*  \param  nNum [GW_U32] The number of the vertex.
*  \param  pVert [FM_volumeVertex*] The vertex.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Set one of the vertexes.
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::SetVertex( GW_U32 nNum, FM_volumeVertex* pVert )
{
	GW_ASSERT( nNum<this->GetNbrVertex() );
	if( this->GetVertex(nNum)!=NULL )
		GW_SmartCounter::CheckAndDelete( this->GetVertex(nNum) );
	vertexVector_[nNum] = pVert;
	if( pVert!=NULL )
	{
		pVert->UseIt();
		pVert->SetID(nNum);
	}
}


/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::GetVolume
/**
*  \param  nNum [GW_U32] The number
*  \return [FM_volume*] The volume.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Get one of the volumes.
*/
/*------------------------------------------------------------------------------*/
FM_volume* FM_volumeMesh::GetVolume( GW_U32 nNum )
{
	GW_ASSERT( nNum<this->GetNbrVolume() );
	return volumeVector_[nNum];
}

const FM_volume* FM_volumeMesh::GetVolume( GW_U32 nNum )const 
{
	GW_ASSERT( nNum<this->GetNbrVolume() );
	return volumeVector_[nNum];
}


/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::GetVertex
/**
*  \param  nNum [GW_U32] The number.
*  \return [FM_volumeVertex*] The vertex.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Get one of the vertex.
*/
/*------------------------------------------------------------------------------*/
FM_volumeVertex* FM_volumeMesh::GetVertex( GW_U32 nNum )
{
	GW_ASSERT( nNum < (this->GetNbrVertex()) );
	return vertexVector_[nNum];
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::GetVertex
/**
*  \param  nNum [GW_U32] The number.
*  \return [FM_volumeVertex*] The vertex.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Get one of the vertex.
*/
/*------------------------------------------------------------------------------*/
const FM_volumeVertex* FM_volumeMesh::GetVertex( GW_U32 nNum ) const 
{
	if( nNum < (this->GetNbrVertex()) )
		return vertexVector_[nNum];
	else
		return NULL;
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::SetNbrVolume
/**
*  \param  nNum [GW_U32] New number of volumes.
*  \author Birgit Stender
*  \date   2-15-2011
* 
*  Resize the mesh.
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::SetNbrVolume( GW_U32 nNum )
{
	GW_U32 nOldSize = this->GetNbrVolume();
	if( nNum<nOldSize )
	{
		/* check if the vertex at the end should be deleted */
		for( GW_U32 i=nNum; i<nOldSize; ++i )
			GW_SmartCounter::CheckAndDelete( this->GetVolume( i ) );
		volumeVector_.resize( nNum );
	}
	if( nNum>nOldSize )
	{
		volumeVector_.resize( nNum );
		/* set to NULL newly appended pointers */
		for( GW_U32 i=nOldSize; i<nNum; ++i )
			this->SetVolume( i, NULL );
	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::SetNbrVertex
/**
*  \param  nNum [GW_U32] New number of vertex.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Resize the mesh.
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::SetNbrVertex( GW_U32 nNum )
{
	GW_U32 nOldSize = this->GetNbrVertex();
	if( nNum<nOldSize )
	{
		/* check if the vertex at the end should be deleted */
		for( GW_U32 i=nNum; i<nOldSize; ++i )
			GW_SmartCounter::CheckAndDelete( this->GetVertex( i ) );
		vertexVector_.resize( nNum );
	}
	if( nNum>nOldSize )
	{
		vertexVector_.resize( nNum );
		/* set to NULL newly appended pointers */
		for( GW_U32 i=nOldSize; i<nNum; ++i )
			this->SetVertex( i, NULL );
	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::BuildConnectivity
/**
*  \author 
*  \date   27-02-2011
* 
*  Call this method when you have set the vertex and the face.
*	This will set up the neighboorhood for each volume
*/
/*------------------------------------------------------------------------------*/
void FM_volumeMesh::BuildConnectivity()
{
#ifdef FM_DEBUG_OUTPUT
	FILE *debug_output;


	if ((debug_output = fopen("DebugFastMarching.dat", "w"))==NULL)
	{
		fprintf(stderr, "error opening debug log file.\n");
		return;
	}
	fprintf(stderr, "Aufbau der Connectivity:\n");
#endif

	GW_U32 NbrVertex = this->GetNbrVertex();
	T_volumeList* vertexToVolumeMap = new T_volumeList[NbrVertex];

	// building up the map vertex->volume
	for( IT_volumeVector it = volumeVector_.begin(); it!=volumeVector_.end(); ++it )
	{
		FM_volume* pVolume = *it;
		if (pVolume!= NULL)
		{
			for( GW_I32 i=0; i<4; ++i )
			{		
				FM_volumeVertex* pVert = pVolume->GetVertex(i);
				if (pVert!=NULL)
				{
					if(pVert->GetID() < NbrVertex)
					{
						GW_U32 vertID = pVert->GetID();
						// for every vertex pointers to the volumes including it
						// are stored in a list referred to by the vertex ID
						vertexToVolumeMap[vertID].push_back( pVolume );
					}
				}
			}
		}
	}


	/* now we can set up connectivity */
	for( IT_volumeVector it=volumeVector_.begin(); it!=volumeVector_.end(); ++it )
	{

		FM_volume* pVolume0 = *it;
		if (pVolume0!=NULL)
		{

#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "vol%u: ", pVolume0->GetID());
#endif

			/* set up the neighboring volumes of the 4 vertices */
			T_volumeList* pVolumeLists[4];		
			for( GW_U32 i=0; i<4; ++i )
			{
				FM_volumeVertex* pVert = pVolume0->GetVertex(i);
				if(pVert!=NULL)
				{
					pVolumeLists[i] = &vertexToVolumeMap[pVert->GetID()];
#ifdef FM_DEBUG_OUTPUT
					fprintf(debug_output, "ver%u: %u (-> vol%u) ", i,  pVert->GetID(), pVert->GetVolume()->GetID());
#endif
				}
			}
			/* compute neighbor in the 4 directions */
			for( GW_U32 i=0; i<4; ++i )
			{
#ifdef FM_DEBUG_OUTPUT
				fprintf(debug_output, "neighbor%u: ", i);
#endif
				FM_volume* pNeighbor = NULL;
				GW_U32 i1 = (i+1)%4;
				GW_U32 i2 = (i+2)%4;
				GW_U32 i3 = (i+3)%4;
				/* we must find the intersection of the surrounding volumes of these 3 vertexes*/
				GW_Bool bFind = GW_False;
				for( IT_volumeList it1 = pVolumeLists[i1]->begin(); it1!=pVolumeLists[i1]->end() && bFind!=GW_True; ++it1 )
				{
					FM_volume* pVolume1 = *it1;
					for( IT_volumeList it2 = pVolumeLists[i2]->begin(); it2!=pVolumeLists[i2]->end() && bFind!=GW_True; ++it2 )
					{
						FM_volume* pVolume2 = *it2;
						for ( IT_volumeList it3 = pVolumeLists[i3]->begin(); it3!=pVolumeLists[i3]->end() && bFind!=GW_True; ++it3 )
						{
							FM_volume* pVolume3 = *it3;
							// we are searching for the neighboring volume corresponding to index i
							// the neighboring volume with index i is the one common among the including volumes
							// of vertexes i1, i2 and i3 
							if( pVolume0!=pVolume1 && pVolume1==pVolume2 && pVolume2==pVolume3 )
							{
#ifdef FM_DEBUG_OUTPUT
								fprintf(debug_output, "%u ", pVolume1->GetID());
#endif
								pNeighbor = pVolume1;
								bFind=GW_True;
							}
						}
					}
				}
				pVolume0->SetVolumeNeighbor( pNeighbor, i );
				/* make some test on the neighbor to assure symmetry
				in the connectivity relationship */
				//if( pNeighbor!=NULL )
				//{
				//	// get the vertex of pNeighbor opposite to the vertexes of the common face between
				//	// pVolume and pNeighbor
				//	GW_I32 nVertexNumber = pNeighbor->GetVertexNumber( *pVolume0->GetVertex(i1),*pVolume0->GetVertex(i2), *pVolume0->GetVertex(i3) );
				//	if (nVertexNumber>=0 )
				//		pNeighbor->SetVolumeNeighbor( pVolume0 , nVertexNumber );
				//}
			}
#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "\n");
#endif
		}
	}

//	// *************************************************************************************
//	// Building up the vertex neighborhood
//	// *************************************************************************************
//
//	GW_DELETEARRAY( vertexNeighbors_ );
//	vertexNeighbors_ = new T_volumeVertexList[NbrVertex];
//
//#ifdef FM_DEBUG_OUTPUT
//	fprintf(debug_output, "\nVertex Neighborhood:");
//#endif
//
//	for( IT_volumeVertexVector it = vertexVector_.begin(); it!=vertexVector_.end(); ++it )
//	{
//		FM_volumeVertex* pVert = *it;
//		GW_U32 vertID = pVert->GetID();
//#ifdef FM_DEBUG_OUTPUT
//		fprintf(debug_output, "\nVertex %u: ", vertID);
//#endif
//
//		for( FM_volumeVertexIterator vertIt = pVert->BeginVolumeVertexIterator(); vertIt!=pVert->EndVolumeVertexIterator(); ++vertIt )
//
//		{
//			FM_volumeVertex* pNeighborVert = *vertIt;
//			vertexNeighbors_[vertID].push_back( pNeighborVert );	
//#ifdef FM_DEBUG_OUTPUT
//			fprintf(debug_output, "%u ", pNeighborVert->GetID());
//#endif
//		}
//	}
//#ifdef FM_DEBUG_OUTPUT
//	fprintf(debug_output, "\n");
//#endif

	// *************************************************************************************
	// Building up the volume neighborhood
	// *************************************************************************************
	if (volumeNeighbors_ != NULL)
	{
		GW_DELETEARRAY( volumeNeighbors_ );
		T_volumeList* volumeNeighbors_ = new T_volumeList[NbrVertex];
	}
	else
	{
		volumeNeighbors_ = new T_volumeList[NbrVertex];
	}
	volumeNeighbors_ = vertexToVolumeMap;


#ifdef FM_DEBUG_OUTPUT
	fprintf(debug_output, "\nVolume Neighborhood:");
#endif

	for( IT_volumeVertexVector it = vertexVector_.begin(); it!=vertexVector_.end(); ++it )
	{
		FM_volumeVertex* pVert = *it;
		GW_U32 vertID = pVert->GetID();
		T_volumeList* pVolumeNeighborList;
		pVolumeNeighborList = &volumeNeighbors_[vertID];

#ifdef FM_DEBUG_OUTPUT
		fprintf(debug_output, "\nVertex %u: ", vertID);
#endif

		for( IT_volumeList volIt = pVolumeNeighborList->begin(); volIt!=pVolumeNeighborList->end(); ++volIt )
		{
			FM_volume* pNeighborVolume = *volIt;
#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "%u ", pNeighborVolume->GetID());
#endif
		}
	}

#ifdef FM_DEBUG_OUTPUT
	fprintf(debug_output, "\n");


#endif


	// *************************************************************************************
	//// Building up the vertex neighborhood
	//// *************************************************************************************
	GW_DELETEARRAY( vertexNeighbors_ );
	vertexNeighbors_ = new T_volumeVertexList[NbrVertex];


#ifdef FM_DEBUG_OUTPUT
	fprintf(debug_output, "\nVertex Neighborhood:");
#endif

	for( IT_volumeVertexVector it = vertexVector_.begin(); it!=vertexVector_.end(); ++it )
	{
		FM_volumeVertex* pVert = *it;
		GW_U32 vertID = pVert->GetID();
		T_volumeList* pVolumeNeighborList;
		pVolumeNeighborList = &volumeNeighbors_[vertID];

		for( IT_volumeList volIt = pVolumeNeighborList->begin(); volIt!=pVolumeNeighborList->end(); ++volIt )
		{
			FM_volumeVertex* pNeighborVert = pVert;
			FM_volume* pNeighborVolume = *volIt;
			for(GW_U32 i = 0; i<2; ++i)
			{
				pNeighborVert = pNeighborVolume->GetNextVertex(*pNeighborVert);
				if (pNeighborVert!=NULL)
				{
					vertexNeighbors_[vertID].push_back( pNeighborVert );

				}
			}
		}
	}

	for( IT_volumeVertexVector it = vertexVector_.begin(); it!=vertexVector_.end(); ++it )
	{
		FM_volumeVertex* pVert = *it;
		GW_U32 vertID = pVert->GetID();
		vertexNeighbors_[vertID].sort();
		vertexNeighbors_[vertID].unique();

#ifdef FM_DEBUG_OUTPUT
		fprintf(debug_output, "\nVertex %u: ", vertID);
#endif

		T_volumeVertexList* pVertexNeighborList;
		pVertexNeighborList = &vertexNeighbors_[vertID];
		T_volumeList* pVolumeNeighborList;


		for( IT_volumeVertexList VertIt = pVertexNeighborList->begin(); VertIt!=pVertexNeighborList->end(); ++VertIt )
		{
			FM_volumeVertex* pNeighborVert = *VertIt;

#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "%u ", pNeighborVert->GetID());
#endif

		}

	}
#ifdef FM_DEBUG_OUTPUT
	fprintf(debug_output, "\n");
#endif


#ifdef FM_DEBUG_OUTPUT
	fclose(debug_output);
#endif


	// free(vertexToVolumeMap); 
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::IterateConnectedComponent_Vertex
/**
*  \param  pCallback [VertexIterate_Callback] The function.
*  \author Gabriel Peyré
*  \date   7-2-2003
* 
*  Iterate a callback on each vertex of a connected component.
*/
/*------------------------------------------------------------------------------*/

//void GW_Mesh::IterateConnectedComponent_volumeVertex( FM_volumeVertex& start_vert, volumeVertexIterate_Callback pCallback )
//{
//	/* march on the voronoi diagram */
//	T_volumeVertexList VertexToProceed;
//	VertexToProceed.push_back( &start_vert );
//	T_volumeVertexMap VertexDone;
//	VertexDone[ start_vert.GetID() ] = &start_vert;
//
//
//	while( !VertexToProceed.empty() )
//	{
//		FM_volumeVertex* pVert = VertexToProceed.front();
//		GW_ASSERT( pVert!=NULL );
//		VertexToProceed.pop_front();
//
//		/* cut the vertex */
//		pCallback( *pVert );
//
//		/* add neighbors */
//		for( FM_volumeVertexIterator it = pVert->BeginVolumeVertexIterator(); it!=pVert->EndVolumeVertexIterator(); ++it )
//		{
//			FM_volumeVertex* pNewVert = (FM_volumeVertex*) *it;
//			if( pNewVert==NULL )
//				break;
//			GW_ASSERT( pNewVert!=NULL );
//			if( VertexDone.find(pNewVert->GetID())==VertexDone.end() )
//			{				
//				VertexToProceed.push_back( pNewVert );
//				VertexDone[ pNewVert->GetID() ] = pNewVert;	// so that it won't be added anymore
//			}
//		}
//	}
//}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::IterateConnectedComponent_vertex
/**
*  \param  pCallback [volumeIterate_Callback] The function.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Iterate a callback on each vertex of a connected component.
*/
/*------------------------------------------------------------------------------*/
//void FM_volumeMesh::IterateConnectedComponent_volume( FM_volume& start_volume, volumeIterate_Callback pCallback )
//{
//	/* march on the voronoi diagram */
//	T_volumeList VolumeToProceed;
//	VolumeToProceed.push_back( &start_volume );
//	T_volumeMap VolumeDone;
//	VolumeDone[ start_volume.GetID() ] = &start_volume;
//
//
//	while( !VolumeToProceed.empty() )
//	{
//		FM_volume* pVolume = VolumeToProceed.front();
//		GW_ASSERT( pVolume!=NULL );
//		VolumeToProceed.pop_front();
//
//		/* cut the volume */
//		pCallback( *pVolume );
//
//		/* add neighbors */
//		for( GW_U32 i=0; i<4; ++i )
//		{
//			FM_volume* pNewVolume = pVolume->GetVolumeNeighbor(i);
//			if( pNewVolume!=NULL && VolumeDone.find(pNewVolume->GetID())==VolumeDone.end() )
//			{				
//				VolumeToProceed.push_back( pNewVolume );
//				VolumeDone[ pNewVolume->GetID() ] = pNewVolume;	// so that it won't be added anymore
//			}
//		}
//	}
//}

void FM_volumeMesh::CheckIntegrity()
{
	for( GW_U32 i=0; i<this->GetNbrVertex(); ++i ) 
	{
		FM_volumeVertex* pVert = this->GetVertex(i); GW_ASSERT( pVert!=NULL );
		FM_volume* pVolume = pVert->GetVolume();	GW_ASSERT( pVolume!=NULL );
		if( pVolume!=NULL && pVolume->GetVertex(0)!=pVert &&
			pVolume->GetVertex(1)!=pVert &&
			pVolume->GetVertex(2)!=pVert )
			GW_ASSERT( GW_False );
	}
	for( GW_U32 i=0; i<this->GetNbrVolume(); ++i )
	{
		FM_volume* pVolume = this->GetVolume(i);	GW_ASSERT( pVolume!=NULL );
		for( GW_U32 k=0; k<4; ++k )
		{
			GW_U32 k1 = (k+1)%4;
			GW_U32 k2 = (k+2)%4;
			GW_U32 k3 = (k+3)%4;
			FM_volume* pNeighVolume = pVolume->GetVolumeNeighbor(k);
			FM_volumeVertex* pV1 = pVolume->GetVertex(k1);	GW_ASSERT( pV1!=NULL );
			FM_volumeVertex* pV2 = pVolume->GetVertex(k2);	GW_ASSERT( pV2!=NULL );
			FM_volumeVertex* pV3 = pVolume->GetVertex(k3);	GW_ASSERT( pV3!=NULL );
			if( pNeighVolume!=NULL )
			{
				GW_ASSERT( pNeighVolume->GetVolumeNeighbor(*pV1, *pV2, *pV3)==pVolume );
				GW_ASSERT( pVolume->GetVolumeNeighbor(*pV1, *pV2, *pV3)==pNeighVolume);
			}
		}
	}
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::CreateNewVertex
/**
*  \return [FM_volumeVertex&] The newly created vertex.
*  \author Birgit Stender
*  \date   03-03-2011
* 
*  Allocate memory for a new vertex. You should overload this 
*  method 
*/
/*------------------------------------------------------------------------------*/
FM_volumeVertex& FM_volumeMesh::CreateNewVertex()
{
	return *(new FM_volumeVertex);
}

/*------------------------------------------------------------------------------*/
// Name : FM_volumeMesh::CreateNewVolume
/**
*  \return [FM_volume&] The newly created volume.
*  \author Bigit Stender
*  \date   03-03-2011
* 
*  Allocate memory for a new face.
*/
/*------------------------------------------------------------------------------*/
FM_volume& FM_volumeMesh::CreateNewVolume()
{
	return *(new FM_volume);
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
