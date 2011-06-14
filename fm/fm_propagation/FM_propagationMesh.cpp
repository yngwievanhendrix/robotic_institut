/*------------------------------------------------------------------------------*/
/** 
*  \file   FM_propagationMesh.cpp
*  \brief  Definition of class \c FM_propagationMesh
*  \author Birgit Stender
*  \date   23-02-2011
*/ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_propagationMesh.h"

using namespace FM;

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh constructor
/**
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Constructor.
*/
/*------------------------------------------------------------------------------*/
FM_propagationMesh::FM_propagationMesh()
:	FM_volumeMesh(),
WeightCallback_		( FM_propagationMesh::BasicWeightCallback ),
VertexInsersionCallback_	( NULL ),
ForceStopCallback_			( NULL ),
NewDeadVertexCallback_		( NULL ),
HeuristicToGoalCallbackFunction_	( NULL ),
bIsMarchingBegin_			( GW_False ),
bIsMarchingEnd_				( GW_False )
{
	/* NOTHING */
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh destructor
/**
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Destructor.
*/
/*------------------------------------------------------------------------------*/
FM_propagationMesh::~FM_propagationMesh()
{
	// nothing 
}


/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::CreateNewVertex
/**
*  \return [FM_volumeVertex&] The newly created vertex.
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Allocate memory for a new vertex. You should overload this 
*  method 
*/
/*------------------------------------------------------------------------------*/
FM_volumeVertex& FM_propagationMesh::CreateNewVertex()
{
	return *(new FM_propagationVertex);
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::CreateNewVolume
/**
*  \return [FM_Volume&] The newly created volume.
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Allocate memory for a new face.
*/
/*------------------------------------------------------------------------------*/
FM_volume& FM_propagationMesh::CreateNewVolume()
{
	return *(new FM_propagationVolume);
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::ResetPropagationMesh
/**
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Reset the flags of all vertex for a new fast marching 
*  computation.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::ResetPropagationMesh()
{

	for( IT_volumeVertexVector it=vertexVector_.begin(); it!=vertexVector_.end(); ++it )
	{
		FM_propagationVertex* pVert = (FM_propagationVertex*) *it;
		pVert->ResetPropagationVertex();
	}
	ActiveVertex_.clear();

}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::AddStartVertex
/**
*  \param  StartVert [FM_propagationVertex&] The new starting point.
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Add a new vertex as starting point.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::AddStartVertex( FM_propagationVertex& StartVert )
{
	StartVert.SetFront( &StartVert );
	StartVert.SetDistance(0);
	StartVert.SetState( FM_propagationVertex::kAlive );
	ActiveVertex_.push_back( &StartVert );
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::PerformFastMarching
/**
*  \param  StartVertex [FM_propagationVertex&] The starting point.
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Compute geodesic distance from a vertex to other one.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::PerformFastMarching( FM_propagationVertex* pStartVertex )
{

#ifdef FM_DEBUG_OUTPUT
	FILE *debug_output;
#endif

	GW_U32 counter = 1;

	this->SetUpFastMarching( pStartVertex );

	// first time : set up the heap
	std::make_heap( ActiveVertex_.begin(), ActiveVertex_.end(), FM_propagationVertex::CompareVertex );
	/* main loop */

#ifdef FM_DEBUG_OUTPUT
	if ((debug_output = fopen("DebugFastMarching.dat", "a"))==NULL)
	{
		fprintf(stderr, "error opening debug log file.\n");
		return;
	}
#endif

	/*GW_SYSTEMTIME st;
	GetSystemTime(&st);
	fprintf(debug_output, "%d:%d:%d:%d\n" ,st.wHour,st.wMinute,st.wSecond,st.wMilliseconds);*/

	do
	{ 

#ifdef FM_DEBUG_OUTPUT
		if ((debug_output = fopen("DebugFastMarching.dat", "a"))==NULL)
		{
			fprintf(stderr, "error opening debug log file.\n");
			return;
		}
		fprintf(debug_output, "\n%u. Durchlauf FM_propagationMesh::PerformFastMarchingOneStep():\n", counter);
#endif

		for( IT_volumeVertexVector it=vertexVector_.begin(); it!=vertexVector_.end(); ++it )
		{
			FM_propagationVertex* pVert = (FM_propagationVertex*) *it;

#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "Vertex %u: Zustand %d, Distance %5.2f\n", pVert->GetID(), pVert->GetState(), pVert->	
				GetDistance());
#endif
		}
		//for( IT_volumeVector it=volumeVector_.begin(); it!=volumeVector_.end(); ++it )
		//{
		//	FM_propagationVolume* pVolume = (FM_propagationVolume*) *it;
		//	fprintf(debug_output, "Volume %u:\n", pVolume->GetID());

		//	for( GW_U32 i=0; i<4; ++i )
		//	{
		//		FM_volumeVertex* pVertInVolume = pVolume->GetVertex(i);
		//		fprintf(debug_output, "ver%u: %u (-> vol%u) ", i,  pVertInVolume->GetID(), pVertInVolume->GetVolume()->GetID());
		//		FM_volume* pVolumeNeighbor = pVolume->GetVolumeNeighbor(*pVertInVolume);
		//		if (pVolumeNeighbor!=NULL)
		//			fprintf(debug_output, "neighbor%u: %u ", i, pVolumeNeighbor->GetID());
		//		else
		//			fprintf(debug_output, "neighbor%u: nicht vorhanden", i);
		//		fprintf(debug_output, "\n");
		//	}
		//}

#ifdef FM_DEBUG_OUTPUT
		fclose(debug_output);
#endif
		counter++;
	}
	while( !this->PerformFastMarchingOneStep() );

	/*fprintf(debug_output, "%d:%d:%d:%d\n" ,st.wHour,st.wMinute,st.wSecond,st.wMilliseconds);*/
	//GetSystemTime(&st);
#ifdef FM_DEBUG_OUTPUT
	fclose(debug_output);
#endif

}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::SetUpFastMarching
/**
*  \param  pStartVertex=NULL [FM_propagationVertex*] A start vertex to add.
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Initialisation before starting the Fast Marching Algorithm
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::SetUpFastMarching( FM_propagationVertex* pStartVertex )
{
	GW_ASSERT( WeightCallback_!=NULL );

	if( pStartVertex!=NULL )
		this->AddStartVertex( *pStartVertex );

	std::make_heap( ActiveVertex_.begin(), ActiveVertex_.end(), FM_propagationVertex::CompareVertex );

	bIsMarchingBegin_ = GW_True;
	bIsMarchingEnd_ = GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::PerformFastMarchingOneStep
/**
*  \return [GW_Bool] Is the marching process finished ?
*  \author Birgit Stender
*  \date   23-02-2011
* 
*  Just one update step of the marching algorithm.
*/
/*------------------------------------------------------------------------------*/
GW_Bool FM_propagationMesh::PerformFastMarchingOneStep()
{
#ifdef FM_DEBUG_OUTPUT
	FILE *debug_output;
	if ((debug_output = fopen("DebugFastMarching.dat", "a"))==NULL)
	{
		fprintf(stderr, "error opening debug log file.\n");
	}
#endif

	// Termination if the TRIAL list ActiveVertex_ is empty
	if( ActiveVertex_.empty() )
		return GW_True;

	// check the state variable set in setUpFastMarching()
	// make_heap for ActiveVertex_ is also called in setUpFastMarching()
	GW_ASSERT( bIsMarchingBegin_ );

	//pCurVert represents the vertex with the larges time value
	FM_propagationVertex* pCurVert = ActiveVertex_.front();
	GW_ASSERT( pCurVert!=NULL );

#ifdef FM_DEBUG_OUTPUT
	fprintf(debug_output, "\nFront Vertex: vertex %u\n\n", pCurVert->GetID());
#endif

	// TRIAL value with the minimum T value is added to KNOWN
	// -> removed from the heap and vector ActiveVertex_
	// -> state of the vertex is set to dead
	std::pop_heap( ActiveVertex_.begin(), ActiveVertex_.end(), FM_propagationVertex::CompareVertex );
	ActiveVertex_.pop_back();
	pCurVert->SetState( FM_propagationVertex::kDead );

	if( NewDeadVertexCallback_!=NULL )
		NewDeadVertexCallback_( *pCurVert );

	GW_U32 curID = pCurVert->GetID();
	T_volumeVertexList* pVertexNeighborList;
	pVertexNeighborList = &vertexNeighbors_[curID];
	T_volumeList* pVolumeNeighborList;


	for( IT_volumeVertexList VertIt = pVertexNeighborList->begin(); VertIt!=pVertexNeighborList->end(); ++VertIt )
	{
		FM_propagationVertex* pNewVert = (FM_propagationVertex*) *VertIt;
		GW_ASSERT( pNewVert!=NULL );
		GW_U32 newID = pNewVert->GetID();
		pVolumeNeighborList = &volumeNeighbors_[newID];
		/* compute it's new value */

		if( (pCurVert->GetIsStoppingVertex() && !pNewVert->GetIsStoppingVertex() && pNewVert->GetState()==FM_propagationVertex::kFar) || (pNewVert->GetState()==FM_propagationVertex::kDead) )
		{
			// stopping vertex is not allowed to set vertexes to alive that are not stopping.
		}
		else
		{
			GW_U32 newVertID = pNewVert->GetID();
			GW_Bool curVertInVol;

#ifdef FM_DEBUG_OUTPUT
			fprintf(debug_output, "\nAktualisierung Vertex %u:\n", newVertID);
#endif


			/* computing the new distance using volume neighborhood information */
			GW_Float rNewDistance = pNewVert->GetDistance();
			for( IT_volumeList VolumeIt=pVolumeNeighborList->begin(); VolumeIt!=pVolumeNeighborList->end(); ++VolumeIt )
			{


				FM_propagationVolume* pVolume= (FM_propagationVolume*) *VolumeIt;
				GW_ASSERT( pVolume!=NULL );
				curVertInVol = pVolume->IsIncluded( *pCurVert );
				if (curVertInVol)
				{



					FM_propagationVertex* pVert1 = (FM_propagationVertex*) pVolume->GetNextVertex( *pNewVert );
					GW_ASSERT( pVert1!=NULL );
					FM_propagationVertex* pVert2 = (FM_propagationVertex*) pVolume->GetNextVertex( *pVert1 );
					GW_ASSERT( pVert2!=NULL );
					FM_propagationVertex* pVert3 = (FM_propagationVertex*) pVolume->GetNextVertex( *pVert2 );
					GW_ASSERT( pVert3!=NULL );

#ifdef FM_DEBUG_OUTPUT
					fprintf(debug_output, "durch Volumen %u mit Vertex %u, %u, %u, %u\n", 
						pVolume->GetID(), pNewVert->GetID(), pVert1->GetID(), pVert2->GetID(), pVert3->GetID());
#endif

					// reordering of the neighboring points to assure T(pVert1)<T(pVert2) and T(pVert1)<T(pVert3)
					if( pVert1->GetDistance() > pVert2->GetDistance())
					{
						if ( pVert3->GetDistance() > pVert2->GetDistance())
						{
							FM_propagationVertex* pTempVert = pVert2;
							pVert2 = pVert3;
							pVert3 = pTempVert;
						}
						else
						{
							FM_propagationVertex* pTempVert = pVert1;
							pVert1 = pVert2;
							pVert2 = pTempVert;
						}
					}


					rNewDistance = GW_MIN( rNewDistance, this->ComputeVolumeDistance( *pNewVert, *pVert1, *pVert2, *pVert3, *pCurVert->GetFront() ) );
					if(rNewDistance==GW_INFINITE)
					{
						rNewDistance=GW_INFINITE;						
					}

#ifdef FM_DEBUG_OUTPUT
					fprintf(debug_output, "neue distanz: %5.2f\n", rNewDistance);

#endif
				}

			}
			// state update
			switch( pNewVert->GetState() ) {
			case FM_propagationVertex::kFar:
				/* ask to the callback if we should update this vertex and add it to the path */
				//if( VertexInsersionCallback_==NULL ||
				//	VertexInsersionCallback_( *pNewVert,rNewDistance ) )
				//{
					pNewVert->SetDistance( rNewDistance );
					/* add the vertex to the heap */
					ActiveVertex_.push_back( pNewVert );
					std::push_heap( ActiveVertex_.begin(), ActiveVertex_.end(), FM_propagationVertex::CompareVertex );
					/* this one can be added to the heap */
					pNewVert->SetState( FM_propagationVertex::kAlive );
					pNewVert->SetFront( pCurVert->GetFront() );
				//}
				break;
			case FM_propagationVertex::kAlive:
				/* just update it's value */
				if( rNewDistance<=pNewVert->GetDistance() )
				{
					pNewVert->SetDistance( rNewDistance );
					pNewVert->SetFront( pCurVert->GetFront() );
					// hum, check if we can correct this (avoid recomputing the whole heap).
					std::make_heap( ActiveVertex_.begin(), ActiveVertex_.end(), FM_propagationVertex::CompareVertex );
				}
				else
				{
				}
				break;
			case FM_propagationVertex::kDead:
				break;
			default:
				GW_ASSERT( GW_False );
			}
		}
	}

	/* have we finished ? */
	bIsMarchingEnd_ = ActiveVertex_.empty();
	/* the user can force ending of the algorithm */
	if( ForceStopCallback_!=NULL && bIsMarchingEnd_==GW_False )
		bIsMarchingEnd_ = ForceStopCallback_(*pCurVert);

#ifdef FM_DEBUG_OUTPUT
	fclose(debug_output);
#endif
	return bIsMarchingEnd_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::PerformFastMarchingFlush
/**
*  \author Birgit Stender
*  \date   24-02-2011
* 
*  Setup and running the algorithm until it terminates
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::PerformFastMarchingFlush()
{
	if( !bIsMarchingBegin_ )
		this->SetUpFastMarching();

	/* main loop */
	while( !this->PerformFastMarchingOneStep() )
	{ }
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::IsFastMarchingFinished
/**
*  \return [GW_Bool] Response.
*  \author Birgit Stender
*  \date  24-02-2011
* 
*  Is the algorhm finished ?
*/
/*------------------------------------------------------------------------------*/
GW_Bool FM_propagationMesh::IsFastMarchingFinished()
{
	return bIsMarchingEnd_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::ComputeVolumeDistance
/**
*  \param  CurrentVertexD [FM_propagationVertex&] The vertex D to update.
*  \param  VertA [FM_propagationVertex&] Vertex A of the tetrahedral face opposite to D
*  \param  VertB [FM_propagationVertex&] Vertex B in the tetrahedral face opposite to D
*  \param  VertC [FM_propagationVertex&] Vertex C in the tetrahedral face opposite to D
*  \return The value of the distance calculated using Fermat's principle within the tetrahedron
*  \author Birgit Stender
*  \date   26-02-2011
* 
*  Compute the update of a vertex from inside of a tetrahedron.
*/
/*------------------------------------------------------------------------------*/
GW_Float FM_propagationMesh::ComputeVolumeDistance(FM_propagationVertex& CurrentVertexD, 
												   FM_propagationVertex& VertA, FM_propagationVertex& VertB, FM_propagationVertex& VertC, 
												   FM_propagationVertex& CurrentFront )
{	
	GW_Float F = this->WeightCallback_( CurrentVertexD );

	// assure that one of the vertexes in face opposite to D does not belong to FAR
	if( VertA.GetState()!=FM_propagationVertex::kFar ||
		VertB.GetState()!=FM_propagationVertex::kFar ||
		VertC.GetState()!=FM_propagationVertex::kFar)
	{

		// take only those vertexes into account with the same predecessor as the predecessor of D why?
		GW_Bool bVertAUsable = VertA.GetState()!=FM_propagationVertex::kFar; //&& VertA.GetFront()==&CurrentFront;
		GW_Bool bVertBUsable = VertB.GetState()!=FM_propagationVertex::kFar; //&& VertB.GetFront()==&CurrentFront;
		GW_Bool bVertCUsable = VertC.GetState()!=FM_propagationVertex::kFar; //&& VertC.GetFront()==&CurrentFront;

		GW_Vector3D EdgeAB = VertB.GetPosition() - VertA.GetPosition();
		GW_Float ab = EdgeAB.Norm();
		EdgeAB /= ab;
		GW_Vector3D EdgeAC = VertC.GetPosition() - VertA.GetPosition();
		GW_Float ac = EdgeAC.Norm();
		EdgeAC /= ac;
		GW_Vector3D EdgeAD = CurrentVertexD.GetPosition() - VertA.GetPosition();
		GW_Float ad = EdgeAD.Norm();
		EdgeAD /= ad;
		GW_Vector3D EdgeBC = VertC.GetPosition() - VertB.GetPosition();
		GW_Float bc = EdgeBC.Norm();
		EdgeBC /= bc;
		GW_Vector3D EdgeBD = CurrentVertexD.GetPosition() - VertB.GetPosition();
		GW_Float bd = EdgeBD.Norm();
		EdgeBD /= bd;
		GW_Vector3D EdgeCD = CurrentVertexD.GetPosition() - VertC.GetPosition();
		GW_Float cd = EdgeCD.Norm();
		EdgeCD /= cd;

		GW_Float dA = VertA.GetDistance();
		GW_Float dB	= VertB.GetDistance();
		GW_Float dC	= VertC.GetDistance();

		if ( !bVertAUsable ) 
		{
			if ( !bVertBUsable )
			{
				if ( !bVertCUsable )
					return GW_INFINITE;
				else
				{
					// update by traveltime along the edge CD
					return dC + cd * F;
				}
			}
			else
			{
				if ( !bVertCUsable )
					// update by traveltime along the edge BD
					return dB + bd * F;
				else
					// first-arrival local update by head wave travelling along face BCD					
					return FM_propagationMesh::ComputeUpdateFace( dB, dC, EdgeBC, EdgeBD, ab, bd, cd, F);
			}
		}
		else
		{
			if ( !bVertBUsable )
			{
				if ( !bVertCUsable )
				{
					// update by traveltime along the edge AD
					return dA + ad * F;
				}
				else
					// first-arrival local update by head wave travelling along face ACD
					return FM_propagationMesh::ComputeUpdateFace( dA, dC, EdgeAC, EdgeAD, ac, ad, cd, F);
			}
			else
				if ( !bVertCUsable )
					// first-arrival local update by head wave travelling along face ABD
					return FM_propagationMesh::ComputeUpdateFace( dA, dB, EdgeAB, EdgeAD, ab, ad, bd, F);
				else
				{
					{
						// first-arrival local update from within the tetrahedron
						// Open Topic: Unfolding of tetrahedrons with obtuse angles
						return FM_propagationMesh::ComputeUpdateTetrahedron( dA, dB, dC, EdgeAB, 
							EdgeAC, EdgeBC, EdgeAD, EdgeCD, ab, ac, bc, ad, bd, cd, F);

					}
				}
		}
	}
	return GW_INFINITE;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::ComputeUpdateTetrahedron
/**
*  \param  dA [GW_Float] Distance value of vertex A.
*  \param  dB [GW_Float] Distance value of vertex B.
*  \param  dC [GW_Float] Distance value of vertex C.
*  \param  EdgeAB [GW_Vector3D] Edge directed from vertex A to vertex B
*  \param  EdgeAC [GW_Vector3D] Edge directed from vertex A to vertex C
*  \param  EdgeAD [GW_Vector3D] Edge directed from vertex A to vertex D
*  \param  ab [GW_Float] Norm of EdgeAB.
*  \param  ac [GW_Float] Norm of EdgeAC.
*  \param  bc [GW_Float] Norm of EdgeBC.
*  \param  ad [GW_Float] Norm of EdgeAD.
*  \param  F [GW_Float] weight of the vertex D
*  \return [GW_Float] The update value of vertex D
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Compute the update using fermat's principle as described in the following paper:
*  "Computing first-arrival seismic traveltimes on unstructured 3-D tetrahedral grids
*  using the Fast Marching Method"
*  Peter G. Lelièvre, Colin G. Farquharson and Charles A. Hurchich
*  Geophys. J. Int. (2011) 184, 885-896
*  doi: 10.1111/j.13655-246X.2010.04880.x
*/
/*------------------------------------------------------------------------------*/
GW_Float FM_propagationMesh::ComputeUpdateTetrahedron( GW_Float dA, GW_Float dB, 
													  GW_Float dC, GW_Vector3D EdgeAB, 
													  GW_Vector3D EdgeAC, GW_Vector3D EdgeBC,
													  GW_Vector3D EdgeAD, GW_Vector3D EdgeCD,
													  GW_Float ab, GW_Float ac, GW_Float bc, 
													  GW_Float ad, GW_Float bd, GW_Float cd,
													  GW_Float F )
{
	GW_Float dD = GW_INFINITE;
	// Projection of EdgeAD onto EdgeAB
	GW_Float xi_0 = EdgeAD*EdgeAB;
	// Projection of EdgeAD onto EdgeAC
	GW_Float zeta_0 = EdgeAD*EdgeAC;
	GW_Float u = dB - dA;
	GW_Float v = dC - dA;

	// rho0: distance from vertex D to the point (xi_0, zeta_0);
	// Part of EdgeAD orthogonal to both EdgeAB and EdgeAC
	GW_Vector3D EdgeXi0Zeta0D  = EdgeAD*ad - EdgeAB*ab*(EdgeAD*EdgeAB) - EdgeAC*ac*(EdgeAD*EdgeAC);
	GW_Float rho_0 = EdgeXi0Zeta0D.Norm();
	// cos_alpha: EdgeAB*EdgeAC (alpha: enclosed angle)
	// d_power2: |EdgeAB|*|EdgeAC|*cos_alpha
	GW_Float d_power2 = ab*ac*(EdgeAB*EdgeAC);
	// sin_alpha = EdgeAB
	// phi: |EdgeAB|*|EdgeAC|*sin_alpha
	GW_Float phi =  ((EdgeAB-EdgeAB*(EdgeAB*EdgeAC))*EdgeAC)*ab*ac;
	GW_Float omega_power2 = F*F*phi*phi-u*u*ac*ac-v*v*ab*ab+2*u*v*d_power2;

	if (omega_power2 > 0 && phi > 0)
	{
		GW_Float omega = sqrt(omega_power2);
		GW_Float beta = u*ac*ac-v*d_power2;
		GW_Float gamma = v*ab*ab-u*d_power2;
		GW_Float xi = - GW_ABS(beta)*1/phi*1/omega*rho_0 + xi_0;
		GW_Float zeta = - GW_ABS(gamma)*1/phi*1/omega*rho_0 + zeta_0;

		if ( 0 < xi && xi < 1 && 0 < zeta && zeta < 1 && 0 < xi+zeta && xi+zeta < 1)
		{
			dD = dA + u*xi_0 + v*zeta_0 + omega*1/phi*rho_0;
		}
		//else
		//{
		//	dD = GW_INFINITE;
		//}


	}
	//else
	//{
	//	dD = -GW_INFINITE;
	//}

	// head wave traveltime face ABD
	GW_Float dD_headABD = FM_propagationMesh::ComputeUpdateFace( dA, dB, EdgeAB, EdgeAD, ab, ad, bd, F );
	// head wave traveltime face BCD
	GW_Float dD_headBCD = FM_propagationMesh::ComputeUpdateFace( dB, dC, EdgeBC, EdgeCD, bc, bd, cd, F );
	// head wave traveltime face ACD
	GW_Float dD_headACD = FM_propagationMesh::ComputeUpdateFace( dA, dC, EdgeAC, EdgeCD, ac, ad, cd, F );
	dD = GW_MIN(dD,dD_headABD);
	dD = GW_MIN(dD,dD_headBCD);
	dD = GW_MIN(dD,dD_headACD);

	return dD;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::ComputeUpdateFace
/**
*  \param  dA [GW_Float] Distance value of vertex A.
*  \param  dB [GW_Float] Distance value of vertex B.
*  \param  EdgeAB [GW_Vector3D] Edge directed from vertex A to vertex B
*  \param  EdgeAC [GW_Vector3D] Edge directed from vertex A to vertex C
*  \param  ab [GW_Float] Norm of EdgeAB.
*  \param  ac [GW_Float] Norm of EdgeAC.
*  \param  bc [GW_Float] Norm of EdgeBC.
*  \param  F [GW_Float] weight of the vertex C
*  \return [GW_Float] The update value of vertex C.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Compute the update using fermat's principle as described in the following paper:
*  "Computing first-arrival seismic traveltimes on unstructured 3-D tetrahedral grids
*  using the Fast Marching Method"
*  Peter G. Lelièvre, Colin G. Farquharson and Charles A. Hurchich
*  Geophys. J. Int. (2011) 184, 885-896
*  doi: 10.1111/j.13655-246X.2010.04880.x
*/
/*------------------------------------------------------------------------------*/
GW_Float FM_propagationMesh::ComputeUpdateFace(  GW_Float dA, GW_Float dB, GW_Vector3D EdgeAB, GW_Vector3D EdgeAC, GW_Float ab, GW_Float ac, GW_Float bc, GW_Float F )
{
	GW_Float dC = GW_INFINITE;
	// Projection of EdgeAC onto EdgeAB
	GW_Float xi_0 = EdgeAC*EdgeAB;
	// order of points assures u >= 0
	GW_Float u = dB - dA;

	// rho0: distance from vertex C to the point xi_0
	// Part of EdgeAC orthogonal to EdgeAB
	GW_Vector3D EdgeXi0C  = EdgeAC*ac - EdgeAB*ab*(EdgeAC*EdgeAB);
	GW_Float rho_0 = EdgeXi0C.Norm();

	GW_Float omega_power2 = F*F*ab*ab-u*u;

	if (omega_power2 > 0 && ab > 0)
	{
		GW_Float omega = sqrt(omega_power2);
		GW_Float xi = xi_0 - 1/omega*1/ab*u*rho_0;

		if ( 0 < xi && xi < 1)
		{
			dC = dA + u*xi_0 + omega*1/ab*rho_0;
		}
		//else
		//{
		//	dC = GW_INFINITE;
		//}


	}

	dC = GW_MIN(dC,dA+F*ac);
	dC = GW_MIN(dC,dB+F*bc);	

	return dC;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::BasicWeightCallback
/**
*  \param  Vert [FM_propagationVertex&] Current vertex.
*  \return [GW_Float] 1
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Just the constant function = 1.
*/
/*------------------------------------------------------------------------------*/
GW_Float FM_propagationMesh::BasicWeightCallback(FM_propagationVertex& Vert)
{
	return 1;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::RegisterForceStopCallbackFunction
/**
*  \param  pFunc [T_FastMarchingCallbackFunction] The function.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Set the function used to test if we should end the fast marching or not.
*	The function return GW_False if the algorithm should be stopped.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::RegisterForceStopCallbackFunction( T_FM_FastMarchingCallbackFunction pFunc )
{
	ForceStopCallback_ = pFunc;
}


/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::RegisterWeightCallbackFunction
/**
*  \param  pFunc [T_WeightCallbackFunction] The function.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Set the function used to define the metric on the mesh.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::RegisterWeightCallbackFunction( T_FM_WeightCallbackFunction pFunc )
{
	GW_ASSERT( pFunc!=NULL );
	WeightCallback_ = pFunc;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::RegisterVertexInsersionCallbackFunction
/**
*  \param  pFunc [T_VertexInsersionCallbackFunction] New function.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Set the function we use when trying to insert a new vertex.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::RegisterVertexInsersionCallbackFunction( T_FM_VertexInsersionCallbackFunction pFunc )
{
	VertexInsersionCallback_ = pFunc;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::RegisterNewDeadVertexCallbackFunction
/**
*  \param  pFunc [T_NewDeadVertexCallbackFunction] New function.
*  \author Birgit Stender
*  \date   27-02-2011
* 
*  Set the function we use when a new dead vertex is set.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::RegisterNewDeadVertexCallbackFunction( T_FM_NewDeadVertexCallbackFunction pFunc )
{
	NewDeadVertexCallback_ = pFunc;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationMesh::RegisterHeuristicToGoalCallbackFunction
/**
*  \param  pFunc [T_HeuristicToGoalCallbackFunction] Callback function.
*  \author Birgit Stender
*  \date   01-03-2011
* 
*  Turn the propagation into an A* like.
*/
/*------------------------------------------------------------------------------*/
void FM_propagationMesh::RegisterHeuristicToGoalCallbackFunction( T_FM_HeuristicToGoalCallbackFunction pFunc )
{
	HeuristicToGoalCallbackFunction_ = pFunc;
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
