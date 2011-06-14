/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_propagationVertex.cpp
 *  \brief  Definition of class \c FM_propagationVertex
 *  \author Birgit Stender
 *  \date   01-03-2011
 */ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_propagationVertex.h"

using namespace FM;

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex constructor
/**
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Constructor.
 */
/*------------------------------------------------------------------------------*/
 FM_propagationVertex::FM_propagationVertex()
:	FM_volumeVertex	(),
	rDistance_	( GW_INFINITE ),
	nState_		( kFar ),
	pFront_		( NULL ),
	bIsStoppingVertex_	( GW_False ),
	bBoundaryReached_	( GW_False )
{
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex destructor
/**
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Destructor.
 */
/*------------------------------------------------------------------------------*/
 FM_propagationVertex::~FM_propagationVertex()
{
	/* NOTHING */
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::GetDistance
/**
 *  \return [GW_Float] Current distance. If the vertex is not "Dead", then it is not reliable.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Get the current distance.
 */
/*------------------------------------------------------------------------------*/
 GW_Float FM_propagationVertex::GetDistance()
{
	return rDistance_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::SetDistance
/**
 *  \param  rDistance [GW_Float] Current distance.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Set the current distance.
 */
/*------------------------------------------------------------------------------*/
 void FM_propagationVertex::SetDistance( GW_Float rDistance )
{
	rDistance_ = rDistance;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::SetState
/**
 *  \param  nState [T_propagationVertexState] The new state.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Set the state of the vertex.
 */
/*------------------------------------------------------------------------------*/
 void FM_propagationVertex::SetState( T_propagationVertexState nState )
{
	nState_ = nState;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::GetState
/**
 *  \return [T_propagationVertexState] The current state.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Return the current state of the vertex.
 */
/*------------------------------------------------------------------------------*/
 FM_propagationVertex::T_propagationVertexState FM_propagationVertex::GetState()
{
	return nState_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::GetFront
/**
 *  \return [FM_propagationVertex*] \c NULL if the vertex hasn't already been reached by a front.
 *  \author Birigt Stender
 *  \date   01-03-2011
 * 
 *  Get the vertex from which the front was started.
 */
/*------------------------------------------------------------------------------*/
 FM_propagationVertex* FM_propagationVertex::GetFront()
{
	return pFront_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::SetFront
/**
 *  \param  pFront [FM_propagationVertex*] The vertex from which the front started.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Set the front to which this vertex belongs.
 */
/*------------------------------------------------------------------------------*/
 void FM_propagationVertex::SetFront( FM_propagationVertex* pFront )
{
	pFront_ = pFront;
}


/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::CompareVertex
/**
 *  \param  pVert1 [FM_propagationVertex*] 1st vertex.
 *  \param  pVert2 [FM_propagationVertex*] 2nd vertex.
 *  \return [GW_Bool] True if distance of the 1st mesh is larger than to the one 
 *  of the 2nd.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Compare the distance of the 2 vertex. Used by the heap sorter.
 */
/*------------------------------------------------------------------------------*/
 GW_Bool FM_propagationVertex::CompareVertex(FM_propagationVertex* pVert1, FM_propagationVertex* pVert2)
{
	return pVert1->GetDistance()>pVert2->GetDistance();
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::ResetPropagationVertex
/**
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Reset the datas for fast marching computations.
 */
/*------------------------------------------------------------------------------*/
 void FM_propagationVertex::ResetPropagationVertex()
{
	rDistance_	= GW_INFINITE;
	nState_		= kFar;
	pFront_		= NULL;
	bIsStoppingVertex_	= GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::SetBoundaryReached
/**
 *  \param  bBoundaryReached [GW_Bool] Was it reached ?
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Is the vertex a boundary one reached by a parametrizing front ?
 */
/*------------------------------------------------------------------------------*/
 void FM_propagationVertex::SetBoundaryReached( GW_Bool bBoundaryReached )
{
	bBoundaryReached_ = bBoundaryReached;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::GetBoundaryReached
/**
 *  \return [GW_Bool] Answer.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Reached on a boundary ?
 */
/*------------------------------------------------------------------------------*/
 GW_Bool FM_propagationVertex::GetBoundaryReached()
{
	return bBoundaryReached_;
}

/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::SetStoppingVertex
/**
 *  \param  bIsStoppingVertex [GW_Bool] Yes/No
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Set on/off stopping criterion.
 */
/*------------------------------------------------------------------------------*/
void FM_propagationVertex::SetStoppingVertex( GW_Bool bIsStoppingVertex )
{
	bIsStoppingVertex_ = bIsStoppingVertex;
}


/*------------------------------------------------------------------------------*/
// Name : FM_propagationVertex::GetIsStoppingVertex
/**
 *  \return [GW_Bool] Stopping state.
 *  \author Birgit Stender
 *  \date   01-03-2011
 * 
 *  Get the stopping vertex state.
 */
/*------------------------------------------------------------------------------*/
GW_Bool FM_propagationVertex::GetIsStoppingVertex()
{
	return bIsStoppingVertex_;
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
