/*------------------------------------------------------------------------------*/
/** 
*  \file   FM_volumeVertexIterator.cpp
*  \brief  Definition of class \c FM_volumeVertexIterator
*  \author Birgit Stender
*  \date   27-02-2011
*/ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_volumeVertexIterator.h"
#include "FM_volume.h"

using namespace FM;

FM_volumeVertexIterator::FM_volumeVertexIterator(  FM_volumeVertex* pStartVert, FM_volumeVertex* pVert, 
												 FM_volumeVertex* pVertNeighbor, GW_Bool pOnVolume,	
												 GW_U32 nNbrIncrement ):
pStartVert_				( pStartVert ),
pVert_					( pVert ),
pVertNeighbor_			( pVertNeighbor ),
pOnVolume_				( pOnVolume ),
nNbrIncrement_			( nNbrIncrement )
{ 
	pVolume_ = NULL;
	if (pStartVert_!=NULL)
	{
		pStartVolume_ = pStartVert_->GetVolume();
		if (pStartVolume_!= NULL && pVert_ != NULL)
		{

		}
		else
		{
			pStartVolume_ = NULL;
			pVolume_ = NULL;
			pStartVert_ = NULL;
			pVert_ = NULL;
		}
	}
	else
	{
			pStartVolume_ = NULL;
			pVolume_ = NULL;
			pVert_ = NULL;
	}
}

/* assignment */
FM_volumeVertexIterator& FM_volumeVertexIterator::operator=( const FM_volumeVertexIterator& it)
{
	this->pStartVert_		= it.pStartVert_;
	this->pVert_			= it.pVert_;
	this->pStartVolume_		= it.pStartVolume_;
	this->pVolume_			= it.pVolume_;
	this->pOnVolume_		= it.pOnVolume_;
	this->pVertNeighbor_	= it.pVertNeighbor_;
	this->nNbrIncrement_	= it.nNbrIncrement_;
	return *this;
}

/* egality */
GW_Bool FM_volumeVertexIterator::operator==( const FM_volumeVertexIterator& it)
{
	return (this->pStartVert_==it.pStartVert_)
		&& (this->pVert_==it.pVert_)
		&& (this->pStartVolume_==it.pStartVolume_)
		&& (this->pVolume_==it.pVolume_)
		&& (this->pOnVolume_==it.pOnVolume_)
		&& (this->pVertNeighbor_==it.pVertNeighbor_);
}

/* egality */
GW_Bool FM_volumeVertexIterator::operator!=( const FM_volumeVertexIterator& it)
{
	return (this->pStartVert_!=it.pStartVert_)
		|| (this->pVert_!=it.pVert_)
		|| (this->pStartVolume_!=it.pStartVolume_)
		|| (this->pVolume_!=it.pVolume_)
		|| (this->pOnVolume_!=it.pOnVolume_)
		|| (this->pVertNeighbor_!=it.pVertNeighbor_);
}


/* egality */
FM_volumeVertex* FM_volumeVertexIterator::operator*(  )
{
	if (pOnVolume_)
		return this->pVert_;
	else
		return this->pVertNeighbor_;
}

void FM_volumeVertexIterator::operator++()
{
	if( this->nNbrIncrement_>3 )
	{
		GW_ASSERT( GW_False );
		(*this) = FM_volumeVertexIterator(NULL,NULL);
		return;
	}

	if (pStartVert_!=NULL && pStartVolume_!=NULL)
	{
		if (pOnVolume_)
		{
			pVolume_ = pStartVolume_->GetVolumeNeighbor(*pVert_);
			if (pVolume_!=NULL)
			{
				pOnVolume_ = GW_False;
				// Neighbor of pStartVert_ is the one vertex on the volume neighbors
				// corresponding to the three other vertexes of pStartVolume_
				// adjacent to pStartVolume_
				GW_I32 pStartVolInNeighborID = pVolume_->GetVertexNumber(*pStartVolume_);
				if (pStartVolInNeighborID>=0)
				{
					pVertNeighbor_ = pVolume_->GetVertex(pStartVolInNeighborID);
					(*this) = FM_volumeVertexIterator( pStartVert_, pVert_, pVertNeighbor_, pOnVolume_, nNbrIncrement_+1 );
					return;
				}
			}
		}
		// pOnVolume_ needs to be set to true to get back to a vertex on pStartVolume_ 
		pOnVolume_ = GW_True;
		pVertNeighbor_ = NULL;
		pVert_ = pStartVolume_->GetNextVertex(*pVert_);
		if (pVert_!=pStartVert_)
		{
			(*this) = FM_volumeVertexIterator( pStartVert_, pVert_, pVertNeighbor_, pOnVolume_, nNbrIncrement_+1 );
			return;
		}
		// we reached again the starting vertex
		else
		{
			(*this) = FM_volumeVertexIterator( NULL, NULL );
			return;
		}
	}
	else
	{
		(*this) = FM_volumeVertexIterator( NULL, NULL );
		return;	
	}
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
