/*------------------------------------------------------------------------------*/
/** 
*  \file   FM_volumeIterator.cpp
*  \brief  Definition of class \c FM_volumeIterator
*  \author Birgit Stender
*  \date   27-02-2011
*/ 
/*------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "FM_volumeIterator.h"
#include "FM_volume.h"

using namespace FM;

FM_volumeIterator::FM_volumeIterator(  FM_volumeVertex* pStartVert, FM_volumeVertex* pVert, FM_volume* pVolume, GW_U32 nNbrIncrement )
:	pStartVert_		( pStartVert ), 
pVert_			( pVert  ), 
pVolume_		(pVolume ),
nNbrIncrement_	( nNbrIncrement )
{//1
	if (pStartVert_!=NULL && pVert_!=NULL)
	{//2
		pStartVolume_ = pStartVert_->GetVolume();
		if (pStartVolume_!= NULL)
		{//3
			if (pVolume_== NULL)
				pVolume_ = pStartVolume_;
		}
		// invalid start volume
		else
		{//3
				(*this) = FM_volumeIterator(NULL,NULL);
				return;
		}//2
	}//1
	// invalid start vertex
	else
	{
		pVert_ = NULL;
		pStartVolume_ = NULL;
		pVolume_ = NULL;
	}

}//0



/* assignement */
FM_volumeIterator& FM_volumeIterator::operator=( const FM_volumeIterator& it)
{
	this->pStartVert_		= it.pStartVert_;
	this->pVert_			= it.pVert_;
	this->pStartVolume_		= it.pStartVolume_;
	this->pVolume_			= it.pVolume_;
	this->nNbrIncrement_	= it.nNbrIncrement_;
	return *this;
}

/* egality */
GW_Bool FM_volumeIterator::operator==( const FM_volumeIterator& it)
{
	return (this->pStartVert_==it.pStartVert_)
		&& (this->pVert_==it.pVert_)
		&& (this->pStartVolume_==it.pStartVolume_)
		&& (this->pVolume_==it.pVolume_);
}

/* egality */
GW_Bool FM_volumeIterator::operator!=( const FM_volumeIterator& it)
{
	return (this->pStartVert_!=it.pStartVert_)
		|| (this->pVert_!=it.pVert_)
		|| (this->pStartVolume_!=it.pStartVolume_)
		|| (this->pVolume_!=it.pVolume_);
}


/* dereference */
FM_volume* FM_volumeIterator::operator*(  )
{
	return this->pVolume_;
}

/* progression : \todo take in acount NULL pointer */
void FM_volumeIterator::operator++()
{
	if( this->nNbrIncrement_>3 )
	{
		GW_ASSERT( GW_False );
		(*this) = FM_volumeIterator(NULL,NULL);
		return;
	}
	if (pStartVert_!=NULL && pStartVolume_!=NULL)
	{
		if (pVert_!=NULL && pVolume_!=NULL)
		{
			// Is pVertex included in pStartVolume_ ?
			if(pStartVolume_->GetVertexNumber(*pVert_)>=0)
			{
				pVert_ = pStartVolume_->GetNextVertex(*pVert_);
				if (pVert_!=pStartVert_)
				{
					if (pVert_!=NULL)
					{
						pVolume_= pStartVolume_->GetVolumeNeighbor(*pVert_);
						// everything is fine
						if (pVolume_!=NULL)
						{
							(*this) = FM_volumeIterator( pStartVert_, pVert_, pVolume_, nNbrIncrement_+1 );
							return;
						}
						// we reached a neighboring volume of pStartVolume_ which has not been initialized
						// try other direction of pStartVolume_ and save change of direction
						else
						{//5
							// try to find a valid vertex referring to a valid volume neighbor of pStartVolume_
							// by marching around pStartVolume as long as the vertex of pStartVolume_ is valid
							while(pVert_!=NULL && pVolume_==NULL)
							{//6
								pVert_= pStartVolume_->GetNextVertex(*pVert_);
								if (pVert_ != NULL)
								{//7
									// we reached the starting vertex
									if (pVert_==pStartVert_)
									{//8
										(*this) = FM_volumeIterator(NULL,NULL);
										return;
									}//7
									pVolume_ = pStartVolume_->GetVolumeNeighbor(*pVert_);
								}//6
							}//5
							// vertex pVert_ is valid and pVolume_ is valid
							// -> success!!!
							if (pVert_!=NULL && pVolume_!=NULL)
							{//6
								(*this) = FM_volumeIterator( pStartVert_, pVert_, pVolume_, nNbrIncrement_+1);
								return;
							}//5
							// we are stuck in a non valid vertex
							// or
							// pVolume_ still remains non valid (how?)
							else
							{
								(*this) = FM_volumeIterator(NULL,NULL);
								return;
							}

						}//4

					}
					// we reached an invalid vertex of pStartVolume_
					// abort
					else
					{
						(*this) = FM_volumeIterator(NULL,NULL);
						return;
					}
				}
				// end: we reached the starting vertex
				else
				{
					(*this) = FM_volumeIterator(NULL,NULL);
					return;
				}
			}
			// pVert_ is not included in pStartVolume_
			else
			{
				(*this) = FM_volumeIterator(NULL,NULL);
				return;
			}
		}
		// pVert_ or pVolume_ do not exist
		else
		{
			(*this) = FM_volumeIterator(NULL,NULL);
			return;
		}
	}
	// pStartVert_ or pStartVolume_ do not exist
	else
	{
		(*this) = FM_volumeIterator(NULL,NULL);
		return;
	}
	(*this) = FM_volumeIterator(NULL,NULL);
	return;
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
