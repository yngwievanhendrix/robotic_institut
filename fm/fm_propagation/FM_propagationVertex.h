
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_propagationVertex.h
 *  \brief  Definition of class \c FM_propagationVertex
 *  \author Birgit Stender
 *  \date   01-03-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_PROPAGATIONVERTEX_H_
#define _FM_PROPAGATIONVERTEX_H_

#include "../../gw/gw_core/GW_Config.h"
#include "../fm_core/FM_volumeVertex.h"

namespace FM {

using namespace GW;

/*------------------------------------------------------------------------------*/
/** 
 *  \class  FM_propagationVertex
 *  \brief  A vertex that contains information to perform fast marching computations.
 *  \author Birgit Stender
 *  \date   01-03-2011
 *
 *  Contains computed distance, and current front.
 */ 
/*------------------------------------------------------------------------------*/

class FM_propagationVertex: public FM_volumeVertex
{

public:

	enum T_propagationVertexState
	{
		kFar,
		kAlive,
		kDead
	};

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    FM_propagationVertex();
    virtual ~FM_propagationVertex();
    //@}

    //-------------------------------------------------------------------------
    /** \name Accessors. */
    //-------------------------------------------------------------------------
    //@{
	GW_Float GetDistance();
	void SetDistance( GW_Float rDistance );
	void SetState( T_propagationVertexState nState );
	T_propagationVertexState GetState();
	FM_propagationVertex* GetFront();
	void SetFront( FM_propagationVertex* pFront );
    //@}

	void ResetPropagationVertex();

	void SetBoundaryReached( GW_Bool bBoundaryReached = GW_True );
	GW_Bool GetBoundaryReached();

	static GW_Bool CompareVertex(FM_propagationVertex* pVert1, FM_propagationVertex* pVert2);

    //-------------------------------------------------------------------------
    /** \name Parametrization helpers. */
    //-------------------------------------------------------------------------
    //@{
	void SetStoppingVertex( GW_Bool bIsStoppingVertex );
	GW_Bool GetIsStoppingVertex();
    //@}

private:

	/** current distance */
	GW_Float rDistance_;
	/** state of the vertex : can be far/alive/dead */
	T_propagationVertexState nState_;
	/** The vertex from which the front this vertex is in started.
	    Can be \c NULL if this vertex hasn't be reached by a front. */
	FM_propagationVertex* pFront_;


    //-------------------------------------------------------------------------
    /** \name specific for parametrization computation. */
    //-------------------------------------------------------------------------
    //@{
	GW_Bool bIsStoppingVertex_;
	GW_Bool bBoundaryReached_;
    //@}
 
};


/*------------------------------------------------------------------------------*/
/** \name a vector of FM_propagationVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class FM_propagationVertex*> T_propagationVertexVector;
typedef T_propagationVertexVector::iterator IT_propagationVertexVector;
typedef T_propagationVertexVector::reverse_iterator RIT_propagationVertexVector;
typedef T_propagationVertexVector::const_iterator CIT_propagationVertexVector;
typedef T_propagationVertexVector::const_reverse_iterator CRIT_propagationVertexVector;
//@}


/*------------------------------------------------------------------------------*/
/** \name a list of FM_propagationVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::list<class FM_propagationVertex*> T_propagationVertexList;
typedef T_propagationVertexList::iterator IT_propagationVertexList;
typedef T_propagationVertexList::reverse_iterator RIT_propaationVertexList;
typedef T_propagationVertexList::const_iterator CIT_propagationVertexList;
typedef T_propagationVertexList::const_reverse_iterator CRIT_propagationVertexList;
//@}

/*------------------------------------------------------------------------------*/
/** \name a map of FM_propagationVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::map<GW_U32, class FM_propagationVertex*> T_propagationVertexMap;
typedef T_propagationVertexMap::iterator IT_propagationVertexMap;
typedef T_propagationVertexMap::reverse_iterator RIT_propagationVertexMap;
typedef T_propagationVertexMap::const_iterator CIT_propagationVertexMap;
typedef T_propagationVertexMap::const_reverse_iterator CRIT_propagationVertexMap;
//@}

} // End namespace FM


#endif // _FM_PROPAGATIONVERTEX_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
