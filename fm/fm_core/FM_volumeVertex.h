
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_volumeVertex.h
 *  \brief  Definition of class \c FM_volumeVertex
 *  \author Birgit Stender
 *  \date   27-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_VOLUMEVERTEX_H_
#define _FM_VOLUMEVERTEX_H_

#include "../../gw/gw_core/GW_Config.h"
#include "../../gw/gw_core/GW_Vertex.h"
#include "FM_volumeVertexIterator.h"
#include "FM_volumeIterator.h"
#include "FM_volume.h"

namespace FM {

	using namespace GW;

	class FM_volumeVertexIterator;
	class FM_volumeIterator;
	class FM_volume;


	/*------------------------------------------------------------------------------*/
	/** 
	*  \class  FM_volumeVertex
	*  \brief  A vertex with a pointer onto the volume belonging to it
	*  \author Birgit Stender
	*  \date   27_02_2011
	*
	*/ 
	/*------------------------------------------------------------------------------*/

	class FM_volumeVertex: public GW_Vertex
	{
	public:

			/*------------------------------------------------------------------------------*/
			/** \name Constructor and destructor */
			/*------------------------------------------------------------------------------*/
			//@{
			FM_volumeVertex();
			virtual ~FM_volumeVertex();
			//@}

			void SetVolume( FM_volume& volume );
			FM_volume* GetVolume();
			const FM_volume* GetVolume() const;

			GW_U32 GetNumberNeighbor();
			//GW_Bool IsBoundaryVertex();

			//-------------------------------------------------------------------------
			/** \name Iterator management. */
			//-------------------------------------------------------------------------
			//@{
			FM_volumeIterator BeginVolumeIterator();
			FM_volumeIterator EndVolumeIterator();

			FM_volumeVertexIterator BeginVolumeVertexIterator();
			FM_volumeVertexIterator EndVolumeVertexIterator();
			//@}

	protected:


	private:
		/** A pointer on the volume which owned the vertex */
		FM_volume* pVolume_;
	};

/*------------------------------------------------------------------------------*/
/** \name a vector of FM_volumeVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class FM_volumeVertex*> T_volumeVertexVector;
typedef T_volumeVertexVector::iterator IT_volumeVertexVector;
typedef T_volumeVertexVector::reverse_iterator RIT_volumeVertexVector;
typedef T_volumeVertexVector::const_iterator CIT_volumeVertexVector;
typedef T_volumeVertexVector::const_reverse_iterator CRIT_volumeVertexVector;
//@}

/*------------------------------------------------------------------------------*/
/** \name a list of FM_volumeVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::list<class FM_volumeVertex*> T_volumeVertexList;
typedef T_volumeVertexList::iterator IT_volumeVertexList;
typedef T_volumeVertexList::reverse_iterator RIT_volumeVertexList;
typedef T_volumeVertexList::const_iterator CIT_volumeVertexList;
typedef T_volumeVertexList::const_reverse_iterator CRIT_volumeVertexList;
//@}

/*------------------------------------------------------------------------------*/
/** \name a map of FM_volumeVertex */
/*------------------------------------------------------------------------------*/
//@{
typedef std::map<GW_U32, class GW_volumeVertex*> T_volumeVertexMap;
typedef T_volumeVertexMap::iterator IT_volumeVertexMap;
typedef T_volumeVertexMap::reverse_iterator RIT_volumeVertexMap;
typedef T_volumeVertexMap::const_iterator CIT_volumeVertexMap;
typedef T_volumeVertexMap::const_reverse_iterator CRIT_volumeVertexMap;
//@}


} // End namespace FM


#endif // _FM_VOLUMEVERTEX_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////