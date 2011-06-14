
/*------------------------------------------------------------------------------*/
/** 
 *  \file   FM_propagationMesh.h
 *  \brief  Definition of class \c FM_propagationMesh
 *  \author Birgit Stender
 *  \date   23-02-2011
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _FM_PROPAGATIONMESH_H_
#define _FM_PROPAGATIONMESH_H_

#include "../../gw/gw_core/GW_Config.h"
#include "../fm_core/FM_Config.h"
#include "../fm_core/FM_volumeMesh.h"
#include "../fm_core/FM_volume.h"
#include "../fm_core/FM_volumeVertex.h"
#include "../fm_core/FM_volumeIterator.h"
#include "../fm_core/FM_volumeVertexIterator.h"
#include "FM_propagationVertex.h"
#include "FM_propagationVolume.h"


namespace FM {

	using namespace GW;


/*------------------------------------------------------------------------------*/
/** 
 *  \class  FM_propagationMesh
 *  \brief  A volume mesh for propagation of interfaces using the fast marching algorithm in volume
 *  \author Birgit Stender
 *  \date   23_02_2011
 *
 */ 
/*------------------------------------------------------------------------------*/

class FM_propagationMesh: public FM_volumeMesh
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    FM_propagationMesh();
    virtual ~FM_propagationMesh();
    //@}

	//-------------------------------------------------------------------------
	/** \name Class factory methods. */
	//-------------------------------------------------------------------------
	//@{
	virtual FM_volumeVertex& CreateNewVertex();
	virtual FM_volume& CreateNewVolume();
	//@}

    //-------------------------------------------------------------------------
    /** \name Fast marching computations. */
    //-------------------------------------------------------------------------
	//@{
	void ResetPropagationMesh();
	void AddStartVertex( FM_propagationVertex& StartVert );
	void PerformFastMarching( FM_propagationVertex* pStartVertex=NULL );
	void SetUpFastMarching( FM_propagationVertex* pStartVertex=NULL );
	GW_Bool PerformFastMarchingOneStep();
	void PerformFastMarchingFlush();
	GW_Bool IsFastMarchingFinished();
    //@}

    //-------------------------------------------------------------------------
    /** \name Callback management. */
    //-------------------------------------------------------------------------
    //@{/*
        typedef GW_Float (*T_FM_WeightCallbackFunction)( FM_propagationVertex& Vert );
        void RegisterWeightCallbackFunction( T_FM_WeightCallbackFunction pFunc );
        typedef GW_Bool (*T_FM_FastMarchingCallbackFunction)( FM_propagationVertex& Vert );
        void RegisterForceStopCallbackFunction( T_FM_FastMarchingCallbackFunction pFunc );
        typedef void (*T_FM_NewDeadVertexCallbackFunction)( FM_propagationVertex& Vert );
        void RegisterNewDeadVertexCallbackFunction( T_FM_NewDeadVertexCallbackFunction pFunc );
        typedef GW_Bool (*T_FM_VertexInsersionCallbackFunction)( FM_propagationVertex& Vert, GW_Float rNewDist );
        void RegisterVertexInsersionCallbackFunction( T_FM_VertexInsersionCallbackFunction pFunc );
        typedef GW_Float (*T_FM_HeuristicToGoalCallbackFunction)( FM_propagationVertex& Vert );
        void RegisterHeuristicToGoalCallbackFunction( T_FM_HeuristicToGoalCallbackFunction pFunc );
	//@}

//        typedef GW_Float (*T_FM_WeightCallbackFunction)( FM_propagationVertex& );
//        void RegisterWeightCallbackFunction( T_FM_WeightCallbackFunction pFunc );
//        typedef GW_Bool (*T_FM_FastMarchingCallbackFunction)( FM_propagationVertex& );
//        void RegisterForceStopCallbackFunction( T_FM_FastMarchingCallbackFunction pFunc );
//        typedef void (*T_FM_NewDeadVertexCallbackFunction)( FM_propagationVertex& );
//        void RegisterNewDeadVertexCallbackFunction( T_FM_NewDeadVertexCallbackFunction pFunc );
//        typedef GW_Bool (*T_FM_VertexInsersionCallbackFunction)( FM_propagationVertex&, GW_Float  );
//        void RegisterVertexInsersionCallbackFunction( T_FM_VertexInsersionCallbackFunction pFunc );
//        typedef GW_Float (*T_FM_HeuristicToGoalCallbackFunction)( FM_propagationVertex&);
//        void RegisterHeuristicToGoalCallbackFunction( T_FM_HeuristicToGoalCallbackFunction pFunc );

	static GW_Float BasicWeightCallback(FM_propagationVertex& Vert);

protected:

	/** should be filled with the starting point of the marching before
	    calling PerformFastMarching */
	T_propagationVertexVector ActiveVertex_;

	/** a function that specify the metric on the mesh */
	T_FM_WeightCallbackFunction WeightCallback_;
	/** the callback function used to test if we should terminate the fast marching or not */
	T_FM_FastMarchingCallbackFunction ForceStopCallback_;
	/** the callback function used to test when a new dead vertex is created */
	T_FM_NewDeadVertexCallbackFunction NewDeadVertexCallback_;
	/** a function called to know if we should insert this vertex */
	T_FM_VertexInsersionCallbackFunction VertexInsersionCallback_;
	/** a function called to give an heuristic for the remaining distance */
	T_FM_HeuristicToGoalCallbackFunction HeuristicToGoalCallbackFunction_;
	
	/** just to controle interactive mode */
	GW_Bool bIsMarchingBegin_;
	GW_Bool bIsMarchingEnd_;


private:

	GW_Float ComputeVolumeDistance( FM_propagationVertex& CurrentVertexD, 
									FM_propagationVertex& VertA, FM_propagationVertex& VertB, FM_propagationVertex& VertC,
									FM_propagationVertex& CurrentFront );

	static GW_Float ComputeUpdateTetrahedron( GW_Float dA, GW_Float dB, 
		 							GW_Float dC, GW_Vector3D EdgeAB, GW_Vector3D EdgeAC, GW_Vector3D EdgeBC,
									GW_Vector3D EdgeAD, GW_Vector3D EdgeCD, GW_Float ab, GW_Float ac, GW_Float bc, 
									GW_Float ad, GW_Float bd, GW_Float cd, GW_Float F );

	static GW_Float ComputeUpdateFace( GW_Float dA, GW_Float dB, GW_Vector3D EdgeAB, GW_Vector3D EdgeAC, 
									GW_Float ab, GW_Float ac, GW_Float bc, GW_Float F );

};

/*------------------------------------------------------------------------------*/
/** \name a vector of FM_propagationMesh */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class FM_propagationMesh*> T_propagationMeshVector;
typedef T_propagationMeshVector::iterator IT_propagationMeshVector;
typedef T_propagationMeshVector::reverse_iterator RIT_propagationMeshVector;
typedef T_propagationMeshVector::const_iterator CIT_propagationMeshVector;
typedef T_propagationMeshVector::const_reverse_iterator CRIT_propagationMeshVector;
//@}

} // End namespace GW

#endif // _FM_PROPAGATIONMESH_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Birgit Stender
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
