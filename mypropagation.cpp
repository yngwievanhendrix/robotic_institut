#include "mypropagation.h"
#include <stdio.h>
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "fm/fm_propagation/FM_propagationMesh.h"
#include "fm/fm_propagation/FM_propagationVertex.h"
using namespace FM;
using namespace GW;



myPropagation::myPropagation(){
}

myPropagation::myPropagation(int nVertex,int nFaces,int nVolumes,int nStartPoints,int nEndPoints){
	this->nverts = nVertex;
	this->vertex = (double*)malloc(sizeof(double) * (3*(nVertex+1)));
	this->nfaces = nFaces;
	this->faces = (int*)malloc(sizeof(int) * (3*(nFaces+1)));
	this->nvolumes = nVolumes;
	this->volumes = (int*)malloc(sizeof(int) * (4*(nVolumes+1)));
	this->nstart = nStartPoints;
	this->start_points = (int*)malloc(sizeof(int) * (nStartPoints+1));
	this->nend = nEndPoints;
	this->end_points = (int*)malloc(sizeof(int) * (nEndPoints+1));
	nbr_iter = 0;
	niter_max = -1;
	dmax = 1e9;


	this->D = (double*)malloc(sizeof(double) * (1+nVertex));
	this->S = (double*)malloc(sizeof(double) * (1+nVertex));
	this->Q = (double*)malloc(sizeof(double) * (1+nVertex));

	Ww = (double*)malloc(sizeof(double) * (nVertex+1));
	for(int i=0;i<nVertex;i++){
		Ww[i] = 1;
	}
}

void myPropagation::setDimensions(int nVertex,int nFaces,int nVolumes,int nStartPoints,int nEndPoints){
	this->nverts = nVertex;
	this->vertex = (double*)malloc(sizeof(double) * (3*(nVertex+1)));
	this->nfaces = nFaces;
	this->faces = (int*)malloc(sizeof(int) * (3*(nFaces+1)));
	this->nvolumes = nVolumes;
	this->volumes = (int*)malloc(sizeof(int) * (4*(nVolumes+1)));
	this->nstart = nStartPoints;
	this->start_points = (int*)malloc(sizeof(int) * (nStartPoints+1));
	this->nend = nEndPoints;
	this->end_points = (int*)malloc(sizeof(int) * (nEndPoints+1));
	nbr_iter = 0;
	niter_max = -1;
	dmax = 1e9;


	this->D = (double*)malloc(sizeof(double) * (1+nVertex));
	this->S = (double*)malloc(sizeof(double) * (1+nVertex));
	this->Q = (double*)malloc(sizeof(double) * (1+nVertex));

	Ww = (double*)malloc(sizeof(double) * (nVertex+1));
	for(int i=0;i<nVertex;i++){
		Ww[i] = 1;
	}
}


myPropagation::~myPropagation(){
	free(this->vertex);
	free(this->faces);
	free(this->volumes);
	free(this->start_points);
	free(this->end_points);
}

//myPropagation::setGeometry(int *inStartPoints, vtkPoints *inVertices, vtkCellArray *inVolumes){
void myPropagation::setGeometry(int *inStartPoints, vtkPoints *inVertices, vtkUnstructuredGrid *inGeometry){
	int i=0;
	int j=0;
	for(i=0;i<this->nverts;i++){
		this->vertex[j++] = inVertices->GetPoint(i)[0];
		this->vertex[j++] = inVertices->GetPoint(i)[1];
		this->vertex[j++] = inVertices->GetPoint(i)[2];
	}
	i=0;
	j=0;
	for(i=0;i<this->nstart;i++){
		this->start_points[i] = inStartPoints[i];
	}
	i=0;

	for(i=0;i<this->nvolumes;i++){
		this->volumes[j++] = inGeometry->GetCell(i)->GetPointIds()->GetId(0);
		this->volumes[j++] = inGeometry->GetCell(i)->GetPointIds()->GetId(1);
		this->volumes[j++] = inGeometry->GetCell(i)->GetPointIds()->GetId(2);
		this->volumes[j++] = inGeometry->GetCell(i)->GetPointIds()->GetId(3);
	}
	saveVertex();
	saveVolumes();
	saveStartPoints();

}

void myPropagation::saveVertex(){
	QString tFileString = "C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\Vertex.txt";
	QFile file(tFileString);
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	QString tString;
	int i=0;
	for(i=0;i<this->nverts;i++){
		out << i << ":";
		out << vertex_(0,i);
		out << ";";
		out << vertex_(1,i);
		out << ";";
		out << vertex_(2,i);		
		out << "\n";
	}
	file.close();
}
void myPropagation::saveVolumes(){
	QString tFileString = "C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\Volumes.txt";
	QFile file(tFileString);
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	QString tString;
	int i=0;
	for(i=0;i<this->nvolumes;i++){
		out << i << ":";
		out << volumes_(0,i);
		out << ";";
		out << volumes_(1,i);
		out << ";";
		out << volumes_(2,i);		
		out << ";";
		out << volumes_(3,i);
		out << "\n";
	}
	file.close();
}
void myPropagation::saveStartPoints(){
	QString tFileString = "C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\StartPoints.txt";
	QFile file(tFileString);
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	QString tString;
	int i=0;
	for(i=10;i<this->nstart;i++){
		out << i << ":";
		out << start_points[i];		
		out << "\n";
	}
	file.close();
}

inline void display_message(const char* mess, int v)
{
	char str[128];
	sprintf(str, mess, v);
	//mexWarnMsgTxt(str);
}

GW_Float WeightCallback(FM_propagationVertex& Vert)
{
	GW_U32 i = Vert.GetID();
	return Ww[i];
}

GW_Bool StopMarchingCallback( FM_propagationVertex& Vert )
{
	// check if the end point has been reached
	GW_U32 i = Vert.GetID();
	//	display_message("ind %d",i );
	//	display_message("dist %f",Vert.GetDistance() );
	if( Vert.GetDistance()>dmax )
		return true;
	//for( int k=0; k<nend; ++k )
	//	if( end_points[k]==i )
	//		return true;
	return false;
}

GW_Bool InsersionCallback( FM_propagationVertex& Vert, GW_Float rNewDist )
{
	// check if the distance of the new point is less than the given distance
	GW_U32 i = Vert.GetID();
	bool doinsersion = nbr_iter<=niter_max;
	if( L!=NULL )
		doinsersion = doinsersion && (rNewDist<L[i]);
	nbr_iter++;
	return doinsersion;
}
GW_Float HeuristicCallback( FM_propagationVertex& Vert )
{
	// return the heuristic distance
	GW_U32 i = Vert.GetID();
	return H[i];
}


void myPropagation::mexFunction() {

	Mesh.SetNbrVertex(nverts);
	for( int i=0; i<nverts; ++i )
	{
		FM_propagationVertex& vert = (FM_propagationVertex&) Mesh.CreateNewVertex();
		vert.SetPosition( GW_Vector3D(vertex_(0,i),vertex_(1,i),vertex_(2,i)) );
		for( int k=0; k<nend; ++k )
			if( end_points[k]==i )
				vert.SetStoppingVertex(GW_True);
		Mesh.SetVertex(i, &vert);

	}
	Mesh.SetNbrVolume(nvolumes);
	for( int i=0; i<nvolumes; ++i )
	{
		FM_propagationVolume& volume = (FM_propagationVolume&) Mesh.CreateNewVolume();
		FM_volumeVertex* v1 = Mesh.GetVertex((int) volumes_(0,i)); GW_ASSERT( v1!=NULL );
		FM_volumeVertex* v2 = Mesh.GetVertex((int) volumes_(1,i)); GW_ASSERT( v2!=NULL );
		FM_volumeVertex* v3 = Mesh.GetVertex((int) volumes_(2,i)); GW_ASSERT( v3!=NULL );
		FM_volumeVertex* v4 = Mesh.GetVertex((int) volumes_(3,i)); GW_ASSERT( v4!=NULL );
		volume.SetVertex( *v1,*v2,*v3, *v4 );
		Mesh.SetVolume(i, &volume);
		v1->SetVolume(volume);
		v2->SetVolume(volume);
		v3->SetVolume(volume);
		v4->SetVolume(volume);
	}

	Mesh.BuildConnectivity();
	//Mesh.CheckIntegrity();

	//// set up fast marching
	Mesh.ResetPropagationMesh();
	for( int i=0; i<nstart; ++i )
	{
		FM_propagationVertex* v = (FM_propagationVertex*) Mesh.GetVertex((GW_U32) start_points[i]);
		GW_ASSERT( v!=NULL );
		Mesh.AddStartVertex( *v );
		// initialize the distance of the starting points if available
		if( values!=NULL)
			v->SetDistance(values[i]);
	}

	//QtConcurrent::run(this, &CoreData::afunc);
	Mesh.SetUpFastMarching();
	Mesh.RegisterWeightCallbackFunction( WeightCallback );
	Mesh.RegisterForceStopCallbackFunction( StopMarchingCallback );
	Mesh.RegisterVertexInsersionCallbackFunction( InsersionCallback );
	if( H!=NULL )
		Mesh.RegisterHeuristicToGoalCallbackFunction( HeuristicCallback );

	// perform fast marching
	//	display_message("itermax=%d", niter_max);
	Mesh.PerformFastMarching();

	// output result
	for( int i=0; i<nverts; ++i )
	{
		FM_propagationVertex* v = (FM_propagationVertex*) Mesh.GetVertex((GW_U32) i);
		GW_ASSERT( v!=NULL );
		D[i] = v->GetDistance();
		S[i] = v->GetState();
		FM_propagationVertex* v1 = v->GetFront();
		if( v1==NULL )
			Q[i] = -1;
		else
			Q[i] = v1->GetID();
	}
	//free(this->vertex);
	//free(this->faces);
	//free(this->volumes);
	//free(this->start_points);
	//free(this->end_points);

}

void myPropagation::saveDistances(QString tFileString){
	QFile file(tFileString);
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	QString tString;
	int i=0;
	for(i=0;i<this->nverts;i++){
		out << i << " ";
		out << this->D[i];
		out << "\n";
	}
	out << "END";
	file.close();

}

double *myPropagation::getDistances(){
	return this->D;
}