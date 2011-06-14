#ifndef MYPROPAGATION_H
#define MYPROPAGATION_H

#include <math.h>
#include "config.h"

#include <algorithm>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>

using std::string;
using std::cerr;
using std::cout;
using std::endl;


//#include "mex.h"
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "fm/fm_propagation/FM_propagationMesh.h"
#include "fm/fm_propagation/FM_propagationVertex.h"

//Vtk libs
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkCell.h"
#include "vtkIdList.h"

#include <QFile>
#include <QTextStream>
#include <QString>

#define faces_(k,i) faces[k+3*i]
#define vertex_(k,i) vertex[k+3*i]
#define volumes_(k,i) volumes[k+4*i]

using namespace GW;
using namespace FM;
class myPropagation{
private:
    double* vertex;
    int nverts;
    int* faces;
    int nfaces;
    int* volumes;
    int nvolumes;
    int* start_points;
    int nstart;
    int* end_points;
    int nend;
//    double* H;	// heuristic
//    double* L;	// bound on current distance
    //double* Ww;	// weight
//    int niter_max;
//    double dmax;
    double* values;
    // outputs
    double* D;	// distance
    double* S;	// state
    double* Q;	// nearest neighbor

	FM_propagationMesh Mesh;

//    int nbr_iter;


public:
	myPropagation();
    myPropagation(int,int,int,int,int);
    ~myPropagation();
	void setDimensions(int,int,int,int,int);
	void saveDistances(QString);
	void saveVertex();
	void saveVolumes();
	void saveStartPoints();
    //void setGeometry(int*,vtkPoints*,vtkCellArray*);
    void setGeometry(int*,vtkPoints*,vtkUnstructuredGrid*);
//    inline void display_message(const char*,int);
//    GW_Float WeightCallback(FM_propagationVertex&);
//    GW_Bool StopMarchingCallback(FM_propagationVertex&);
//    GW_Bool InsersionCallback(FM_propagationVertex&,GW_Float);
//    GW_Float HeuristicCallback(FM_propagationVertex&);
    void mexFunction();
	double *getDistances();


};
static double* Ww;	// weight: ja
static double dmax; //ja
static int nbr_iter; //ja
static int niter_max; // 
static double* H;	// heuristic :ja
static double* L;	// bound on current distance : ja
inline void display_message(const char*,int);
static GW_Float WeightCallback(FM_propagationVertex&);
static GW_Bool StopMarchingCallback(FM_propagationVertex&);
static GW_Bool InsersionCallback(FM_propagationVertex&,GW_Float);
static GW_Float HeuristicCallback(FM_propagationVertex&);

#endif
