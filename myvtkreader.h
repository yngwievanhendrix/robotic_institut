#ifndef MYVTKREADER_H
#define MYVTKREADER_H

#include <QFile>
#include <QString>
#include <QTextStream>
#include <QStringList>
#include "stdio.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkTetra.h"
#include "vtkDoubleArray.h"

class myVtkReader
{
private:
    double maxReadScalarValue;
    int total_points;

public:
    myVtkReader();
    double getMaxReadScalarValue();
    int getTotalPoints();
    void readScalars(vtkDoubleArray *inArray);
    void readMesh(vtkUnstructuredGrid *inGrid,vtkPoints *inPoints);
};


#endif // MYVTKREADER_H
