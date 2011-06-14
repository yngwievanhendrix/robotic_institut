#ifndef MYVTKANIMATOR_H
#define MYVTKANIMATOR_H

#include "vtkPointData.h"
#include "vtkDoubleArray.h"

class myVtkAnimator
{
private:
    vtkDoubleArray *refAnimArray;
    vtkDoubleArray *varAnimArray;
    double maxReadScalarValue;
    int total_points;
    int processed_points;
    double processed_time;
    int processing_time;
public:
    myVtkAnimator();
    vtkDoubleArray *getRefAnimArray();

    void setMaxReadScalarValue(double);
    vtkDoubleArray *getVarAnimArray();
	void setRefAnimArray(double*);
    void setTotalPoints(int);
    int getTotalPoints();
    int getProcessedPoints();
    void resetVarAnimArray();
    void updateVarAnimArray();
    void increaseProcessedTime();
    double getProcessedTime();
    void reset();
	void updateMaxScalarValue();
    void setProcessingTime(int);
	double getMaxScalarValue();
	void setMaxTempScalarValue(double);
	void setDefRefAnimArray(double);
protected:

};

#endif // MYVTKANIMATOR_H
