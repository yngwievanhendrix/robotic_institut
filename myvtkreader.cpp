#include "myvtkreader.h"

myVtkReader::myVtkReader()
{
    maxReadScalarValue = 0.0;
}

double myVtkReader::getMaxReadScalarValue(){
    return maxReadScalarValue;
}

int myVtkReader::getTotalPoints(){
    return this->total_points;
}

void myVtkReader::readScalars(vtkDoubleArray *inArray){
    QFile file("C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\Scalars.txt");
    QString line;
    file.open(QIODevice::ReadOnly);
    QTextStream ts(&file);
    int tempDouble;
    line = ts.readLine();
    while(line != "END"){
        tempDouble = line.toDouble();
        if(tempDouble > maxReadScalarValue){
            maxReadScalarValue = tempDouble;
        }
        inArray->InsertNextValue(tempDouble);
        line = ts.readLine();
    }
}

void myVtkReader::readMesh(vtkUnstructuredGrid *inGrid,vtkPoints *points){

    QFile file("C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\phantom_vtkMesh");
    QString line;
    QStringList list;
    file.open(QIODevice::ReadOnly);
    QTextStream ts(&file);
    bool read = true;
    this->total_points;

    /* Reading untill POINTS */
    while(read){
        line = ts.readLine();
        list = line.split(" ");
        if(list.at(0) == "POINTS"){
            read = false;
            total_points = QString(list.at(1)).toInt();

        }
    }
    /* Read the Points */
    for(int i=0;i<total_points;i++){
        line = ts.readLine();
        list = line.split(" ");
        points->InsertNextPoint(QString(list.at(0)).toDouble(),QString(list.at(1)).toDouble(),QString(list.at(2)).toDouble());
    }
    inGrid->SetPoints(points);

    /* Reading line that start with 3 and ignoring them*/
    read = true;
    while(read){
        line = ts.readLine();
        list = line.split(" ");
        if(list.at(0) == "4"){
            read = false;

            vtkIdType ptIds[] = {QString(list.at(1)).toInt(),
                             QString(list.at(2)).toInt(),
                         QString(list.at(3)).toInt(),
                     QString(list.at(4)).toInt()};
            inGrid->InsertNextCell(VTK_TETRA,4,ptIds);
        }
    }

    /* Reading the line that start with 4*/
    read = true;
    int count = 0;
    while(read){
        line = ts.readLine();
        list = line.split(" ");
        if(list.at(0) != "4"){
        //if(count == 350){
            read = false;
            break;
        }
        else {
            vtkIdType ptIds[] = {QString(list.at(1)).toInt(),
                             QString(list.at(2)).toInt(),
                         QString(list.at(3)).toInt(),
                     QString(list.at(4)).toInt()};
             inGrid->InsertNextCell(VTK_TETRA,4,ptIds);
        }
        count++;
    }
}
