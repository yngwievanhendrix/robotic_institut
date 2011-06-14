#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include "vtkEventQtSlotConnect.h"
#include "vtkCommand.h"
#include "vtkRenderWindow.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "myvtkanimator.h"
#include "myvtkreader.h"
#include "mypropagation.h"
#include "vtkDataSetMapper.h"
#include "vtkPointPicker.h"
#include "vtkKdTreePointLocator.h"

#define TIMER_TIME 100

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    vtkRenderWindow *renWin;
    vtkUnstructuredGrid *myGeometry;
    vtkPoints *points;
    vtkLookupTable *lut;
	vtkDataSetMapper *mapper;
	vtkKdTreePointLocator *pointLocator;
	int selectedPointId;

    int qTimerTime;
    QTimer timer;

    myVtkReader myReader;
    myVtkAnimator myAnimator;
	myPropagation myPropagate;

public slots:
    void onClickedMessageFromMyVtkInteractor(vtkObject*);
    void onClickedAnimateButton(void);
    void onClickedResetButton(void);
    void onTimerTimeOut(void);
	void onClickedMexFunctionButton(void);
protected:
    //vtkEventQtSlotConnect* vtkQtConnection;

private slots:
    
    void on_timeSlider_sliderMoved(int position);
};

#endif // MAINWINDOW_H
