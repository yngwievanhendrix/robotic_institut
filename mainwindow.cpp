#include "mainwindow.h"
#include "myvtkreader.h"
#include "mypropagation.h"
//#include "myvtkinteractor.h"

#include "ui_mainwindow.h"

#include <QFile>
#include <QString>
#include <QTextStream>
#include <QStringList>

#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkTetra.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkPointData.h"
#include "vtkLookupTable.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRendererCollection.h"

#include "myvtkanimator.h"

MainWindow::MainWindow(QWidget *parent) :QMainWindow(parent),ui(new Ui::MainWindow) {
    
    QString tString;
    ui->setupUi(this);
    //setup time
    this->qTimerTime = 100;
	//setup selectedPoint
	this->selectedPointId = -1;


    //Declare RenderWindow
    renWin = vtkRenderWindow::New();

    
    myReader = myVtkReader();
    myAnimator = myVtkAnimator();
	//Aufruf von dem leeren myPropagation() konstruktor
	//verursacht die Haltepunkte 
	//myPropagate = myPropagation();

	//Declare myGeomety
    myGeometry = vtkUnstructuredGrid::New();
    
	//Declare points

    points = vtkPoints::New();
    

    //Read the Mesh from a file that in the Class specified
    myReader.readMesh(myGeometry,points);

	//Creating Point locator
	pointLocator = vtkKdTreePointLocator::New();
	pointLocator->SetDataSet(myGeometry);
	pointLocator->BuildLocator();

	this->ui->animButton->setDisabled(true);
	this->ui->mexFunctionButton->setDisabled(true);
	this->ui->resetButton->setDisabled(true);

	int nVertex = this->myGeometry->GetNumberOfPoints();	
	int nVolumes = this->myGeometry->GetNumberOfCells();
	int nFaces = nVolumes*4;
	int nStartPoints = 1;
	int nEndPoints = 1;
	myPropagate.setDimensions(nVertex,nFaces,nVolumes,nStartPoints,nEndPoints);
    
	int inStartPoint[] = {1609};
	vtkPoints *inVertices = vtkPoints::New();
	inVertices = this->myGeometry->GetPoints();
	myPropagate.setGeometry(inStartPoint,inVertices,this->myGeometry);


	//myReader.readScalars(myAnimator.getRefAnimArray());
	
	myAnimator.setTotalPoints(myReader.getTotalPoints());
	myAnimator.setDefRefAnimArray(1.0);
	myAnimator.setMaxTempScalarValue(1.0);
    
    // das soll Protected sein wird irgendwann intern ausgeführt
    myAnimator.reset();


    //Output the max. read Scalar value
    //tString.setNum(myAnimator.getMaxScalarValue());
    //ui->textEdit->append(tString);
    //assign read Scalars to myGeometry
    myGeometry->GetPointData()->SetScalars(myAnimator.getVarAnimArray());

    //myGeometry->GetPointData()->GetScalars()->SetLookupTable();
    // LookUpTable
   
	lut = vtkLookupTable::New();
    //lut->SetNumberOfTableValues((int)myReader.getMaxReadScalarValue());
	lut->SetNumberOfTableValues(1);
    lut->Build();




    /* Displaying the result */

    mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection(myGeometry->GetProducerPort());
    //mapper->SetScalarRange(0,(int)myReader.getMaxReadScalarValue());
	mapper->SetScalarRange(0,1);
    mapper->SetLookupTable(lut);



    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);

    vtkRenderer *renderer = vtkRenderer::New();
    renderer->AddActor(actor);

    this->ui->qvtkWidget->SetRenderWindow(renWin);
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(renderer);


    /*Ende Testen*/
	
    vtkEventQtSlotConnect *myConnection = vtkEventQtSlotConnect::New();
    myConnection->Connect(ui->qvtkWidget->GetRenderWindow()->GetInteractor(),
                          vtkCommand::LeftButtonPressEvent,
                          this,
                          SLOT(onClickedMessageFromMyVtkInteractor(vtkObject*)));
    QObject::connect(this->ui->animButton,SIGNAL(clicked()),this,SLOT(onClickedAnimateButton()));
    QObject::connect(&timer,SIGNAL(timeout()),this,SLOT(onTimerTimeOut()));
    QObject::connect(this->ui->resetButton,SIGNAL(clicked()),this,SLOT(onClickedResetButton()));
	QObject::connect(this->ui->mexFunctionButton,SIGNAL(clicked()),this,SLOT(onClickedMexFunctionButton()));
    int TSPos = this->ui->timeSlider->sliderPosition();
    tString.setNum(TSPos);
    this->ui->textEdit->setText(tString);
	
    //this->ui->pushButton->;
}

MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::onClickedAnimateButton(){
	
    //int pos = this->ui->timeSlider->pos().x;
    //int pos = this->ui->timeSlider->cursor().pos().x();

    QString tString;
    int pos = this->ui->timeSlider->sliderPosition();
    tString.setNum(pos*this->qTimerTime);

    this->ui->textEdit->setText(tString);
    myAnimator.setProcessingTime(pos*this->qTimerTime);
    this->timer.start(qTimerTime);
}

void MainWindow::onClickedMessageFromMyVtkInteractor(vtkObject *inObj){
	if(this->ui->checkBox->isChecked()) {
		vtkRenderWindowInteractor *tIntActor = vtkRenderWindowInteractor::SafeDownCast(inObj);
		int tEvent_pos[2];
		tIntActor->GetEventPosition(tEvent_pos);
		tIntActor->GetPicker()->Pick(tEvent_pos[0],tEvent_pos[1],0,
			this->ui->qvtkWidget->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
	
		double tPicked_point[3];
		tIntActor->GetPicker()->GetPickPosition(tPicked_point);
		this->selectedPointId = this->pointLocator->FindClosestPoint(tPicked_point);
		this->ui->checkBox->setChecked(false);
		this->ui->animButton->setDisabled(false);
		this->ui->mexFunctionButton->setDisabled(false);
		this->ui->resetButton->setDisabled(false);

		int inStartPoint[] = {0};
		inStartPoint[0] = this->selectedPointId;
		vtkPoints *inVertices = vtkPoints::New();
		inVertices = this->myGeometry->GetPoints();
		myPropagate.setGeometry(inStartPoint,inVertices,this->myGeometry);

		QString tString;
		tString.setNum(this->selectedPointId);
		tString = "selected Point ID:" + tString;
		this->ui->textEdit->setText(tString);
	}
}
void MainWindow::onTimerTimeOut(){

    QString tString1;
    QString tString2;
    tString1.setNum(this->myAnimator.getProcessedPoints());
    tString1.append(";");
    tString2.setNum(this->myAnimator.getTotalPoints());
    tString1.append(tString2);
    this->ui->textEdit_2->append(tString1);

	if((this->myAnimator.getProcessedPoints() < this->myAnimator.getTotalPoints()) && ((this->myAnimator.getProcessedTime()/1000.0) < this->myAnimator.getMaxScalarValue())){
        this->myAnimator.increaseProcessedTime();
        this->myAnimator.updateVarAnimArray();
        myGeometry->GetPointData()->SetScalars(myAnimator.getVarAnimArray());
        myGeometry->Modified();
        this->ui->qvtkWidget->GetRenderWindow()->Render();



        tString1.setNum(this->myAnimator.getProcessedTime()/1000.0);
        tString1.append("/");
		tString2.setNum(this->myAnimator.getMaxScalarValue());
		tString1.append(tString2);
        this->ui->textEdit->setText(tString1);
        timer.start(qTimerTime);
    }
    else{
        timer.stop();

    }
}

void MainWindow::onClickedResetButton(){
    myAnimator.reset();
    myGeometry->GetPointData()->SetScalars(myAnimator.getVarAnimArray());
    myGeometry->Modified();
    this->ui->qvtkWidget->GetRenderWindow()->Render();
}

void MainWindow::on_timeSlider_sliderMoved(int position)
{
    QString tString;
    tString.setNum(position*this->qTimerTime);
    myAnimator.setProcessingTime(position*this->qTimerTime);
    //this->ui->textEdit->setText(tString);
}

void MainWindow::onClickedMexFunctionButton()
{
	
	QString tString;
	myPropagate.mexFunction();
	this->ui->textEdit_2->setText("mexFunction done!");
	myPropagate.saveDistances("C:\\Dokumente und Einstellungen\\awada.ROB\\Desktop\\vtk_Reader_extended_1\\D.txt");
	this->myAnimator.setRefAnimArray(myPropagate.getDistances());
	this->ui->textEdit->setText("RefAnim done!");
	this->mapper->SetScalarRange(0,myAnimator.getMaxScalarValue());
	//->SetScalarRange(0,(int)myAnimator.getMaxScalarValue());
	this->lut->SetNumberOfTableValues((int)myAnimator.getMaxScalarValue());
    this->lut->Build();
	tString.setNum((int)myAnimator.getMaxScalarValue());
	this->ui->textEdit->setText(tString);
	this->myGeometry->GetPointData()->SetScalars(myAnimator.getVarAnimArray());
	this->myGeometry->Modified();
}


//git@github.com:yngwievanhendrix/RoboticProject.git