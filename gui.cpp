 //This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.
//Author: Sean Mann

//CLA: mesh1 filename, mesh2 filename, correspondence file, render every nth cylinder

//VTK
#include <vtkVRMLImporter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkActor.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkLineSource.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkMath.h>
#include <vtkPiecewiseFunction.h>
#include <vtkAxesActor.h>

#include <fstream>
#include <math.h>
#include <vector>

#define DRAWAXES 1

void get_centroid_translation(vtkPolyData *points,vtkTransform *trans) {
  double c[3];
  c[0]=0;  c[1]=0;  c[2]=0;
  double np=points->GetNumberOfPoints();
  for (int i=0;i<np;i++) {
    double *p=points->GetPoint(i);
    c[0]+=p[0];    c[1]+=p[1];    c[2]+=p[2];
  }
  c[0]/=np;  c[1]/=np;  c[2]/=np;

  trans->Translate(-c[0],-c[1],-c[2]);
}

struct  pointWeight
{
  double weight;
  double point1[3];
  double point2[3];
};
 

int main(int argc, char *argv[]) 
{
  const char *mesh1File=argv[1];
  const char *mesh2File=argv[2];
  const char *corFile;
  int nthLine; 
  bool showLines;
  float r1,g1,b1,r2,g2,b2,r_bg,g_bg,b_bg;
  float lineHue=1;
  float lineValue=1;
  float lineWidth=5;

  // 1st data set color:
  r1=.5;  g1=0.7;  b1=0.5;
  // 2nd data set color:
  r2=1.0;  g2=1.0;  b2=0.0;
  // background color:
  r_bg=0.0;  g_bg=0.0;  b_bg=0.0;

  vtkRenderer *renderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  vtkDataArray *array = (vtkDataArray*) vtkDataArray::New();
  vtkDataArray *array2 = (vtkDataArray*) vtkDataArray::New();

  if (argc==3) {
    showLines=false;
  }
  else if (argc==5) {
    corFile = argv[3];
    nthLine = atoi(argv[4]);
    showLines=true;
  }  
  else {
    std::cerr << "Wrong number of arguments" << std::endl;
  }

  //read in mesh 1
  vtkVRMLImporter *vrmlImporter = vtkVRMLImporter::New();
  vrmlImporter->SetFileName(mesh1File);
  vrmlImporter->Update();
  vtkPolyDataMapper *pdm = vtkPolyDataMapper::New();
  pdm = (vtkPolyDataMapper*) vrmlImporter->GetVRMLDEFObject("_Faces");
  vtkPolyData* coords1 = pdm->GetInput();

  //read in mesh 2
  vtkVRMLImporter *vrmlImporter2 = vtkVRMLImporter::New();
  vrmlImporter2->SetFileName(mesh2File);
  vrmlImporter2->Update();
  vtkPolyDataMapper *pdm2 = vtkPolyDataMapper::New();
  pdm2 = (vtkPolyDataMapper*) vrmlImporter2->GetVRMLDEFObject("_Faces");
  vtkPolyData* coords2 = pdm2->GetInput();

 // translate them so the centroid of the first mesh is at the origin
  vtkTransform *centroid1 = vtkTransform::New();
  get_centroid_translation(coords1,centroid1);
  vtkTransformPolyDataFilter *translater1 = vtkTransformPolyDataFilter::New();
  translater1->SetInput( coords1 );
  translater1->SetTransform( centroid1 );
  translater1->Update();
  vtkPolyData* t1 = translater1->GetOutput();
  pdm->SetInput(t1);

  vtkTransformPolyDataFilter *translater2 = vtkTransformPolyDataFilter::New();
  translater2->SetInput( coords2 );
  translater2->SetTransform( centroid1 );
  translater2->Update();
  vtkPolyData* t2 = translater2->GetOutput();
  pdm2->SetInput(t2);

  //put mesh points into array
  array = pdm->GetInput()->GetPoints()->GetData();
  array2 = pdm2->GetInput()->GetPoints()->GetData();

  vtkActor *meshActor1 = vtkActor::New();
  vtkActor *meshActor2 = vtkActor::New();
  meshActor1->SetMapper(pdm);
  meshActor1->GetProperty()->SetColor(r1, g1, b1 );
  meshActor2->SetMapper(pdm2);
  meshActor1->GetProperty()->SetColor(r2, g2, b2 );
  renderer->SetBackground(r_bg, g_bg, b_bg);
  renderer->AddActor(meshActor1);
  renderer->AddActor(meshActor2);

#ifdef DRAWAXES
  vtkAxesActor *axes = vtkAxesActor::New();
  axes->SetTotalLength(100,100,100);
  axes->SetTipType(0);
  axes->SetShaftType(1);
  renderer->AddActor(axes);
#endif

  if (showLines) { 
    //read in correspondence file
    std::vector<pointWeight> pointWeights;
    std::fstream corr(corFile);
    double maxWeight=-1e+10;
    double minWeight=1e+10;

    while(corr.peek()!=std::ifstream::traits_type::eof()) {
      pointWeight pW;
      corr >> pW.point1[0];  
      corr >> pW.point1[1];  
      corr >> pW.point1[2];  
      corr >> pW.point2[0];  
      corr >> pW.point2[1];  
      corr >> pW.point2[2];  
      corr >> pW.weight;
      if (pW.weight<minWeight) { minWeight=pW.weight; }
      if (pW.weight>maxWeight) { maxWeight=pW.weight; }
      pointWeights.push_back(pW);
    }
    int numPoints=pointWeights.size();

    // set up color mapping for lines based on weights from file
    vtkPiecewiseFunction *lineColorFunction = vtkPiecewiseFunction::New();
    lineColorFunction->Initialize();
    lineColorFunction->AddSegment(minWeight,0,maxWeight,1);
    lineColorFunction->ClampingOn();
    vtkMath *math = vtkMath::New();

    double x, x1, x0, y, y1, y0, z, z1, z0;
    int pointIndex = 0; //used to index points, based on nthLine
    for (int i = 0; i < numPoints / nthLine; i++)
      {
	vtkActor *lineActor = vtkActor::New(); //allocate a new line
	vtkLineSource *line= vtkLineSource::New();
	vtkDataSetMapper *lineMapper = vtkDataSetMapper::New(); 
	vtkTransformFilter *transFilter = vtkTransformFilter::New();
	
	x0 = pointWeights[pointIndex].point1[0];
	y0 = pointWeights[pointIndex].point1[1];
	z0 = pointWeights[pointIndex].point1[2];

	x1 = pointWeights[pointIndex].point2[0];
	y1 = pointWeights[pointIndex].point2[1];
	z1 = pointWeights[pointIndex].point2[2];
	
	x = x1 - x0; y = y1 - y0; z = z1 - z0; //get differnece

	line->SetPoint1(x0,y0,z0);
	line->SetPoint2(x1,y1,z1);
	
	transFilter->SetTransform(centroid1);
	transFilter->SetInput(line->GetOutput());
	transFilter->Update();
	lineMapper->SetInput(transFilter->GetOutput()); 
	
	lineActor->SetMapper(lineMapper);
	float r,g,b;

	math->HSVToRGB(lineHue,lineColorFunction->GetValue(pointWeights[pointIndex].weight),lineValue,&r,&g,&b);
	lineActor->GetProperty()->SetColor(r,g,b);
	lineActor->GetProperty()->SetLineWidth(lineWidth);
	lineActor->GetProperty()->SetRepresentationToSurface();
	renderer->AddActor(lineActor);
	pointIndex += nthLine;
      }
  }


  renWin->AddRenderer(renderer);
  iren->SetRenderWindow(renWin);
  iren->SetInteractorStyle(style);
 
  renWin->SetSize(400, 400);
  vtkCamera *currentCamera = renderer->GetActiveCamera();
  currentCamera->Zoom(.8); 
  currentCamera->SetPosition(200,200,200);
  currentCamera->SetFocalPoint(0,0,0);
  currentCamera->SetClippingRange(.01,1000);
  renWin->Render();
  iren->Initialize();
  iren->Start();


  //Delete:
  vrmlImporter->Delete();
  vrmlImporter2->Delete();
  meshActor1->Delete();
  meshActor2->Delete();
  renderer->Delete();
  iren->Delete();
  style->Delete();
  renWin->Delete();
  return 0;
  
}
