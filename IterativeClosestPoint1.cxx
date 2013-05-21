/*=========================================================================
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: IterativeClosestPoint1.cxx,v $
  Language:  C++
  Date:      $Date: 2004/04/27 21:42:25 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifdef _WIN32
#pragma warning ( disable : 4786 )
#define _CRT_SECURE_NO_DEPRECATE
#endif

// Software Guide : BeginLatex
//
// This example illustrates how to perform Iterative Closest Point (ICP) 
// registration in ITK. The main class featured in this section is the 
// \doxygen{IterativeClosestPointMetric}.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet

// Standard Includes:
#include <iostream>
#include <fstream>

// ITK includes:
#include "itkVersorRigid3DTransform.h"
#include "itkIterativeClosestPointMetricIndexed.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSet.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkExceptionObject.h"
#include "itkCommand.h" 

// VTK includes:
#include <vtkVRMLImporter.h>
#include <vtkVRMLExporter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCell.h>

// Typedefs for the metric, transform, etc:

#define DIMENSION 3
typedef itk::PointSet< float, DIMENSION >   PointSetType;
typedef itk::IterativeClosestPointMetricIndexed< PointSetType,PointSetType> MetricType;
typedef PointSetType::PointType     PointType;
typedef PointSetType::PointsContainer  PointsContainer;
typedef MetricType::TransformType                 TransformBaseType;
typedef TransformBaseType::ParametersType         ParametersType;
typedef TransformBaseType::JacobianType           JacobianType;
typedef itk::VersorRigid3DTransform< double >      TransformType;
typedef itk::LevenbergMarquardtOptimizer OptimizerType;
typedef itk::PointSetToPointSetRegistrationMethod<PointSetType,PointSetType > RegistrationType;

// This function calls metric->GetClosestPoints() to compute all of the closest fixed mesh points for all of the moving mesh vertices.  Each of the closest point relationships is printed out to file.

void OutputCPFile(char *output_cp_fname,TransformType::Pointer transform, MetricType::Pointer  metric) {
  ofstream output_cp_file(output_cp_fname);
  
  // Write out the closest-point distance for each point
  std::vector<MetricType::CPInfo> cpInfos;
  metric->GetClosestPoints(transform->GetParameters(),cpInfos);
  std::vector<MetricType::CPInfo>::iterator cpInfos_i;
  for (cpInfos_i=cpInfos.begin();
       cpInfos_i!=cpInfos.end();
       cpInfos_i++) {
    MetricType::FixedPointSetPointType& fp=(*cpInfos_i).fixedPoint;
    MetricType::Superclass::OutputPointType& mp= (*cpInfos_i).movingPoint;
    for (int i=0;i<DIMENSION;i++) { output_cp_file << fp[i] << " ";  }
    for (int i=0;i<DIMENSION;i++) { output_cp_file << mp[i] << " ";  }
    output_cp_file << (*cpInfos_i).distance << std::endl;
  }
}

// This function uses the transform argument to transform a mesh of points stored in a vtkPolyData; the transformed mesh is then written out to file.

void OutputVRMLFile(const char *output_fname,vtkPolyData* coords,TransformType::Pointer transform) {

  std::ofstream output_vrml(output_fname);

  // VRML file preamble:
  output_vrml << "#VRML V2.0 utf8   \n\n        Shape { \n          geometry DEF _Faces IndexedFaceSet { \n\n  coord DEF _Points Coordinate {\n   point [\n";
  vtkPoints* points=coords->GetPoints();
  int np=coords->GetNumberOfPoints();

  double pf[DIMENSION];
  // Print out each point:

  for (int i=0;i<np;i++) {
    TransformType::InputPointType p,p2;
    points->GetPoint(i,pf);
    for (int d=0;d<DIMENSION;d++) { p[d]=pf[d]; }
    // Transform the point:
    p2=transform->TransformPoint(p);
    // Print it to file:
    for(int j=0;j<DIMENSION;j++) { output_vrml << p2[j] << " ";  }
    output_vrml <<  "," << std::endl;
  } 
  output_vrml << "              ]\n            }\n            coordIndex  [\n";

  // Print out mesh faces:

  int ncells=coords->GetNumberOfCells();
  for (int i=0;i<ncells;i++) {
    vtkCell *thisCell = coords->GetCell(i);
    if (thisCell->GetCellType() == VTK_TRIANGLE) {
      vtkIdList* pids=thisCell->GetPointIds();
      output_vrml << pids->GetId(0) << ", " << pids->GetId(1) << ", " << pids->GetId(2) << ", -1," << std::endl;
    }
   }
  
  output_vrml <<  "         ]\n          }\n        }\n";
}

// This computes the center of mass of a set of points:

void get_centroid(vtkPolyData *points, double *c, int Dim) {
  for (int i=0;i<Dim;i++) { c[i]=0; }
  double np=points->GetNumberOfPoints();

  for (int i=0;i<np;i++) {
    double *p=points->GetPoint(i);
    for (int j=0;j<Dim;j++) { c[j]+=p[j]; }
  }
  for (int j=0;j<Dim;j++) { c[j]/=np; }
}

// This class listens to the Metric and catches the IterationEvents that it throws.  When the Metric issues an IterationEvent, the Execute() function is called, the current transformation parameters are read from the metric, and are printed out.

class CommandIterationUpdate : public itk::Command 
{ 
public: 
  typedef CommandIterationUpdate Self; 
  typedef itk::Command Superclass; 
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self ); 
protected: 
  CommandIterationUpdate() {}; 
  typedef const MetricType *MetricPointer;
  
  virtual void Execute(itk::Object *caller, const itk::EventObject & event) 
  { 
    Execute( (const itk::Object *)caller, event); 
  } 

  virtual void Execute(const itk::Object * object, const itk::EventObject & event) 
  { 
    MetricPointer metric = dynamic_cast< MetricPointer >(object ); 
    if( typeid( event) !=typeid( itk::IterationEvent )) 
      { 
	return; 
      } 
    // print out the transformation parameters:
    std::cout << metric->GetTransformConst()->GetParameters() << std::endl; 
  }
};

int main(int argc, char * argv[] )
{ try {

    // Allocate the fixed and moving point sets:
  PointSetType::Pointer fixedPointSet  = PointSetType::New();
  PointSetType::Pointer movingPointSet = PointSetType::New();

  PointsContainer::Pointer fixedPointContainer  = PointsContainer::New();
  PointsContainer::Pointer movingPointContainer = PointsContainer::New();

  PointType fixedPoint;
  PointType movingPoint;

  // Use the vtkVRMLImporter to read the file containing coordinates of fixed points:
  const char *input1_fname=argv[1];
  vtkVRMLImporter *importer1 = vtkVRMLImporter::New();
  importer1->SetFileName( input1_fname );
  importer1->Update();
  vtkPolyDataMapper* pdm1 = (vtkPolyDataMapper*) importer1->GetVRMLDEFObject("_Faces");
  vtkPolyData* coords1 = pdm1->GetInput();
  double np1=coords1->GetNumberOfPoints();

  // Convert the fixed mesh points from the vtkPolyData representation to the FixedPointSet representation:

  for (int i=0;i<np1;i++) {
    MetricType::FixedPointSetPointType p; 
    p.CastFrom(itk::Point<double,DIMENSION>(coords1->GetPoint(i)));
    fixedPointContainer->InsertElement( i, p );
  }
  fixedPointSet->SetPoints( fixedPointContainer );
 std::cout << "Number of fixed Points = " << fixedPointSet->GetNumberOfPoints() << std::endl;

  // Use the vtkVRMLImporter to read the file containing coordinates of moving points:
  const char *input2_fname=argv[2];
  vtkVRMLImporter *importer2 = vtkVRMLImporter::New();
  importer2->SetFileName( input2_fname );
  importer2->Update();
  vtkPolyDataMapper* pdm2 = (vtkPolyDataMapper*) importer2->GetVRMLDEFObject("_Faces");
  vtkPolyData* coords2 = pdm2->GetInput();
  int np=coords2->GetNumberOfPoints();

  // Convert the moving mesh points from the vtkPolyData representation to the MovingPointSet representation:
  for (int i=0;i<np;i++) {
    MetricType::MovingPointSetPointType p; 
    p.CastFrom(itk::Point<double,DIMENSION>(coords2->GetPoint(i)));
    movingPointContainer->InsertElement( i, p );
  }
  movingPointSet->SetPoints( movingPointContainer );
 std::cout << "Number of moving Points = " << movingPointSet->GetNumberOfPoints() << std::endl;

//-----------------------------------------------------------
// Set up  the Metric
//-----------------------------------------------------------


  MetricType::Pointer  metric = MetricType::New();

  //  ----------  Here is how you would plug in your very own DistanceFunction to solve question 2 :

#if(1)
    MetricType::DistanceFunctionYOURS* distanceFunction = new MetricType::DistanceFunctionYOURS;
    metric->SetDistanceFunction(distanceFunction);
#else
    MetricType::DistanceFunctionL2* distanceFunction = new MetricType::DistanceFunctionL2;
    metric->SetDistanceFunction(distanceFunction);
#endif

  //  ----------  Here is how you would plug in your very own ClosestPointFunction to solve question 3 :

#if(1)
  MetricType::ClosestPointFunctionYOURS* closestPointFunction = new MetricType::ClosestPointFunctionYOURS;
  metric->SetClosestPointFunction(closestPointFunction);
#else
  MetricType::ClosestPointFunctionVertex* closestPointFunction = new MetricType::ClosestPointFunctionVertex;
  metric->SetClosestPointFunction(closestPointFunction);
#endif

//-----------------------------------------------------------
// Set up a Transform
//-----------------------------------------------------------

  TransformType::Pointer transform = TransformType::New();
  // these lines only really make sense for the VersorRigid3DTransform...
  double cent[3];  get_centroid(coords2,cent,DIMENSION);
  TransformType::InputPointType centPoint;
  for (int i=0;i<DIMENSION;i++) { centPoint[i]=cent[i]; }
  transform->SetCenter(centPoint);

  // Optimizer Type
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  optimizer->SetUseCostFunctionGradient(false);

  // Optimizer observer
  CommandIterationUpdate::Pointer observer =CommandIterationUpdate::New(); 
  metric->AddObserver( itk::IterationEvent(), observer ); 

  // Registration Method


  RegistrationType::Pointer   registration  = RegistrationType::New();

  // Scale the components of the Transform in the Optimizer
  OptimizerType::ScalesType scales( transform->GetNumberOfParameters() );
  scales.Fill( 1 );
  
  unsigned long   numberOfIterations =  1000;

  double          gradientTolerance  =  1e-6;   // convergence criterion
  double          valueTolerance     =  1e-6;   // convergence criterion
  double          epsilonFunction    =  1e-3;   // convergence criterion

  optimizer->SetScales( scales );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetValueTolerance( valueTolerance );
  optimizer->SetGradientTolerance( gradientTolerance );
  optimizer->SetEpsilonFunction( epsilonFunction );

  // Start the transformation from the initial guess provided on the command line:
  transform->SetIdentity();
  TransformType::OutputVectorType trans; 
  trans[0]=atof(argv[4]);  trans[1]=atof(argv[5]); trans[2]=atof(argv[6]);
  transform->SetTranslation(trans);
  TransformType::VersorType rotx,roty,rotz,rot; 
  rotx.SetRotationAroundX(atof(argv[7]) * 3.14159 / 180 );
  roty.SetRotationAroundY(atof(argv[8]) * 3.14159 / 180 );
  rotz.SetRotationAroundZ(atof(argv[9]) * 3.14159 / 180 );
  rot=rotx*roty*rotz;
  transform->SetRotation(rot);
  registration->SetInitialTransformParameters( transform->GetParameters() );

  //------------------------------------------------------
  // Connect all the components required for Registration
  //------------------------------------------------------
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetFixedPointSet( fixedPointSet );
  registration->SetMovingPointSet(   movingPointSet   );
  metric->SetFixedPointSet( fixedPointSet );
  metric->SetMovingPointSet( movingPointSet );
  metric->SetTransform( transform );

  // Output the initial vrml file and set of closest points ("before"):
  const char* output_base=argv[3];
  char output_vrml_name[1024],output_cp_name[1024];
  sprintf(output_vrml_name,"%s_before.wrl",output_base);
  OutputVRMLFile(output_vrml_name,coords2,transform);
  sprintf(output_cp_name,"%s_before.cp",output_base);
  OutputCPFile(output_cp_name,transform,metric);
  std::cout << optimizer << std::endl;

  // Run the registration:

  try 
    {
//    registration->StartRegistration();
    registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Solution = " << transform->GetParameters() << std::endl;
  std::cout << "Solution = " << transform << std::endl;

  // Output the final vrml file and set of closest points ("after"):
  sprintf(output_vrml_name,"%s_after.wrl",output_base);
  OutputVRMLFile(output_vrml_name,coords2,transform);
  sprintf(output_cp_name,"%s_after.cp",output_base);
  OutputCPFile(output_cp_name,transform,metric);
  
// Software Guide : EndCodeSnippet

  return EXIT_SUCCESS;

  } catch( itk::ExceptionObject & exp ) 
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << exp << std::endl;
    }
}

