/*=========================================================================
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkIterativeClosestPointMetric.h,v $
  Language:  C++
  Date:      $Date: 2004/04/28 20:06:26 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkIterativeClosestPointMetricIndexed_h
#define __itkIterativeClosestPointMetricIndexed_h

#include "itkPointSetToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkImage.h"
#include "itkArray.h"
#include <vector>
namespace itk
{
/**
itk::IterativeClosestPointMetricIndexed

This class computes the closest-point metric between two point sets a la the Iterative Closest Point (ICP).  The metric is flexible enough to allow the user to specify different functions for computing point-to-point distance and determining closest points.  It inherits from the itk::PointSetToPointSetMetric base class that comes with ITK.
 */

template < class TFixedPointSet, class TMovingPointSet >
class ITK_EXPORT IterativeClosestPointMetricIndexed : 
    public PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet>
{
public:

  /** Standard class typedefs. */
  typedef IterativeClosestPointMetricIndexed    Self;
  typedef PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet >  Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(IterativeClosestPointMetricIndexed, Object);
 
  /** Types transferred from the base class */
  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::TransformParametersType    TransformParametersType;
  typedef typename Superclass::TransformJacobianType      TransformJacobianType;

  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::FixedPointSetType          FixedPointSetType;

  typedef typename FixedPointSetType::PointType           FixedPointSetPointType;
  typedef typename Superclass::MovingPointSetType         MovingPointSetType;
  typedef typename MovingPointSetType::PointType           MovingPointSetPointType;

  typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
  typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;

  typedef typename Superclass::PointIterator              PointIterator;
  typedef typename Superclass::PointDataIterator          PointDataIterator;
  
  // The CPInfo structure stores closest-point information for one moving mesh
  // point and its closest fixed mesh point.  The positions of the fixed and 
  // moving points are stored as well as the computed distance between them.

  struct CPInfo {
    double distance;
    FixedPointSetPointType fixedPoint;
    typename Superclass::OutputPointType movingPoint;
  };

  //  DistanceFunction is the base class for classes that compute the distance between a point on the moving mesh and a point on the fixed mesh.  Its member function, Evaluate(), computes the distance.  Evaluate is purely virtual in the base class, meaning that the user must allocate a derived class to compute distances.
  class DistanceFunction
  {
  public:
    virtual double Evaluate(const FixedPointSetPointType&,const typename Superclass::OutputPointType&)=0;
  };
  
  // DistanceFunctionL2 computes the square of the traditional L2 (Euclidean) distance between the fixed and moving points
 class DistanceFunctionL2 : public DistanceFunction
  {
  public:
    virtual double Evaluate(const FixedPointSetPointType&,const typename Superclass::OutputPointType&);
  };

  // This is the DistanceFunction that you will define for Question 2 of the homework.  You must define the Evaluate() function to compute distances the way you want.  Feel free to add other member functions and/or state variables that will provide your DistanceFunction with parameters that will help it do what you want.
 class DistanceFunctionYOURS : public DistanceFunction
  {
  public:
    virtual double Evaluate(const FixedPointSetPointType&,const typename Superclass::OutputPointType&);
  };

// ClosestPointFunction is the base class for classes that determine, for a given moving point, its closest point on the fixed mesh.  Its function Evaluate() determines the closest fixed mesh point; in the base class, Evaluate is purely virutal, meaning that you need to allocate a derived class to determine closest points.  Note that m_DistanceFunction used by the ClosestPointFunction to determine closest points need not be the same as the m_DistanceFunction used to assign goodness-of-fit values to each moving point; see the assignment handout for an explanation of why that is.

  class ClosestPointFunction  {
  protected:
    DistanceFunction* m_DistanceFunction;
  public:
    virtual void Evaluate(FixedPointSetConstPointer,const typename Superclass::OutputPointType&,CPInfo&)=0;
    void SetDistanceFunction(DistanceFunction* df) { m_DistanceFunction=df; }
    DistanceFunction* GetDistanceFunction()  { return m_DistanceFunction; }
    ClosestPointFunction() { m_DistanceFunction=new DistanceFunctionL2; }
  };

  //  ClosestPointFunctionVertex implements the classical vertex-to-vertex closest point function used traditionally by ICP: the closest fixed mesh vertex is the closest fixed mesh surface point to the moving mesh point.

 class ClosestPointFunctionVertex : public ClosestPointFunction {
  public:
    virtual void Evaluate(FixedPointSetConstPointer,const typename Superclass::OutputPointType &,CPInfo&);
  };

  // ClosestPointFunctionYOURS is where you will write your code that will hopefully compute closest points a little smarter than ClosestPointFunctionVertex does.  You must implement your own Evaluate() function that does this.  Feel free to add other member functions or state variables that store parameters that help Evaluate() do the closest point calculations effectively.
 
 class ClosestPointFunctionYOURS : public ClosestPointFunction {
  public:
    virtual void Evaluate(FixedPointSetConstPointer,const typename Superclass::OutputPointType &,CPInfo&);
  };

  /** Get the number of values */
  unsigned int GetNumberOfValues() const;
  
  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & Derivative ) const;
  
  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  //  This function uses the ClosestPointFunction to compute, for each moving mesh vertex, its closest point on the fixed mesh.  These closest points are stored in a vector of CPInfo structures, one per moving mesh point

  void GetClosestPoints( const TransformParametersType & parameters,std::vector<CPInfo>&) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& Derivative ) const;

  // Set the ClosestPointFunction:

  void SetClosestPointFunction(ClosestPointFunction* cpf) { 
    m_ClosestPointFunction=cpf;
 }

  // Get the closest point function
  itkGetMacro(ClosestPointFunction,ClosestPointFunction*);

  // Get/Set the DistanceFunction used to assign goodness-of-fit values for each point
  void SetDistanceFunction(DistanceFunction* df) { 
    m_DistanceFunction=df; 
  }

  itkGetMacro(DistanceFunction,DistanceFunction*);

  // Get a const pointer to the transform that the metric uses to transform the moving mesh.  This is here so that the observer can print out the transformation parameters after each Metric evaluation.

  virtual const TransformType* GetTransformConst() const { return this->m_Transform; }

protected:
  IterativeClosestPointMetricIndexed();
  virtual ~IterativeClosestPointMetricIndexed() {};

  /** PrintSelf funtion */
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  IterativeClosestPointMetricIndexed(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ClosestPointFunction* m_ClosestPointFunction;  // Pointer to the ClosestPointFunction that finds closest fixed mesh points for each moving mesh point
  DistanceFunction* m_DistanceFunction;  // Computes distances between fixed and moving mesh points
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIterativeClosestPointMetricIndexed.txx"
#endif

#endif
