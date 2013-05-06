/*=========================================================================
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkIterativeClosestPointMetric.txx,v $
  Language:  C++
  Date:      $Date: 2004/04/28 20:06:26 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkIterativeClosestPointMetric_txx
#define _itkIterativeClosestPointMetric_txx

#include "itkIterativeClosestPointMetricIndexed.h"
#include "itkExceptionObject.h"

namespace itk
{

/** Constructor */
template <class TFixedPointSet, class TMovingPointSet> 
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::IterativeClosestPointMetricIndexed()
{
  m_ClosestPointFunction = 0;
  m_DistanceFunction = 0;
}

/** Return the number of values, i.e the number of points in the moving set */
template <class TFixedPointSet, class TMovingPointSet>  
unsigned int
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>  
::GetNumberOfValues() const
{
 MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

 if( !movingPointSet ) 
    {
    itkExceptionMacro( << "Moving point set has not been assigned" );
    }

 return  movingPointSet->GetPoints()->Size();
}

/** Get the match Measure */
template <class TFixedPointSet, class TMovingPointSet>  
typename IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>::MeasureType
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::GetValue( const TransformParametersType & parameters ) const
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet ) 
    {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
    }

  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

  if( !movingPointSet ) 
    {
    itkExceptionMacro( << "Moving point set has not been assigned" );
    }

  PointIterator pointItr = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd = movingPointSet->GetPoints()->End();

  MeasureType measure;
  measure.set_size(movingPointSet->GetPoints()->Size());

  this->SetTransformParameters( parameters );

  unsigned int id = 0;
  while( pointItr != pointEnd )
    {
    typename Superclass::InputPointType  inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint = 
      this->m_Transform->TransformPoint( inputPoint );
    CPInfo cp;
    if (m_ClosestPointFunction) {
      m_ClosestPointFunction->Evaluate(fixedPointSet,transformedPoint,cp);
    } 
    else {
      itkExceptionMacro( << "Closest Point Function not set.");
    }
    if (m_DistanceFunction) { 
      measure.put(id,m_DistanceFunction->Evaluate(cp.fixedPoint,transformedPoint));
    }
    else {
       itkExceptionMacro( << "Distance Function not set.");
    }
    ++pointItr;
    id++;
    }
  
  this->InvokeEvent( IterationEvent() );
  return measure;
}

template <class TFixedPointSet, class TMovingPointSet>  
void IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::GetClosestPoints( const TransformParametersType & parameters , std::vector<CPInfo>& cpInfos) const
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet ) 
    {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
    }

  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

  if( !movingPointSet ) 
    {
    itkExceptionMacro( << "Moving point set has not been assigned" );
    }

  PointIterator pointItr = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd = movingPointSet->GetPoints()->End();

  MeasureType measure;
  measure.set_size(movingPointSet->GetPoints()->Size());

  this->SetTransformParameters( parameters );

  unsigned int id = 0;
  while( pointItr != pointEnd )
    {
    typename Superclass::InputPointType  inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint = 
      this->m_Transform->TransformPoint( inputPoint );
    CPInfo cp;
    m_ClosestPointFunction->Evaluate(fixedPointSet,transformedPoint,cp);
    cp.distance=m_DistanceFunction->Evaluate(cp.fixedPoint,transformedPoint);
    cpInfos.push_back(cp);
    ++pointItr;
    id++;
    }
}

/** Get the Derivative Measure */
template <class TFixedPointSet, class TMovingPointSet>
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
                 DerivativeType & itkNotUsed(derivative) ) const
{

}

/** Get both the match Measure and theDerivative Measure  */
template <class TFixedPointSet, class TMovingPointSet>  
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::GetValueAndDerivative(const TransformParametersType & parameters, 
                        MeasureType & value, DerivativeType  & derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters,derivative);
}

/** PrintSelf method */
template <class TFixedPointSet, class TMovingPointSet>  
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

}
template <class TFixedPointSet, class TMovingPointSet>
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::ClosestPointFunctionVertex::Evaluate(FixedPointSetConstPointer pointSet,const typename Superclass::OutputPointType &p,CPInfo& cp) {
   // Go trough the list of fixed point and find the closest distance
    PointIterator pointItr = pointSet->GetPoints()->Begin();
    PointIterator pointEnd = pointSet->GetPoints()->End();
    cp.distance = NumericTraits<double>::max();
    cp.fixedPoint=pointItr.Value();

    if (!this->m_DistanceFunction) { 
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"Distance Function not set.");
      exception.SetDescription(msg);
      char loc[1024];
      sprintf(loc,"IterativeClosestPointMetricIndexed::ClosestPointFunctionVertex::Evaluate.");
      exception.SetLocation(loc);
      throw exception;
    }

    cp.movingPoint=p;
    while( pointItr != pointEnd )
      {
        double dist = this->m_DistanceFunction->Evaluate(pointItr.Value(),p);
	
        if(dist<cp.distance)
          {
	    cp.distance = dist;
	    cp.fixedPoint=pointItr.Value();
          }
        pointItr++;
      }
}

template <class TFixedPointSet, class TMovingPointSet>
double
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::DistanceFunctionL2::Evaluate(const FixedPointSetPointType& p1,const typename Superclass::OutputPointType& p2) {
  return p1.SquaredEuclideanDistanceTo(p2);
}

  //    ------------   Your code for question 2 here ---------------
template <class TFixedPointSet, class TMovingPointSet>
double
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::DistanceFunctionYOURS::Evaluate(const FixedPointSetPointType& p1,const typename Superclass::OutputPointType& p2) {

  // Put your code here for answering the second question on the assignment
}

//    ------------   Your code for question 3 here ---------------

template <class TFixedPointSet, class TMovingPointSet>
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::ClosestPointFunctionYOURS::Evaluate(FixedPointSetConstPointer pointSet,const typename Superclass::OutputPointType &p,CPInfo& cp) {
  // Put your code for answering the third question on the assignment here
}
} // end namespace itk


#endif
