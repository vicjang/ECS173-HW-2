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

    // The following part is basically what was in the original
    // SquaredEuclideanDistanceTo function
    int dimension = p1.GetPointDimension();
    double sum = 0;
    for ( unsigned int i = 0; i < dimension; i++ )
        {
        const double component = (double)( p2[i] );
        const double difference = (double)( p1[i] ) - component;
        sum += difference * difference;
        }

    // Truncated distance function
    // value 700 was determined through observation of sum values printed.
    return ( sum > 700 ) ? 700 : sum;


}

//    ------------   Your code for question 3 here ---------------

template <class TFixedPointSet, class TMovingPointSet>
void
IterativeClosestPointMetricIndexed<TFixedPointSet,TMovingPointSet>
::ClosestPointFunctionYOURS::Evaluate(FixedPointSetConstPointer pointSet,const typename Superclass::OutputPointType &p,CPInfo& cp) {

    // K is the number of points on the fixed mesh to be taken into consideration for robust point matching
    const int K = 3;
    int it = 0;
    std::list< double >                   shortestDistances;
    std::list< FixedPointSetPointType >   thePoints;
    shortestDistances.push_back( 0 );
    shortestDistances.back();

// This part is the code from the original closest point function
   // Go trough the list of fixed point and find the closest distance
    PointIterator pointItr = pointSet->GetPoints()->Begin();
    PointIterator pointEnd = pointSet->GetPoints()->End();
    cp.distance = NumericTraits<double>::max();
    cp.fixedPoint=pointItr.Value();

    if (!this->m_DistanceFunction)
    { 
        itk::ExceptionObject exception(__FILE__, __LINE__);
        char msg[1024];
        sprintf(msg,"Distance Function not set.");
        exception.SetDescription(msg);
        char loc[1024];
        sprintf(loc,"IterativeClosestPointMetricIndexed::ClosestPointFunctionVertex::Evaluate.");
        exception.SetLocation(loc);
        throw exception;
    }
// end of original code

    double temp_min_distance = NumericTraits<double>::max();
    FixedPointSetPointType temp_point;


    cp.movingPoint = p;
    while( it < K )
    {
        while( pointItr != pointEnd )
        {
            double dist = this->m_DistanceFunction->Evaluate( pointItr.Value(), p );
        
            if( ( dist < temp_min_distance ) && ( dist > shortestDistances.back() + 0.001 ) )
            {
                temp_min_distance = dist;
                temp_point = pointItr.Value();
            }
            pointItr++;
        }

        // when the minimum is found, store it in a vector
        shortestDistances.push_back( temp_min_distance );
        thePoints.push_back( temp_point );


        // reset the varaibl and restart from beginning until we have the K minimums
        temp_min_distance = NumericTraits<double>::max();
        pointItr = pointSet->GetPoints()->Begin();
        it++;
    }


    // pop the 0 that was initially pushed for iteratively finding minimums
    shortestDistances.pop_front();
    std::list< double >::iterator dis_iter = shortestDistances.begin();
    typename std::list< FixedPointSetPointType >::iterator pt_iter =  thePoints.begin();

    double sum = 0;
    itk::Vector<float, 3> vec;
    vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 0.0;

    for( ; dis_iter != shortestDistances.end() ; dis_iter++, pt_iter++ )
    {
        // convert itk::Point to itk::Vector to perform mathematical computation
        vec += ( pt_iter->GetVectorFromOrigin() ) * ( 1 / sqrt(*dis_iter) );
        sum += ( 1 / sqrt(*dis_iter) );
    }
    // the closest point afte RPM
    vec = vec * ( 1 / sum );

    FixedPointSetPointType matched_point;
    matched_point[0] = 0.0;
    matched_point[1] = 0.0;
    matched_point[2] = 0.0;
    // convert vector back to point by adding a vector from origin
    matched_point += vec;

    // store final return value
    cp.distance = this->m_DistanceFunction->Evaluate( matched_point, p );
    cp.fixedPoint = matched_point;

}

} // end namespace itk


#endif
