/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaskedSpatialObjectToImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/07/19 15:21:11 $
  Version:   $Revision: 1.00 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskedSpatialObjectToImageFilter_h
#define __itkMaskedSpatialObjectToImageFilter_h

#include "itkImageSource.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkConceptChecking.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkResampleImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"

namespace itk
{

/** \class MaskedSpatialObjectToImageFilter
 * \brief Base class for filters that take a SpatialObject
 *        as input and produce an image as output, using a Mask to
 *        help decrease computation time.
 *  The user can set a pre-defined mask using SetMask(), else this
 *  filter will automatically compute a mask at a lower resolution.
 *
 *  See "Using a Mask to Decrease Computation Time for SpatialObject
 *  to Image Conversions", Insight Journal, 2006 July-December.
 */
template <class TInputSpatialObject, class TOutputImage>
class ITK_EXPORT MaskedSpatialObjectToImageFilter :
  public SpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef MaskedSpatialObjectToImageFilter Self;
  typedef SpatialObjectToImageFilter<TInputSpatialObject, TOutputImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef typename TOutputImage::SizeType   SizeType;
  typedef typename TOutputImage::PointType   PointType;
  typedef TOutputImage  OutputImageType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputImageType::Pointer OutputImagePointer;
  typedef typename OutputImageType::ValueType  ValueType;
  typedef typename OutputImageType::SpacingType SpacingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MaskedSpatialObjectToImageFilter, SpatialObjectToImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Some typedefs for input SpatialObject. */
  typedef TInputSpatialObject InputSpatialObjectType;
  typedef typename InputSpatialObjectType::Pointer        InputSpatialObjectPointer;
  typedef typename InputSpatialObjectType::ConstPointer   InputSpatialObjectConstPointer;
  typedef typename TInputSpatialObject::ChildrenListType  ChildrenListType;

  /** Some convenient typedefs. */
  typedef double                                          MaskCoordRep;
  typedef unsigned char                                   MaskPixelType;
  typedef itk::Image<MaskPixelType, itkGetStaticConstMacro(OutputImageDimension)> MaskImageType;
  typedef typename MaskImageType::Pointer                 MaskImagePointer;
  typedef typename MaskImageType::ConstPointer            MaskImageConstPointer;

  /** Set/Get the image input of this process object.  */
  virtual void SetInput( const InputSpatialObjectType *object);
  const InputSpatialObjectType *GetInput(void);

  /** Set/Get the input mask image to help speed up processing.
    * If a mask is not specified, one will be automatically generated
    * using GetMaskResampleFactor(). */
  virtual void SetMask (const MaskImageType *mask);
  const MaskImageType *GetMask(void);

  /** If the Mask image is not set, this filter will automatically
   *  generate a mask by creating low-resolution version of the SpatialObject.
   *  This parameter sets the resampling factor for all dimensions.
   *  Example: 4.0 means that the mask is generated at 1/4 the output size.
   *  A large the factor means faster mask generation and less accuracy.
   *  A smaller factor means slower mask generation and more accuarcy.
   *  Default = 4.0 */
  itkSetMacro(MaskResampleFactor, double);
  itkGetMacro(MaskResampleFactor, double);

  /** Sets/Gets the size of the binary dilation structuring element used to
   *  generate the Mask.
   *  The larger the value, the more accurate the mask but the slower the generation.
   *  Default = 2. */
  itkSetMacro(MaskDilationSize, unsigned int);
  itkGetMacro(MaskDilationSize, unsigned int);

protected:
  MaskedSpatialObjectToImageFilter();
  ~MaskedSpatialObjectToImageFilter();

  virtual void GenerateData();
  virtual void GenerateDataUsingMask();
  virtual void SetupOutputImage();
  virtual void GenerateMask();

  bool         m_UserSuppliedMask;
  double       m_MaskResampleFactor;
  unsigned int m_MaskDilationSize;

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

private:
  MaskedSpatialObjectToImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMaskedSpatialObjectToImageFilter.hxx"
#endif

#endif
