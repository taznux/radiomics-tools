/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkMaskedSpatialObjectToImageFilter.cxx,v $
Language:  C++
Date:      $Date: 2006/07/19 20:06:17 $
Version:   $Revision: 1.00 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskedSpatialObjectToImageFilter_cxx
#define __itkMaskedSpatialObjectToImageFilter_cxx

#include "itkMaskedSpatialObjectToImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
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

/** Constructor */
template <class TInputSpatialObject, class TOutputImage>
MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::MaskedSpatialObjectToImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->m_ChildrenDepth = 1;
  this->m_Size.Fill(0);

  for (unsigned int i = 0; i < itkGetStaticConstMacro(OutputImageDimension); i++)
  {
    this->m_Spacing[i] = 1.0;
    this->m_Origin[i] = 0.0;
  }

  this->m_InsideValue = 0;
  this->m_OutsideValue = 0;
  this->m_UseObjectValue = false;

  this->m_UserSuppliedMask = false;
  this->m_MaskResampleFactor = 4.0;
  this->m_MaskDilationSize = 2;
}

/** Destructor */
template <class TInputSpatialObject, class TOutputImage>
MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::~MaskedSpatialObjectToImageFilter()
{
}

/** Set the Input SpatialObject */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::SetInput(const InputSpatialObjectType *input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0,
                                   const_cast<InputSpatialObjectType *>(input));
}

/** Set the Input Mask image */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::SetMask(const MaskImageType *mask)
{
  if (mask)
  {
    this->m_UserSuppliedMask = true;

    // Const-correct the mask image
    this->ProcessObject::SetNthInput(1,
                                     const_cast<MaskImageType *>(mask));
  }
}

/** Get the input Spatial Object */
template <class TInputSpatialObject, class TOutputImage>
const typename MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::InputSpatialObjectType *
MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
  {
    return 0;
  }

  return static_cast<const TInputSpatialObject *>(this->ProcessObject::GetInput(0));
}

/** Get the input Mask image */
template <class TInputSpatialObject, class TOutputImage>
const typename MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::MaskImageType *
MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::GetMask(void)
{
  if (this->GetNumberOfInputs() < 2)
  {
    return 0;
  }

  return static_cast<const MaskImageType *>(this->ProcessObject::GetInput(1));
}

/** GenerateData */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::GenerateData(void)
{
  this->SetupOutputImage();

  if (!this->m_UserSuppliedMask)
    //We must generate the mask ourselves
    this->GenerateMask();

  this->GenerateDataUsingMask();
}

/** SetupOutputImage */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::SetupOutputImage(void)
{
  unsigned int i;

  // Get the input, and output pointers
  const InputSpatialObjectType *InputObject = this->GetInput();
  OutputImagePointer OutputImage = this->GetOutput();

  // Generate the image
  SizeType size;

#if ITK_VERSION_MAJOR == 5
  InputObject->ComputeFamilyBoundingBox(this->m_ChildrenDepth);
  for (i = 0; i < this->ObjectDimension; i++)
  {
    size[i] = static_cast<SizeValueType>(
        InputObject->GetFamilyBoundingBoxInWorldSpace()->GetMaximum()[i] - InputObject->GetFamilyBoundingBoxInWorldSpace()->GetMinimum()[i]);
  }
#else
  InputObject->ComputeBoundingBox();
  for (unsigned int i = 0; i < itkGetStaticConstMacro(ObjectDimension); i++)
  {
    size[i] = (long unsigned int)(InputObject->GetBoundingBox()->GetMaximum()[i] - InputObject->GetBoundingBox()->GetMinimum()[i]);
  }
#endif

  typename OutputImageType::IndexType index;
  index.Fill(0);
  typename OutputImageType::RegionType region;

  // If the size of the output has been explicitly specified, the filter
  // will set the output size to the explicit size, otherwise the size from the
  // spatial
  // object's bounding box will be used as default.

  bool specified = false;
  for (i = 0; i < this->OutputImageDimension; i++)
  {
    if (this->m_Size[i] != 0)
    {
      specified = true;
      break;
    }
  }

  if (specified)
  {
    region.SetSize(this->m_Size);
  }
  else
  {
    region.SetSize(size);
  }
  region.SetIndex(index);

  OutputImage->SetLargestPossibleRegion(region); //
  OutputImage->SetBufferedRegion(region);        // set the region
  OutputImage->SetRequestedRegion(region);       //
  OutputImage->SetSpacing(this->m_Spacing);      // set spacing
  OutputImage->SetOrigin(this->m_Origin);        //   and origin
  OutputImage->SetDirection(this->m_Direction);
  OutputImage->Allocate(); // allocate the image
}

/** Generate Mask */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::GenerateMask(void)
{
  // NOTE: The mask is store in Inputs(1)
  MaskImagePointer maskImage = nullptr;

  //Convert SpatialObject to low-resolution mask
  //NOTE: We set the InsideValue and Outside values to max() and zero, and
  //      tell the filter to not use the object value to ensure that
  //      grey-scale SpatialObjects are converted to a binary mask.
  using SpatialObjectToImageFilterType = itk::SpatialObjectToImageFilter<TInputSpatialObject, MaskImageType>;
  typename SpatialObjectToImageFilterType::Pointer filterConvertMask =
      SpatialObjectToImageFilterType::New();
  filterConvertMask->SetInput(this->GetInput());
  typename SpatialObjectToImageFilterType::SizeType maskSize;       //Set Size
  typename SpatialObjectToImageFilterType::SpacingType maskSpacing; //Set Spacing
  for (unsigned int dim = 0; dim < itkGetStaticConstMacro(OutputImageDimension); dim++)
  {
    maskSize[dim] = (unsigned int)(this->m_Size[dim] / this->m_MaskResampleFactor);
    maskSpacing[dim] = this->m_Spacing[dim] * ((double)this->m_Size[dim] / (double)maskSize[dim]);
  }
  //std::cout << this->m_Size << " " << maskSize << std::endl;
  //std::cout << (SpacingType)this->m_Spacing << " " << maskSpacing << std::endl;
  filterConvertMask->SetOrigin(this->m_Origin);
  filterConvertMask->SetSize(maskSize);
  filterConvertMask->SetSpacing(maskSpacing);
  filterConvertMask->SetInsideValue(itk::NumericTraits<MaskPixelType>::max());
  filterConvertMask->SetOutsideValue(itk::NumericTraits<MaskPixelType>::Zero);
  filterConvertMask->SetUseObjectValue(false);
  filterConvertMask->Update();

  //Create structuring element
  using BinaryBallElementType = itk::BinaryBallStructuringElement<MaskPixelType, itkGetStaticConstMacro(OutputImageDimension)>;
  BinaryBallElementType kernel;
  kernel.SetRadius(this->m_MaskDilationSize);
  kernel.CreateStructuringElement();

  //Binary dilate the resampled Mask image
  using BinaryDilateImageFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, BinaryBallElementType>;
  typename BinaryDilateImageFilterType::Pointer filterDilateMask = BinaryDilateImageFilterType::New();
  filterDilateMask->SetInput(filterConvertMask->GetOutput());
  filterDilateMask->SetKernel(kernel);
  filterDilateMask->SetDilateValue(itk::NumericTraits<MaskPixelType>::max());
  filterDilateMask->Update();

  //Resample mask to full size
  using ResampleFilterType = itk::ResampleImageFilter<MaskImageType, MaskImageType, MaskCoordRep>;
  using TransformType = itk::CenteredAffineTransform<MaskCoordRep, itkGetStaticConstMacro(OutputImageDimension)>;
  using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<MaskImageType, MaskCoordRep>;

  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  typename ResampleFilterType::Pointer filterMaskResample = ResampleFilterType::New();
  filterMaskResample->SetTransform(transform);
  filterMaskResample->SetInterpolator(interpolator);
  filterMaskResample->SetInput(filterDilateMask->GetOutput());
  double finalMaskOrigin[itkGetStaticConstMacro(OutputImageDimension)]; //Set Origin
  typename ResampleFilterType::SizeType finalMaskSize;                  //Set Size
  typename ResampleFilterType::SpacingType finalMaskSpacing;            //Set Spacing
  for (unsigned int dim = 0; dim < itkGetStaticConstMacro(OutputImageDimension); dim++)
  {
    finalMaskOrigin[dim] = this->m_Origin[dim];
    finalMaskSize[dim] = this->m_Size[dim];
    finalMaskSpacing[dim] = this->m_Spacing[dim];
  }
  //std::cout << filterDilateMask->GetOutput()->GetLargestPossibleRegion().GetSize() << " " << finalMaskSize << std::endl;
  //std::cout << filterDilateMask->GetOutput()->GetSpacing() << " " << finalMaskSpacing << std::endl;

  filterMaskResample->SetOutputOrigin(finalMaskOrigin);
  filterMaskResample->SetSize(finalMaskSize);
  filterMaskResample->SetOutputSpacing(finalMaskSpacing);
  filterMaskResample->Update();

  //Set mask image
  maskImage = filterMaskResample->GetOutput();
  this->SetMask(maskImage);
}

/** GenerateData (using Mask) */
template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::GenerateDataUsingMask(void)
{
  itkDebugMacro(<< "MaskedSpatialObjectToImageFilter::Update() called");

  //Get the input, mask and output pointers
  const InputSpatialObjectType *InputObject = this->GetInput();
  OutputImagePointer OutputImage = this->GetOutput();
  const MaskImageType *MaskImage = this->GetMask();

  //Get output region
  typename OutputImageType::RegionType region =
      OutputImage->GetLargestPossibleRegion();

  //Setup Progress reporter
  ProgressReporter progress(this, 0, region.GetNumberOfPixels());

  using OutputIteratorType = itk::ImageRegionIteratorWithIndex<OutputImageType>;
  using MaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;

  OutputIteratorType itOutput(OutputImage, region);
  MaskIteratorType itMask(MaskImage, region);

  itk::Point<double, itkGetStaticConstMacro(ObjectDimension)> point;

  while (!itOutput.IsAtEnd())
  {
    //Check if the mask is "ON" for this point
    if (itMask.Get() > NumericTraits<MaskPixelType>::Zero)
    {
      // ValueAt requires the point to be in physical coordinate i.e
      for (unsigned int i = 0; i < itkGetStaticConstMacro(ObjectDimension); i++)
      {
        point[i] = (itOutput.GetIndex()[i] * this->m_Spacing[i]) + this->m_Origin[i];
      }
      double val = 0.0;
#if ITK_VERSION_MAJOR == 5
      InputObject->ValueAtInObjectSpace(point, val, TInputSpatialObject::MaximumDepth);
#else
      InputObject->ValueAt(point, val, TInputSpatialObject::MaximumDepth);
#endif

      if (this->m_InsideValue != 0 || this->m_OutsideValue != 0)
      {
        if (val)
        {
          if (this->m_UseObjectValue)
          {
            itOutput.Set(static_cast<ValueType>(val));
          }
          else
          {
            itOutput.Set(this->m_InsideValue);
          }
        }
        else
        {
          itOutput.Set(this->m_OutsideValue);
        }
      }
      else
      {
        itOutput.Set(static_cast<ValueType>(val));
      }
    }
    else
    {
      //The mask is "OFF" for this region - set as outside
      itOutput.Set(this->m_OutsideValue);
    }
    ++itOutput;
    ++itMask;
    progress.CompletedPixel();
  }

  itkDebugMacro(<< "MaskedSpatialObjectToImageFilter::Update() finished");

} // end update function

template <class TInputSpatialObject, class TOutputImage>
void MaskedSpatialObjectToImageFilter<TInputSpatialObject, TOutputImage>::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Size : " << this->m_Size << std::endl;
  os << indent << "Spacing : " << (SpacingType)this->m_Spacing << std::endl;
  os << indent << "Children depth : " << this->m_ChildrenDepth << std::endl;
  os << indent << "Inside Value : " << this->m_InsideValue << std::endl;
  os << indent << "Outside Value : " << this->m_OutsideValue << std::endl;
  if (this->m_UseObjectValue)
  {
    os << indent << "Using Object Value : ON" << std::endl;
  }
  else
  {
    os << indent << "Using Object Value : OFF" << std::endl;
  }
  if (this->m_UserSuppliedMask)
  {
    os << indent << "User Supplied Mask : YES" << std::endl;
  }
  else
  {
    os << indent << "User Supplied Mask : NO" << std::endl;
  }
}

} // end namespace itk

#endif
