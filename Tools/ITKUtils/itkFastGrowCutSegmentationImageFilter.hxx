#ifndef _itkFastGrowCutSegmentationImageFilter_txx_
#define _itkFastGrowCutSegmentationImageFilter_txx_


#include "itkFastGrowCutSegmentationImageFilter.h"


#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkImageFileWriter.h"
#include "itkLabelMap.h"
#include "itkShapeLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

#include "itkProgressReporter.h"
# include "itkIterationReporter.h"

#include <vcl_map.h>
#include <vcl_string.h>
#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_utility.h>

#include <iostream>
#include <ctime>

#include "itkTimeProbe.h"

#include "FastGrowCutSegmenter.h"


// #define LOG_OUTPUT
// #define DEBUG_LEVEL_0

#define USE_LABELSHAPEFILTER

namespace itk
{

template< class TInputImage, class TOutputImage>
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::FastGrowCutSegmentationImageFilter()
{
    this->SetNumberOfRequiredInputs(2);

    m_LabelImage = OutputImageType::New();

    m_fastGC = NULL;
}


/**
 * Standard PrintSelf method.
 */
template <class TInputImage, class TOutputImage>
void
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
}

template <class TInputImage, class TOutputImage>
void
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
    Superclass::GenerateInputRequestedRegion();
    if ( this->GetInput() )
    {
        InputImagePointer image = const_cast< InputImageType * >( this->GetInput() );
        image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template <class TInputImage, class TOutputImage>
void FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
    Superclass::EnlargeOutputRequestedRegion(output);
    output->SetRequestedRegionToLargestPossibleRegion();
}

template< class TInputImage, class TOutputImage>
const typename FastGrowCutSegmentationImageFilter< TInputImage, TOutputImage>::OutputImagePointer
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::GetLabelImage()
{

    //typename OutputImageType::Pointer labelImage = OutputImageType::New();
    //labelImage->Graft(this->ProcessObject::GetInput(1));
    //return labelImage;
    return m_LabelImage;
    //return const_cast< const OutputImageType*> (this->ProcessObject::GetInput(1));
}

template< class TInputImage, class TOutputImage>
const typename FastGrowCutSegmentationImageFilter< TInputImage, TOutputImage>::SeedImagePointer
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::GetSeedImage()
{
    return const_cast< SeedImageType *> (this->GetInput(1));
}

template< class TInputImage, class TOutputImage>
const typename FastGrowCutSegmentationImageFilter< TInputImage, TOutputImage>::InputImagePointer
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::GetInputImage()
{
    return const_cast< InputImageType *> (this->GetInput(0));
}

template <class TInputImage, class TOutputImage>
void
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::Initialize( OutputImageType *output)
{

    // is the output allocated already
    if (output->GetBufferedRegion().GetNumberOfPixels() == 0)
    {
        // allocate memory for the output buffer
        output->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
        output->Allocate();
    }


    if (m_LabelImage->GetBufferedRegion().GetNumberOfPixels() == 0)
    {
        m_LabelImage->CopyInformation( output );
        m_LabelImage->SetBufferedRegion( output->GetBufferedRegion() );
        m_LabelImage->Allocate();
    }
}

template<class TInputImage, class TOutputImage>
void FastGrowCutSegmentationImageFilter< TInputImage, TOutputImage>::
ComputeLabelVolumes(TOutputImage *outputImage, vcl_vector< unsigned > &volumes,
                    vcl_vector< unsigned > &physicalVolumes)
{
    //volumes.resize(3);

#ifndef USE_LABELSHAPEFILTER

    vcl_map<unsigned short, unsigned int>labelMap;
    vcl_vector< unsigned short > labels;
    unsigned int index = 0;
    ImageRegionConstIteratorWithIndex< OutputImageType > label( outputImage,
            outputImage->GetBufferedRegion());

    for (label.GoToBegin(); !label.IsAtEnd(); ++label)
    {

        unsigned short pix = static_cast<unsigned short>(label.Get());
        vcl_map<unsigned short, unsigned int>::iterator it = labelMap.find(pix);
        if (it == labelMap.end())
        {
            labelMap.insert(vcl_pair<unsigned short, unsigned int>(pix, index));
            volumes.push_back(1);
            labels.push_back(pix);
            ++index;
        }
        else
        {
            int i = it->second;
            ++volumes[i];
        }
    }
    // std::cout<<" label \t "<<" volume "<<std::endl;
    // for (unsigned int i = 0; i < index; i++)
    //   {
    //   vcl_map<unsigned short, unsigned int>::iterator it = labelMap.find(labels[i]);
    //   std::cout<<labels[i]<<"\t"<<volumes[it->second]<<std::endl;
    //   }

#else
    using LabelFilterType = LabelImageToShapeLabelMapFilter< OutputImageType >;
    typename LabelFilterType::Pointer labelFilter = LabelFilterType::New();

    labelFilter->SetInput(outputImage);
    labelFilter->SetBackgroundValue(0);
    labelFilter->Update();

    unsigned long numObjects = labelFilter->GetOutput()->GetNumberOfLabelObjects();
    volumes.resize((const int)numObjects);
    physicalVolumes.resize((const int)numObjects);
    for (unsigned n = 0; n < numObjects; n++)
    {
#if ITK_VERSION_MAJOR < 4
        volumes[n] = labelFilter->GetOutput()->
                     GetNthLabelObject(n)->GetSize();
        physicalVolumes[n] = (unsigned int)
                             (labelFilter->GetOutput()->
                              GetNthLabelObject(n)->GetSize() / 1000.0);
#else
        volumes[n] = labelFilter->GetOutput()->
                     GetNthLabelObject(n)->GetNumberOfPixels();
        physicalVolumes[n] = (unsigned int)
                             (labelFilter->GetOutput()->
                              GetNthLabelObject(n)->GetNumberOfPixels() / 1000.0);
#endif
    }
#endif
}




template <class TInputImage, class TOutputImage>
void
FastGrowCutSegmentationImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
    IterationReporter iterate(this, 0, 1);

    itk::TimeProbe timer;

    m_InitializationFlag = false;
    if (m_fastGC == NULL)
    {
        m_fastGC = new FGC::FastGrowCut<InputPixelType, OutputPixelType>();
    }

    typename OutputImageType::Pointer output = this->GetOutput();
    this->Initialize(output);


    unsigned int ndims = static_cast< unsigned int>(output->GetImageDimension());

    typename InputImageType::Pointer inputImage = InputImageType::New();
    inputImage->Graft( this->ProcessObject::GetInput(0));

    typename SeedImageType::Pointer seedImage = SeedImageType::New();
    seedImage->Graft( this->ProcessObject::GetInput(1));

    timer.Start();

    // Find ROI
    if (!m_InitializationFlag)
    {
        FGC::FindITKImageROI<InputImageType>(inputImage, m_imROI);
        std::cout << "image ROI = [" << m_imROI[0] << "," << m_imROI[1] << "," << m_imROI[2] << ";"  \
                  << m_imROI[3] << "," << m_imROI[4] << "," << m_imROI[5] << "]" << std::endl;
        FGC::ExtractITKImageROI<InputPixelType>(inputImage, m_imROI, m_imSrcVec);
    }

    FGC::ExtractITKImageROI<OutputPixelType>(seedImage, m_imROI, m_imSeedVec);

    // Initialize FastGrowCut
    std::vector<long> imSize(3);
    for (int i = 0; i < 3; i++)
    {
        imSize[i] = m_imROI[i + 3] - m_imROI[i];
    }

    m_fastGC->SetSourceImage(m_imSrcVec);
    m_fastGC->SetSeedlImage(m_imSeedVec);
    m_fastGC->SetImageSize(imSize);
    m_fastGC->SetWorkMode(m_InitializationFlag);

    this->UpdateProgress(0.1);

    // Do Segmentation
    m_fastGC->DoSegmentation();
    //m_fastGC->GetForegroundmage(m_imLabVec);
    m_fastGC->GetLabeImage(m_imLabVec);

    this->UpdateProgress(0.9);

    // Update result
    FGC::UpdateITKImageROI<OutputPixelType>(m_imLabVec, m_imROI, m_LabelImage);

    timer.Stop();

    if (!m_InitializationFlag)
        std::cout << "Initial fast GrowCut segmentation time: " << timer.GetMean() << " seconds\n";
    else
        std::cout << "Adaptive fast GrowCut segmentation time: " << timer.GetMean() << " seconds\n";

    this->UpdateProgress(1.0);

    this->GraftOutput(m_LabelImage);

    vcl_vector< unsigned > labelVolumes;
    vcl_vector< unsigned > physicalVolumes;

    this->ComputeLabelVolumes(m_LabelImage, labelVolumes, physicalVolumes);

    // std::cout<<"...................................................."<<std::endl;
    // for (unsigned n = 0; n < labelVolumes.size(); n++)
    //   {
    //   std::cout<<" Object "<<n<<" volume "<<labelVolumes[n]<<" Physical Volume "<<
    //     physicalVolumes[n]<<std::endl;
    //   }

    /*
    using WriterType = ImageFileWriter < OutputImageType >;
    typename WriterType::Pointer writer = WriterType::New();

    vcl_string outfilename = "segmented_img.mhd";
    writer->SetFileName(outfilename.c_str());
    writer->SetInput(m_LabelImage);
    writer->Update();
    */
    // std::cout<<" Time taken for segmentation "<<iter<<
    //   " : "<<difftime(end, start)<<" secs"<<std::endl;

    if (m_fastGC != NULL)
    {
        delete m_fastGC;
    }

}

}//namespace itk

#endif
