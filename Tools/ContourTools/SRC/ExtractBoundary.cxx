#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkCommand.h>

#include <itkRegionOfInterestImageFilter.h>
#include <itkExtractImageFilter.h>

// Operations on Images

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <vnl/vnl_math.h>

// For Feature
#include <itkRescaleIntensityImageFilter.h>


#include <itkMinimumMaximumImageCalculator.h>
#include <itkStatisticsImageFilter.h>

#include <itkScalarImageToTextureFeaturesFilter.h>
#include <itkScalarImageToRunLengthFeaturesFilter.h>

#include <itkLabelGeometryImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>

// For Label Object Representation
#include <itkLabelObject.h>
#include <itkLabelMap.h>
#include <itkShapeLabelObject.h>
#include <itkStatisticsLabelObject.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkShapeLabelMapFilter.h>
#include <itkStatisticsLabelMapFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkBinaryImageToStatisticsLabelMapFilter.h>

// For island removing
#include <itkBinaryShapeKeepNObjectsImageFilter.h>

#include <itkBinaryFillholeImageFilter.h>
#include <itkSliceBySliceImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkAndImageFilter.h>

#include "itkFlatStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"

// For file operation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// For threshold
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include "ITKUtils.h"


// Definition of Feature Extraction
typedef itk::MinimumMaximumImageCalculator< InputImageType > MaxMinFilterType;

using namespace std;

int main( int argc, char *argv[] )
{


    if ( argc < 3)
    {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage = " << argv[0];
        cerr << " InputLabelImage OutputLabelImage ThicknessOfBoundary={2} Mode={2D}" << endl;
        return EXIT_FAILURE;
    }


    double radius = 2;
    if (argc > 3)
    {
        radius = atoi(argv[3]);
    }

    char mode3D = 0;
    if (argc > 4)
    {
        cout << argv[4] << endl;
        if(strcmp(argv[4],"3D")==0)
            mode3D = 1;
        else if(strcmp(argv[4],"3D1")==0)
            mode3D = 2;
    }

    bool isTumorOnly = 0;
    bool isExpand = 0;
    if (argc > 5)
    {
        cout << argv[5] << endl;
        if(strcmp(argv[5],"T")==0)
            isTumorOnly = 1;
        if(strcmp(argv[5],"TE")==0)
        {
            isTumorOnly = 1;
            isExpand = 1;
        }
    }


    string inputLabelImageName = argv[1];
    string outputLabelImageName = argv[2];
    

    cout << "Input Label Image Name = " << inputLabelImageName << endl;
    cout << "Output Label Image Name = " << outputLabelImageName << endl;
    

    ///////// To Read the Mask image ////////////////////
    cout << "Get the mask from Input" << endl;
    LabelImageType::Pointer inputLabelImage = ReadImageFile<LabelImageType>(inputLabelImageName);
    
    // To check coordinate of the input label image
    SpacingType inputImageSpacing = inputLabelImage->GetSpacing();
    OriginType  inputImageOrigin  = inputLabelImage->GetOrigin();
    RegionType  inputImageRegion  = inputLabelImage->GetLargestPossibleRegion();
    SizeType    inputImageSize    = inputImageRegion.GetSize();


    /////////////////////////////

    cout << "Input Image Spacing = " << inputImageSpacing << endl;
    cout << "Input Image Origin = " << inputImageOrigin << endl;
    cout << "Input Image Size = " << inputImageSize << endl << endl << endl;


    // Image ROI
    MaskImageType::RegionType maskRegion;

    MaskImageType::Pointer normalMaskImage;
    MaskImageType::Pointer tumorMaskImage;

    normalMaskImage = GetMaskImage(inputLabelImage, 306); // normal
    tumorMaskImage = GetMaskImage(inputLabelImage, 307); // tumor
    maskRegion = GetRoi(tumorMaskImage);
    ExpandRoi(tumorMaskImage, maskRegion);
    BoundingCheck(normalMaskImage, inputLabelImage, maskRegion, inputImageRegion);
    BoundingCheck(tumorMaskImage, inputLabelImage, maskRegion, inputImageRegion);

    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;

    RegionToIndex(maskRegion, roiStart, roiEnd);

    MaskImageType::Pointer tumorMaskImageROI = ApplyRoi<MaskImageType>(tumorMaskImage, maskRegion);
    MaskImageType::Pointer normalMaskImageROI = ApplyRoi<MaskImageType>(normalMaskImage, maskRegion);


    //Definition of Morphological operation sturcture elements
    typedef itk::FlatStructuringElement<Dimension> StructuringElementType;
    typedef itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType> DilateFilterType;
    typedef itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType> ErodeFilterType;

    typedef itk::AndImageFilter <MaskImageType> AndImageFilterType;
    typedef itk::SubtractImageFilter <MaskImageType, MaskImageType> SubtractImageFilterType;

    typedef itk::Image<MaskPixelType, Dimension-1> MaskImage2DType;
    typedef itk::FlatStructuringElement<Dimension-1> StructuringElement2DType;
    typedef itk::BinaryDilateImageFilter<MaskImage2DType, MaskImage2DType, StructuringElement2DType> DilateFilter2DType;
    typedef itk::BinaryErodeImageFilter<MaskImage2DType, MaskImage2DType, StructuringElement2DType> ErodeFilter2DType;

    typedef itk::SliceBySliceImageFilter<MaskImageType, MaskImageType> SliceBySliceFilterType;

    StructuringElementType::RadiusType elementRadius;
    elementRadius.Fill(static_cast<unsigned int>(round(radius/inputImageSpacing[0])));
    elementRadius[2] = static_cast<unsigned int>(round(radius/inputImageSpacing[2]));
    if(mode3D == 2)
        elementRadius[2] = 1;

    StructuringElement2DType::RadiusType elementRadius2D;
    elementRadius2D.Fill(static_cast<unsigned int>(round(radius/inputImageSpacing[0])));

    StructuringElementType structuringElement = StructuringElementType::Ball(elementRadius, true);
    StructuringElement2DType structuringElement2D = StructuringElement2DType::Ball(elementRadius2D, true);
 

    MaskImageType::Pointer outputNormalMaskImageROI;
    MaskImageType::Pointer outputTumorMaskImageROI;
    {
        // Nb = N&dilate(T,size), 306
        AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
        itk::ImageToImageFilter<MaskImageType, MaskImageType>::Pointer filterPointer;
        if(mode3D==0)
        {
            cout << elementRadius2D << endl;
            SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

            DilateFilter2DType::Pointer dilate2DFilter = DilateFilter2DType::New();
            dilate2DFilter->SetKernel(structuringElement2D);
            dilate2DFilter->SetDilateValue(maskValue);

            sliceBySliceFilter->SetInput(tumorMaskImageROI);
            sliceBySliceFilter->SetFilter(dilate2DFilter);

            filterPointer = sliceBySliceFilter;
        }
        else
        {
            cout << elementRadius << endl;
            DilateFilterType::Pointer dilateFilter = DilateFilterType::New();

            dilateFilter->SetInput(tumorMaskImageROI);
            dilateFilter->SetKernel(structuringElement);
            dilateFilter->SetDilateValue(maskValue);

            filterPointer = dilateFilter;
        }
        ShowProgressObject progressWatch(filterPointer);
        itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
        command = itk::SimpleMemberCommand<ShowProgressObject>::New();
        command->SetCallbackFunction(&progressWatch,
                                     &ShowProgressObject::ShowProgress);

        filterPointer->AddObserver(itk::ProgressEvent(), command);


        filterPointer->Update();

        if(isTumorOnly)
        {
            outputNormalMaskImageROI = MaskImageType::New();
            outputNormalMaskImageROI->CopyInformation( normalMaskImageROI );
            outputNormalMaskImageROI->SetSpacing(inputImageSpacing);
            outputNormalMaskImageROI->SetBufferedRegion( normalMaskImageROI->GetBufferedRegion() );
            outputNormalMaskImageROI->Allocate();
            outputNormalMaskImageROI->FillBuffer(0);
            outputTumorMaskImageROI = filterPointer->GetOutput();
        }
        else
        {
            andFilter->SetInput1(normalMaskImageROI);
            andFilter->SetInput2(filterPointer->GetOutput());
            andFilter->Update();

            outputNormalMaskImageROI = andFilter->GetOutput();
        }
        cout << endl;
    }
 
    if(isTumorOnly == 0)
    {
        // Tb = T & dilate(N,size), 307
        AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
        itk::ImageToImageFilter<MaskImageType, MaskImageType>::Pointer filterPointer;
        if(mode3D==0)
        {
            cout << elementRadius2D << endl;
            SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

            DilateFilter2DType::Pointer dilate2DFilter = DilateFilter2DType::New();
            dilate2DFilter->SetKernel(structuringElement2D);
            dilate2DFilter->SetDilateValue(maskValue);

            sliceBySliceFilter->SetInput(normalMaskImageROI);
            sliceBySliceFilter->SetFilter(dilate2DFilter);

            filterPointer = sliceBySliceFilter;
        }
        else
        {
            DilateFilterType::Pointer dilateFilter = DilateFilterType::New();

            dilateFilter->SetInput(normalMaskImageROI);
            dilateFilter->SetKernel(structuringElement);
            dilateFilter->SetDilateValue(maskValue);

            filterPointer = dilateFilter;
        }
        ShowProgressObject progressWatch(filterPointer);
        itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
        command = itk::SimpleMemberCommand<ShowProgressObject>::New();
        command->SetCallbackFunction(&progressWatch,
                                     &ShowProgressObject::ShowProgress);

        filterPointer->AddObserver(itk::ProgressEvent(), command);


        filterPointer->Update();

        if(isTumorOnly==0) // intersection
        {
            andFilter->SetInput1(tumorMaskImageROI);
            andFilter->SetInput2(filterPointer->GetOutput());
            andFilter->Update();
            outputTumorMaskImageROI = andFilter->GetOutput();
        }

        cout << endl;
    }

    cout << "Complete Extraction!" << endl;

    LabelImageType::Pointer outputLabelImage = LabelImageType::New();
    // allocate outputLabelImage first
    outputLabelImage->CopyInformation( inputLabelImage );
    outputLabelImage->SetSpacing(inputImageSpacing);
    outputLabelImage->SetBufferedRegion( inputLabelImage->GetBufferedRegion() );
    outputLabelImage->Allocate();
    outputLabelImage->FillBuffer(0);


    if(isTumorOnly == 0)
    {
        itk::ImageRegionIterator< MaskImageType > normalOut(outputNormalMaskImageROI, outputNormalMaskImageROI->GetBufferedRegion());
        itk::ImageRegionIterator< MaskImageType > tumorOut(outputTumorMaskImageROI, outputTumorMaskImageROI->GetBufferedRegion());
        itk::ImageRegionIterator< LabelImageType > out(outputLabelImage, inputImageRegion);

        for (normalOut.GoToBegin(), tumorOut.GoToBegin(), out.GoToBegin(); !normalOut.IsAtEnd(); ++normalOut, ++tumorOut, ++out)
        {
            if(normalOut.Get() > 0 )
            {
                //cout << "N ";
                out.Set(306);
            }
            if(tumorOut.Get() > 0 )
            {
                //cout << "T ";
                out.Set(307);
            }
        }
    }
    else
    {
        itk::ImageRegionIterator< MaskImageType > tumor(tumorMaskImageROI, tumorMaskImageROI->GetBufferedRegion());
        itk::ImageRegionIterator< MaskImageType > tumorOut(outputTumorMaskImageROI, outputTumorMaskImageROI->GetBufferedRegion());
        itk::ImageRegionIterator< LabelImageType > out(outputLabelImage, inputImageRegion);

        for (tumor.GoToBegin(), tumorOut.GoToBegin(), out.GoToBegin(); !tumor.IsAtEnd(); ++tumor, ++tumorOut, ++out)
        {
            if((isExpand == 0 && tumor.Get() ==0 && tumorOut.Get() > 0)
                || (isExpand==1 && tumorOut.Get() > 0))
            {
                //cout << "T ";
                out.Set(307);
            }
        }
    }

        


    std::cout << "writeImage_spacing = " << outputLabelImage->GetSpacing() << std::endl ;
    std::cout << "writeImage_origin  = "  << outputLabelImage->GetOrigin() << std::endl ;
    std::cout << "writeImage_LargestPossibleRegion = " << outputLabelImage->GetLargestPossibleRegion() << std::endl ;


    ////////To write Output images /////////////////////
    try
    {
        WriteImageFile<LabelImageType>(outputLabelImage,outputLabelImageName);
    }
    catch (itk::ExceptionObject &excp)
    {
        std::cerr << "Exception thrown while writing the series " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }


    std::cout << "Saved the normal tissue and tumor boundary image" << std::endl;

    return EXIT_SUCCESS;
}
