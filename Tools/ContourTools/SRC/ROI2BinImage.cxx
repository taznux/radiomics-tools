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

using namespace std;

int main( int argc, char *argv[] )
{
    if ( argc < 3)
    {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage = " << argv[0];
        cerr << " OutputLabelImage ReferenceImage [ROIfile or start(x,y,z) size(x,y,z)]" << endl;
        return EXIT_FAILURE;
    }


    string outputLabelImageName = argv[1];
    string refImageName = argv[2];
    string roiName = argv[3];
    
    cout << "Output Label Image Name = " << outputLabelImageName << endl;
    cout << "Reference Image 1 Name = " << refImageName << endl;
    
    ///////// To Read the Reference image ////////////////////
    cout << "Get the Reference Image" << endl;
    InputImageType::Pointer refImage = ReadImageFile<InputImageType>(refImageName);
    

    // To check coordinate of the input label image
    SpacingType refImageSpacing = refImage->GetSpacing();
    OriginType  refImageOrigin  = refImage->GetOrigin();
    RegionType  refImageRegion  = refImage->GetLargestPossibleRegion();
    SizeType    refImageSize    = refImageRegion.GetSize();


    /////////////////////////////

    cout << "Reference Image Spacing = " << refImageSpacing << endl;
    cout << "Reference Image Origin = " << refImageOrigin << endl;
    cout << "Reference Image Size = " << refImageSize << endl << endl << endl;


    // Read ROI
    MaskImageType::RegionType roiRegion;
    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;

    roiRegion = ReadROI(roiName);
    RegionToIndex(roiRegion, roiStart, roiEnd);

    

    LabelImageType::Pointer outputLabelImage = LabelImageType::New();
    // allocate outputLabelImage first
    outputLabelImage->CopyInformation( refImage );
    outputLabelImage->SetSpacing(refImageSpacing);
    outputLabelImage->SetBufferedRegion( refImageRegion );
    outputLabelImage->Allocate();
    outputLabelImage->FillBuffer(0);

    itk::ImageRegionIteratorWithIndex< LabelImageType > out(outputLabelImage, refImageRegion);

    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        bool isIn = true;
        LabelImageType::IndexType idx = out.GetIndex();
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (idx[i] < roiStart[i] || idx[i] > roiEnd[i]) isIn = false;
        }

        if (isIn) out.Set(1);
    }


    std::cout << "writeImage_spacing = " << outputLabelImage->GetSpacing() << std::endl ;
    std::cout << "writeImage_origin  = "  << outputLabelImage->GetOrigin() << std::endl ;
    std::cout << "writeImage_LargestPossibleRegion = " << outputLabelImage->GetLargestPossibleRegion() << std::endl ;


    ////////To write Output images /////////////////////
    WriteImageFile<LabelImageType>(outputLabelImage,outputLabelImageName);

    std::cout << "Saved the normal tissue and tumor boundary image" << std::endl;

    return EXIT_SUCCESS;
}
