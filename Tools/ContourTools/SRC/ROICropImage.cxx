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
        cerr << " OutputLabelImage inputImage [ROIfile or start(x,y,z) size(x,y,z)]" << endl;
        return EXIT_FAILURE;
    }


    string outputImageName = argv[1];
    string inputImageName = argv[2];
    string roiName = argv[3];
    
    cout << "Output Label Image Name = " << outputImageName << endl;
    cout << "Input Image 1 Name = " << inputImageName << endl;
    
    ///////// To Read the Reference image ////////////////////
    cout << "Get the Input Image" << endl;
    InputImageType::Pointer inputImage = ReadImageFile<InputImageType>(inputImageName);
    

    // To check coordinate of the input label image
    SpacingType inputImageSpacing = inputImage->GetSpacing();
    OriginType  inputImageOrigin  = inputImage->GetOrigin();
    RegionType  inputImageRegion  = inputImage->GetLargestPossibleRegion();
    SizeType    inputImageSize    = inputImageRegion.GetSize();


    /////////////////////////////

    cout << "Input Image Spacing = " << inputImageSpacing << endl;
    cout << "Input Image Origin = " << inputImageOrigin << endl;
    cout << "Input Image Size = " << inputImageSize << endl << endl << endl;


    // Read ROI
    InputImageType::RegionType roiRegion;
    InputImageType::IndexType roiStart;
    InputImageType::IndexType roiEnd;
    InputImageType::Pointer outputImage;

    roiRegion = ReadROI(roiName);
    cout << roiRegion << endl;
    RegionToIndex(roiRegion, roiStart, roiEnd);
    outputImage = ApplyRoi<InputImageType>(inputImage, roiRegion);

    cout << "writeImage_spacing = " << outputImage->GetSpacing() << endl ;
    cout << "writeImage_origin  = "  << outputImage->GetOrigin() << endl ;
    cout << "writeImage_LargestPossibleRegion = " << outputImage->GetLargestPossibleRegion() << endl ;


    ////////To write Output images /////////////////////
    WriteImageFile<InputImageType>(outputImage,outputImageName);


    std::cout << "Saved the normal tissue and tumor boundary image" << std::endl;

    return EXIT_SUCCESS;
}
