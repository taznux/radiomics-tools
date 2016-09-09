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
        cerr << " Output InputLabelImage1 InputLabelImage2 InputLabelImage3 ..." << endl;
        return EXIT_FAILURE;
    }


    string outputName = argv[1];
    cout << "Output Name = " << outputName << endl;

    ///////// To Read the Mask images ////////////////////
    int numInputs = argc - 2;
    vector<string> inputLabelImageNames;
    vector<LabelImageType::Pointer> inputLabelImages;
    cout << "The number of input imeages = " << numInputs << endl;
    for (int i = 0; i < numInputs; i++)
    {
        inputLabelImageNames.push_back(argv[i + 2]);
        cout << "Input Label Image " << i << " Name = " << inputLabelImageNames[i] << endl;
        inputLabelImages.push_back(ReadImageFile<LabelImageType>(inputLabelImageNames[i]));
        cout << "Get the mask from Input " << i << endl;
    }

    // To check coordinate of the input label image
    SpacingType inputImageSpacing = inputLabelImages[0]->GetSpacing();
    OriginType  inputImageOrigin  = inputLabelImages[0]->GetOrigin();
    RegionType  inputImageRegion  = inputLabelImages[0]->GetLargestPossibleRegion();
    SizeType    inputImageSize    = inputImageRegion.GetSize();


    /////////////////////////////

    cout << "Input Image Spacing = " << inputImageSpacing << endl;
    cout << "Input Image Origin = " << inputImageOrigin << endl;
    cout << "Input Image Size = " << inputImageSize << endl << endl << endl;


    // Image ROI
    MaskImageType::RegionType roiRegion;
    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;

    roiStart.Fill(1000);
    roiEnd.Fill(0);

    for (int i = 0; i < numInputs; i++)
    {
        MaskImageType::RegionType roiRegion1;
        MaskImageType::IndexType roiStart1;
        MaskImageType::IndexType roiEnd1;
    
        roiRegion1 = GetRoi(inputLabelImages[i]);
        ExpandRoi(inputLabelImages[i], roiRegion1, 10);
        RegionToIndex(roiRegion1, roiStart1, roiEnd1);

        for (unsigned i = 0; i < Dimension; i++)
        {
            if (roiStart[i] > roiStart1[i]) roiStart[i] = roiStart1[i];
            if (roiEnd[i] < roiEnd1[i]) roiEnd[i] = roiEnd1[i];
        }
    }

    roiRegion = RoiIndexToRegion(roiStart, roiEnd);

    cout << endl << "ROI start " << roiStart << endl;
    cout << "ROI end " << roiEnd << endl;
    cout << "ROI size " << roiRegion.GetSize() << endl;

    WriteROI(roiRegion,outputName);

    return EXIT_SUCCESS;
}
