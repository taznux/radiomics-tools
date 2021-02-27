#if defined(_MSC_VER)
# pragma warning ( disable : 4786 )
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
// #include "itkScalarImageToRunLengthFeaturesFilter1.h"
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

unsigned int NumberOfBins = 64;

void
decodeCharactheristicGLCM(int x, char name[])
{
    switch (x) {
        case 0:
            strcpy(name, "Energy");
            break;
        case 1:
            strcpy(name, "Entropy");
            break;
        case 2:
            strcpy(name, "Correlation");
            break;
        case 3:
            strcpy(name, "InverseDifferenceMoment");
            break;
        case 4:
            strcpy(name, "Inertia");
            break;
        case 5:
            strcpy(name, "ClusterShade");
            break;
        case 6:
            strcpy(name, "ClusterProminence");
            break;
        case 7:
            strcpy(name, "HaralickCorrelation");
            break;
    }
}

void
decodeCharactheristicGLRM(int x, char name[])
{
    switch (x) {
        case 0:
            strcpy(name, "ShortRunEmphasis");
            break;
        case 1:
            strcpy(name, "LongRunEmphasis");
            break;
        case 2:
            strcpy(name, "GreyLevelNonuniformity");
            break;
        case 3:
            strcpy(name, "RunLengthNonuniformity");
            break;
        case 4:
            strcpy(name, "LowGreyLevelRunEmphasis");
            break;
        case 5:
            strcpy(name, "HighGreyLevelRunEmphasis");
            break;
        case 6:
            strcpy(name, "ShortRunLowGreyLevelEmphasis");
            break;
        case 7:
            strcpy(name, "ShortRunHighGreyLevelEmphasis");
            break;
        case 8:
            strcpy(name, "LongRunLowGreyLevelEmphasis");
            break;
        case 9:
            strcpy(name, "LongRunHighGreyLevelEmphasis");
            break;
    }
} // decodeCharactheristicGLRM

void
CalcIntensityFeatures(InputImageType::Pointer inputImage, MaskImageType::Pointer maskImage, InputPixelType &inputMin,
  InputPixelType &inputMax, ostream &outfile)
{
    using LabelStatisticsFilterType = itk::LabelStatisticsImageFilter<InputImageType, MaskImageType>;

    LabelStatisticsFilterType::Pointer LabelStatisticsfilter = LabelStatisticsFilterType::New();
    LabelStatisticsfilter->UseHistogramsOn();
    LabelStatisticsfilter->SetHistogramParameters(3000, -1024, 2000);

    LabelStatisticsfilter->SetInput(inputImage);
    LabelStatisticsfilter->SetLabelInput(maskImage);

    LabelStatisticsfilter->Update();
    cout << "inputImage features from the LabelStatisticsImageFilter:" << endl;
    // outfile << endl<< "inputImage features from the LabelStatisticsImageFilter:" << endl;

    cout << "Minimum = " << LabelStatisticsfilter->GetMinimum(maskValue) << endl;
    cout << "Maximum = " << LabelStatisticsfilter->GetMaximum(maskValue) << endl;
    cout << "Median  = " << LabelStatisticsfilter->GetMedian(maskValue) << endl;
    cout << "Mean    = " << LabelStatisticsfilter->GetMean(maskValue) << endl;
    cout << "Sigma   = " << LabelStatisticsfilter->GetSigma(maskValue) << endl;
    cout << "Variance = " << LabelStatisticsfilter->GetVariance(maskValue) << endl;
    cout << "Sum     = " << LabelStatisticsfilter->GetSum(maskValue) << endl;
    // cout << "Bounding box = " << LabelStatisticsfilter->GetBoundingBox( maskValue )  << endl;
    cout << "Region  = " << LabelStatisticsfilter->GetRegion(maskValue) << endl;
    // cout << "HistogramPointer  = " << LabelStatisticsfilter->GetHistogram( maskValue )   << endl;

    outfile << "Minimum = " << LabelStatisticsfilter->GetMinimum(maskValue) << endl;
    outfile << "Maximum = " << LabelStatisticsfilter->GetMaximum(maskValue) << endl;
    outfile << "Median  = " << LabelStatisticsfilter->GetMedian(maskValue) << endl;
    outfile << "Mean = " << LabelStatisticsfilter->GetMean(maskValue) << endl;
    outfile << "Sigma = " << LabelStatisticsfilter->GetSigma(maskValue) << endl;
    outfile << "Variance = " << LabelStatisticsfilter->GetVariance(maskValue) << endl;
    outfile << "Sum = " << LabelStatisticsfilter->GetSum(maskValue) << endl;
    // outfile << "Bounding box = " << LabelStatisticsfilter->GetBoundingBox( maskValue )  << endl;
    // outfile << "Region  = " << LabelStatisticsfilter->GetRegion( maskValue )   << endl;
    // cout << "HistogramPointer  = " << LabelStatisticsfilter->GetHistogram( maskValue )   << endl;


    inputMin = static_cast<InputPixelType>(LabelStatisticsfilter->GetMinimum(maskValue));
    inputMax = static_cast<InputPixelType>(LabelStatisticsfilter->GetMaximum(maskValue));
} // CalcIntensityFeatures

void
CalcGeometryFeatures(InputImageType::Pointer inputImage, MaskImageType::Pointer maskImage, ostream &outfile)
{
    SpacingType spacing = inputImage->GetSpacing();

    using LabelGeometryType = itk::LabelGeometryImageFilter<MaskImageType, InputImageType>;
    LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();

    // Set up a connected components filter to label the binary objects.
    using ConnectedComponentType = itk::ConnectedComponentImageFilter<MaskImageType, MaskImageType>;
    ConnectedComponentType::Pointer connectedComponentFilter = ConnectedComponentType::New();
    connectedComponentFilter->SetInput(maskImage);

    // Relabel the components in order of size.
    using RelabelType = itk::RelabelComponentImageFilter<MaskImageType, MaskImageType>;
    RelabelType::Pointer relabeler = RelabelType::New();
    relabeler->SetInput(connectedComponentFilter->GetOutput());
    labelGeometryFilter->SetInput(relabeler->GetOutput());
    // labelGeometryFilter->SetInput(maskImage);

    // These generate optional outputs.
    labelGeometryFilter->CalculatePixelIndicesOn();
    labelGeometryFilter->CalculateOrientedBoundingBoxOn();
    labelGeometryFilter->CalculateOrientedLabelRegionsOn();
    labelGeometryFilter->CalculateOrientedIntensityRegionsOn();

    labelGeometryFilter->ReleaseDataFlagOn();

    labelGeometryFilter->SetIntensityInput(inputImage); //


    labelGeometryFilter->Update();
    cout << "Image features from the LabelGeometryImageFilter:" << endl;
    // outfile << "inputImage features from the LabelGeometryImageFilter:" << endl;

    outfile << "Volume = " << labelGeometryFilter->GetVolume(maskValue) * spacing[0] * spacing[1] * spacing[2] / 1000
            << endl;
    outfile << "IntegratedIntensity = " << labelGeometryFilter->GetIntegratedIntensity(maskValue) << endl;
    outfile << "Centroid = " << labelGeometryFilter->GetCentroid(maskValue) << endl;
    outfile << "WeightedCentroid = " << labelGeometryFilter->GetWeightedCentroid(maskValue) << endl;
    outfile << "AxesLength = " << labelGeometryFilter->GetAxesLength(maskValue) << endl;
    outfile << "MajorAxisLength = " << labelGeometryFilter->GetMajorAxisLength(maskValue) << endl;
    outfile << "MinorAxisLength = " << labelGeometryFilter->GetMinorAxisLength(maskValue) << endl;
    outfile << "Eccentricity = " << labelGeometryFilter->GetEccentricity(maskValue) << endl;
    outfile << "Elongation = " << labelGeometryFilter->GetElongation(maskValue) << endl;
    outfile << "Orientation = " << labelGeometryFilter->GetOrientation(maskValue) << endl;
    // outfile << "Boundingbox = " << labelGeometryFilter->GetBoundingBox(maskValue) << endl;

    outfile << "BoundingBoxVolume = " << labelGeometryFilter->GetBoundingBoxVolume(maskValue) * spacing[0]
        * spacing[1] * spacing[2] / 1000 << endl;
    outfile << "BoundingBoxSize = " << labelGeometryFilter->GetBoundingBoxSize(maskValue) << endl;
    outfile << "OrientedBoundingBoxVolume = " << labelGeometryFilter->GetOrientedBoundingBoxVolume(maskValue)
        * spacing[0] * spacing[1] * spacing[2] / 1000 << endl;
    outfile << "OrientedBoundingBoxSize = " << labelGeometryFilter->GetOrientedBoundingBoxSize(maskValue) << endl;

    // outfile << "Eigenvalues = " << labelGeometryFilter->GetEigenvalues(maskValue) << endl;
    outfile << "Eigenvalues= [" << labelGeometryFilter->GetEigenvalues(maskValue)[0] << ", "
            << labelGeometryFilter->GetEigenvalues(maskValue)[1] << ", "
            << labelGeometryFilter->GetEigenvalues(maskValue)[2] << "]" << endl;
    outfile << "Eigenvectors = " << labelGeometryFilter->GetEigenvectors(maskValue) << endl;
    outfile << "RotationMatrix = " << labelGeometryFilter->GetRotationMatrix(maskValue) << endl;
    // outfile << "\n\n";
} // CalcGeometryFeatures

void
CalcShapeAndIntensityFeatures(InputImageType::Pointer inputImage, MaskImageType::Pointer maskImage, ostream &outfile)
{
    using LabelObjectTypeS = itk::StatisticsLabelObject<LabelPixelType, Dimension>;
    using LabelMapTypeS = itk::LabelMap<LabelObjectTypeS>;

    // converting binary image to Statistics label map
    using I2LSType = itk::BinaryImageToStatisticsLabelMapFilter<MaskImageType, InputImageType, LabelMapTypeS>;
    I2LSType::Pointer i2ls = I2LSType::New();
    // itk::SimpleFilterWatcher watcher1( i2ls );

    i2ls->FullyConnectedOn();
    int inputForegroundValue  = 1;
    bool computeFeretDiameter = 1;
    unsigned int outputBackgroundValue = 0;
    bool computePerimeter     = 1;
    bool computeHistogram     = 1;
    unsigned int numberOfBins = NumberOfBins;
    i2ls->SetInputForegroundValue(inputForegroundValue);
    i2ls->SetOutputBackgroundValue(outputBackgroundValue);
    i2ls->SetComputeFeretDiameter(computeFeretDiameter);
    i2ls->ComputeFeretDiameterOff();
    i2ls->ComputeFeretDiameterOn();
    i2ls->SetComputePerimeter(computePerimeter);
    i2ls->ComputePerimeterOn();
    i2ls->SetComputeHistogram(computeHistogram);
    i2ls->ComputeHistogramOn();
    i2ls->SetNumberOfBins(numberOfBins);


    LabelMapTypeS::Pointer labelMapS;
    LabelObjectTypeS * labelObjectS;
    i2ls->SetFeatureImage(inputImage);
    i2ls->SetInput(maskImage);

    i2ls->Update();
    cout << "inputImage features from the BinaryImageToStatisticsLabelMapFilter:" << endl;
    // outfile << "inputImage features from the BinaryImageToStatisticsLabelMapFilter:" << endl;


    labelMapS    = i2ls->GetOutput();
    labelObjectS = labelMapS->GetLabelObject(maskValue);


    //    cout << "BoundingBox = " << labelObjectS->GetBoundingBox() << endl;
    cout << "CenterOfGravity = " << labelObjectS->GetCenterOfGravity() << endl;
    cout << "Centroid = " << labelObjectS->GetCentroid() << endl;
    cout << "Elongation = " << labelObjectS->GetElongation() << endl;
    //  cout << "EquivalentEllipsoidDiameter = " << labelObjectS->GetEquivalentEllipsoidDiameter() << endl;
    //  cout << "EquivalentSphericalPerimeter = " << labelObjectS->GetEquivalentSphericalPerimeter() << endl;
    //       cout << "EquivalentSphericalRadius = " << labelObjectS->GetEquivalentSphericalRadius() << endl;
    cout << "FeretDiameter = " << labelObjectS->GetFeretDiameter() << endl;
    // cout << "Histogram = " << labelObjectS->GetHistogram() << endl;
    //  cout << "Index = " << labelObjectS->GetIndex() << endl;
    cout << "Kurtosis = " << labelObjectS->GetKurtosis() << endl;
    cout << "Flatness = " << labelObjectS->GetFlatness() << endl;
    cout << "Label = " << labelObjectS->GetLabel() << endl;
    //  cout << "Line = " << labelObjectS->GetLine() << endl;
    //        cout << "LineContainer = " << labelObjectS->GetLineContainer() << endl;
    cout << "Maximum = " << labelObjectS->GetMaximum() << endl;
    // cout << "MaximumIndex = " << labelObjectS->GetMaximumIndex() << endl;
    cout << "Mean = " << labelObjectS->GetMean() << endl;
    cout << "Median = " << labelObjectS->GetMedian() << endl;
    cout << "Minimum = " << labelObjectS->GetMinimum() << endl;
    // cout << "MinimumIndex = " << labelObjectS->GetMinimumIndex() << endl;
    cout << "NameOfClass = " << labelObjectS->GetNameOfClass() << endl;
    cout << "NumberOfLines = " << labelObjectS->GetNumberOfLines() << endl;
    //        cout << "NumberOfPixels = " << labelObjectS->GetNumberOfPixels() << endl;
    //  cout << "NumberOfPixelsOnBorder = " << labelObjectS->GetNumberOfPixelsOnBorder() << endl;
    cout << "Perimeter = " << labelObjectS->GetPerimeter() << endl;
    //  cout << "PerimeterOnBorder = " << labelObjectS->GetPerimeterOnBorder() << endl;
    //        cout << "PerimeterOnBorderRatio = " << labelObjectS->GetPerimeterOnBorderRatio() << endl;
    cout << "PhysicalSize = " << labelObjectS->GetPhysicalSize() << endl;
    cout << "PrincipalAxes = " << labelObjectS->GetPrincipalAxes() << endl;
    cout << "PrincipalMoments = " << labelObjectS->GetPrincipalMoments() << endl;
    cout << "ReferenceCount = " << labelObjectS->GetReferenceCount() << endl;
    cout << "RegionElongation = " << labelObjectS->GetElongation() << endl; // RegionElongation is disappeared
    cout << "Roundness = " << labelObjectS->GetRoundness() << endl;
    cout << "SizeRegionRatio = " << 0 /*labelObjectS->GetSizeRegionRatio()*/ << endl; // This property is no longer provide
    cout << "Skewness = " << labelObjectS->GetSkewness() << endl;
    cout << "StandardDeviation = " << labelObjectS->GetStandardDeviation() << endl;
    cout << "Sum = " << labelObjectS->GetSum() << endl;
    cout << "Variance = " << labelObjectS->GetVariance() << endl;
    //  cout << "WeightedElongation = " << labelObjectS->GetWeightedElongation() << endl;
    //        cout << "WeightedFlatness = " << labelObjectS->GetWeightedFlatness() << endl;
    //  cout << "WeightedPrincipalAxes = " << labelObjectS->GetWeightedPrincipalAxes() << endl;
    //        cout << "WeightedPrincipalMoments = " << labelObjectS->GetWeightedPrincipalMoments() << endl;

    // outfile << "BoundingBox = " << labelObjectS->GetBoundingBox() << endl;
    outfile << "CenterOfGravity = " << labelObjectS->GetCenterOfGravity() << endl;
    outfile << "Centroid = " << labelObjectS->GetCentroid() << endl;
    outfile << "Elongation = " << labelObjectS->GetElongation() << endl;
    outfile << "EquivalentEllipsoidDiameter = " << labelObjectS->GetEquivalentEllipsoidDiameter() << endl;
    outfile << "EquivalentSphericalPerimeter = " << labelObjectS->GetEquivalentSphericalPerimeter() << endl;
    outfile << "EquivalentSphericalRadius = " << labelObjectS->GetEquivalentSphericalRadius() << endl;
    outfile << "FeretDiameter = " << labelObjectS->GetFeretDiameter() << endl;
    // outfile << "Histogram = " << labelObjectS->GetHistogram() << endl;
    // outfile << "Index = " << labelObjectS->GetIndex() << endl;
    outfile << "Kurtosis = " << labelObjectS->GetKurtosis() << endl;
    outfile << "Flatness = " << labelObjectS->GetFlatness() << endl;
    // outfile << "Label = " << labelObjectS->GetLabel() << endl;
    // outfile << "Line = " << labelObjectS->GetLine() << endl;
    // outfile << "LineContainer = " << labelObjectS->GetLineContainer() << endl;
    outfile << "Maximum = " << labelObjectS->GetMaximum() << endl;
    // outfile << "Maximum Index = " << labelObjectS->GetMaximumIndex() << endl;
    outfile << "Mean = " << labelObjectS->GetMean() << endl;
    outfile << "Median = " << labelObjectS->GetMedian() << endl;
    outfile << "Minimum = " << labelObjectS->GetMinimum() << endl;
    // outfile << "MinimumIndex = " << labelObjectS->GetMinimumIndex() << endl;
    // outfile << "NameOfClass = " << labelObjectS->GetNameOfClass() << endl;
    outfile << "NumberOfLines = " << labelObjectS->GetNumberOfLines() << endl;
    outfile << "NumberOfPixels = " << labelObjectS->GetNumberOfPixels() << endl;
    outfile << "NumberOfPixelsOnBorder = " << labelObjectS->GetNumberOfPixelsOnBorder() << endl;
    outfile << "Perimeter = " << labelObjectS->GetPerimeter() << endl;
    outfile << "PerimeterOnBorder = " << labelObjectS->GetPerimeterOnBorder() << endl;
    outfile << "PerimeterOnBorderRatio = " << labelObjectS->GetPerimeterOnBorderRatio() << endl;
    outfile << "PhysicalSize = " << labelObjectS->GetPhysicalSize() << endl;
    outfile << "PrincipalAxes = " << labelObjectS->GetPrincipalAxes() << endl;
    outfile << "PrincipalMoments = " << labelObjectS->GetPrincipalMoments() << endl;
    outfile << "ReferenceCount = " << labelObjectS->GetReferenceCount() << endl;
    outfile << "RegionElongation = " << labelObjectS->GetElongation() << endl; // deleted
    outfile << "Roundness = " << labelObjectS->GetRoundness() << endl;
    // outfile << "SizeRegionRatio = " << labelObjectS->GetSizeRegionRatio() << endl; // deleted
    outfile << "Skewness = " << labelObjectS->GetSkewness() << endl;
    outfile << "StandardDeviation = " << labelObjectS->GetStandardDeviation() << endl;
    outfile << "Sum = " << labelObjectS->GetSum() << endl;
    outfile << "Variance = " << labelObjectS->GetVariance() << endl;
    outfile << "WeightedElongation = " << labelObjectS->GetWeightedElongation() << endl;
    outfile << "WeightedFlatness = " << labelObjectS->GetWeightedFlatness() << endl;
    outfile << "WeightedPrincipalAxes = " << labelObjectS->GetWeightedPrincipalAxes() << endl;
    outfile << "WeightedPrincipalMoments = " << labelObjectS->GetWeightedPrincipalMoments() << endl;
} // CalcShapeAndIntensityFeatures

void
CalcShapeAndIntensityFeatures(InputImage2DType::Pointer inputImage, MaskImage2DType::Pointer maskImage,
  ostream &outfile)
{
    using LabelObject2DTypeS = itk::StatisticsLabelObject<LabelPixelType, Dimension - 1>;
    using LabelMap2DTypeS = itk::LabelMap<LabelObject2DTypeS>;

    // converting binary image to Statistics label map
    using I2LS2DType = itk::BinaryImageToStatisticsLabelMapFilter<MaskImage2DType, InputImage2DType, LabelMap2DTypeS>;
    I2LS2DType::Pointer i2ls = I2LS2DType::New();
    // itk::SimpleFilterWatcher watcher1( i2ls );

    i2ls->FullyConnectedOn();
    int inputForegroundValue  = 1;
    bool computeFeretDiameter = 1;
    unsigned int outputBackgroundValue = 0;
    bool computePerimeter     = 1;
    bool computeHistogram     = 1;
    unsigned int numberOfBins = NumberOfBins;
    i2ls->SetInputForegroundValue(inputForegroundValue);
    i2ls->SetOutputBackgroundValue(outputBackgroundValue);
    i2ls->SetComputeFeretDiameter(computeFeretDiameter);
    i2ls->ComputeFeretDiameterOff();
    i2ls->ComputeFeretDiameterOn();
    i2ls->SetComputePerimeter(computePerimeter);
    i2ls->ComputePerimeterOn();
    i2ls->SetComputeHistogram(computeHistogram);
    i2ls->ComputeHistogramOn();
    i2ls->SetNumberOfBins(numberOfBins);


    LabelMap2DTypeS::Pointer labelMapS;
    LabelObject2DTypeS * labelObjectS;
    i2ls->SetFeatureImage(inputImage);
    i2ls->SetInput(maskImage);

    i2ls->Update();
    cout << "inputImage features from the BinaryImageToStatisticsLabelMapFilter:" << endl;
    // outfile << "inputImage features from the BinaryImageToStatisticsLabelMapFilter:" << endl;


    labelMapS    = i2ls->GetOutput();
    labelObjectS = labelMapS->GetLabelObject(maskValue);


    //    cout << "2D BoundingBox = " << labelObjectS->GetBoundingBox() << endl;
    cout << "2D CenterOfGravity = " << labelObjectS->GetCenterOfGravity() << endl;
    cout << "2D Centroid = " << labelObjectS->GetCentroid() << endl;
    cout << "2D Elongation = " << labelObjectS->GetElongation() << endl;
    cout << "2D EquivalentEllipsoidDiameter = " << labelObjectS->GetEquivalentEllipsoidDiameter() << endl;
    cout << "2D EquivalentSphericalPerimeter = " << labelObjectS->GetEquivalentSphericalPerimeter() << endl;
    cout << "2D EquivalentSphericalRadius = " << labelObjectS->GetEquivalentSphericalRadius() << endl;
    cout << "2D FeretDiameter = " << labelObjectS->GetFeretDiameter() << endl;
    // cout << "2D Histogram = " << labelObjectS->GetHistogram() << endl;
    //  cout << "2D Index = " << labelObjectS->GetIndex() << endl;
    cout << "2D Kurtosis = " << labelObjectS->GetKurtosis() << endl;
    cout << "2D Flatness = " << labelObjectS->GetFlatness() << endl;
    // cout << "2D Label = " << labelObjectS->GetLabel() << endl;
    //  cout << "2D Line = " << labelObjectS->GetLine() << endl;
    //        cout << "2D LineContainer = " << labelObjectS->GetLineContainer() << endl;
    cout << "2D Maximum = " << labelObjectS->GetMaximum() << endl;
    // cout << "2D MaximumIndex = " << labelObjectS->GetMaximumIndex() << endl;
    cout << "2D Mean = " << labelObjectS->GetMean() << endl;
    cout << "2D Median = " << labelObjectS->GetMedian() << endl;
    cout << "2D Minimum = " << labelObjectS->GetMinimum() << endl;
    // cout << "2D MinimumIndex = " << labelObjectS->GetMinimumIndex() << endl;
    cout << "2D NameOfClass = " << labelObjectS->GetNameOfClass() << endl;
    cout << "2D NumberOfLines = " << labelObjectS->GetNumberOfLines() << endl;
    //        cout << "2D NumberOfPixels = " << labelObjectS->GetNumberOfPixels() << endl;
    //  cout << "2D NumberOfPixelsOnBorder = " << labelObjectS->GetNumberOfPixelsOnBorder() << endl;
    cout << "2D Perimeter = " << labelObjectS->GetPerimeter() << endl;
    cout << "2D PerimeterOnBorder = " << labelObjectS->GetPerimeterOnBorder() << endl;
    cout << "2D PerimeterOnBorderRatio = " << labelObjectS->GetPerimeterOnBorderRatio() << endl;
    cout << "2D PhysicalSize = " << labelObjectS->GetPhysicalSize() << endl;
    cout << "2D PrincipalAxes = " << labelObjectS->GetPrincipalAxes() << endl;
    cout << "2D PrincipalMoments = " << labelObjectS->GetPrincipalMoments() << endl;
    cout << "2D ReferenceCount = " << labelObjectS->GetReferenceCount() << endl;
    cout << "2D RegionElongation = " << labelObjectS->GetElongation() << endl; // RegionElongation is disappeared
    cout << "2D Roundness = " << labelObjectS->GetRoundness() << endl;
    cout << "2D SizeRegionRatio = " << 0 /*labelObjectS->GetSizeRegionRatio()*/ << endl; // This property is no longer provide
    cout << "2D Skewness = " << labelObjectS->GetSkewness() << endl;
    cout << "2D StandardDeviation = " << labelObjectS->GetStandardDeviation() << endl;
    cout << "2D Sum = " << labelObjectS->GetSum() << endl;
    cout << "2D Variance = " << labelObjectS->GetVariance() << endl;
    //  cout << "2D WeightedElongation = " << labelObjectS->GetWeightedElongation() << endl;
    //        cout << "2D WeightedFlatness = " << labelObjectS->GetWeightedFlatness() << endl;
    //  cout << "2D WeightedPrincipalAxes = " << labelObjectS->GetWeightedPrincipalAxes() << endl;
    //        cout << "2D WeightedPrincipalMoments = " << labelObjectS->GetWeightedPrincipalMoments() << endl;

    // outfile << "2DBoundingBox = " << labelObjectS->GetBoundingBox() << endl;
    outfile << "2DCenterOfGravity = " << labelObjectS->GetCenterOfGravity() << endl;
    outfile << "2DCentroid = " << labelObjectS->GetCentroid() << endl;
    outfile << "2DElongation = " << labelObjectS->GetElongation() << endl;
    outfile << "2DEquivalentEllipsoidDiameter = " << labelObjectS->GetEquivalentEllipsoidDiameter() << endl;
    outfile << "2DEquivalentSphericalPerimeter = " << labelObjectS->GetEquivalentSphericalPerimeter() << endl;
    outfile << "2DEquivalentSphericalRadius = " << labelObjectS->GetEquivalentSphericalRadius() << endl;
    outfile << "2DFeretDiameter = " << labelObjectS->GetFeretDiameter() << endl;
    // outfile << "2DHistogram = " << labelObjectS->GetHistogram() << endl;
    // outfile << "2DIndex = " << labelObjectS->GetIndex() << endl;
    outfile << "2DKurtosis = " << labelObjectS->GetKurtosis() << endl;
    outfile << "2DFlatness = " << labelObjectS->GetFlatness() << endl;
    // outfile << "2DLabel = " << labelObjectS->GetLabel() << endl;
    // outfile << "2DLine = " << labelObjectS->GetLine() << endl;
    // outfile << "2DLineContainer = " << labelObjectS->GetLineContainer() << endl;
    outfile << "2DMaximum = " << labelObjectS->GetMaximum() << endl;
    // outfile << "2DMaximum Index = " << labelObjectS->GetMaximumIndex() << endl;
    outfile << "2DMean = " << labelObjectS->GetMean() << endl;
    outfile << "2DMedian = " << labelObjectS->GetMedian() << endl;
    outfile << "2DMinimum = " << labelObjectS->GetMinimum() << endl;
    // outfile << "2DMinimumIndex = " << labelObjectS->GetMinimumIndex() << endl;
    // outfile << "2DNameOfClass = " << labelObjectS->GetNameOfClass() << endl;
    outfile << "2DNumberOfLines = " << labelObjectS->GetNumberOfLines() << endl;
    outfile << "2DNumberOfPixels = " << labelObjectS->GetNumberOfPixels() << endl;
    outfile << "2DNumberOfPixelsOnBorder = " << labelObjectS->GetNumberOfPixelsOnBorder() << endl;
    outfile << "2DPerimeter = " << labelObjectS->GetPerimeter() << endl;
    outfile << "2DPerimeterOnBorder = " << labelObjectS->GetPerimeterOnBorder() << endl;
    outfile << "2DPerimeterOnBorderRatio = " << labelObjectS->GetPerimeterOnBorderRatio() << endl;
    outfile << "2DPhysicalSize = " << labelObjectS->GetPhysicalSize() << endl;
    outfile << "2DPrincipalAxes = " << labelObjectS->GetPrincipalAxes() << endl;
    outfile << "2DPrincipalMoments = " << labelObjectS->GetPrincipalMoments() << endl;
    outfile << "2DReferenceCount = " << labelObjectS->GetReferenceCount() << endl;
    outfile << "2DRegionElongation = " << labelObjectS->GetElongation() << endl; // deleted
    outfile << "2DRoundness = " << labelObjectS->GetRoundness() << endl;
    // outfile << "2DSizeRegionRatio = " << labelObjectS->GetSizeRegionRatio() << endl; // deleted
    outfile << "2DSkewness = " << labelObjectS->GetSkewness() << endl;
    outfile << "2DStandardDeviation = " << labelObjectS->GetStandardDeviation() << endl;
    outfile << "2DSum = " << labelObjectS->GetSum() << endl;
    outfile << "2DVariance = " << labelObjectS->GetVariance() << endl;
    outfile << "2DWeightedElongation = " << labelObjectS->GetWeightedElongation() << endl;
    outfile << "2DWeightedFlatness = " << labelObjectS->GetWeightedFlatness() << endl;
    outfile << "2DWeightedPrincipalAxes = " << labelObjectS->GetWeightedPrincipalAxes() << endl;
    outfile << "2DWeightedPrincipalMoments = " << labelObjectS->GetWeightedPrincipalMoments() << endl;
} // CalcShapeAndIntensityFeatures

void
CalcGlcmFeatures(InputImageType::Pointer inputImage, InputImageType::Pointer maskImage, InputPixelType inputMin,
  InputPixelType inputMax, ostream &outfile)
{
    using TextureFilterType = itk::Statistics::ScalarImageToTextureFeaturesFilter<InputImageType>;
    TextureFilterType::Pointer glcm = TextureFilterType::New();

    using TextureFeaturesFilterType = TextureFilterType::TextureFeaturesFilterType;

    TextureFilterType::FeatureNameVectorPointer requestedFeatures =
      TextureFilterType::FeatureNameVector::New();

#if ((ITK_VERSION_MAJOR == 5) && (ITK_VERSION_MINOR >= 1))
    using TextureFeatureEnums = itk::Statistics::HistogramToTextureFeaturesFilterEnums::TextureFeature;
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::Energy));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::Entropy));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::Correlation));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::InverseDifferenceMoment));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::Inertia));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::ClusterShade));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::ClusterProminence));
    requestedFeatures->push_back(static_cast<uint8_t>(TextureFeatureEnums::HaralickCorrelation));
#else
    requestedFeatures->push_back(TextureFeaturesFilterType::Energy);
    requestedFeatures->push_back(TextureFeaturesFilterType::Entropy);
    requestedFeatures->push_back(TextureFeaturesFilterType::Correlation);
    requestedFeatures->push_back(TextureFeaturesFilterType::InverseDifferenceMoment);
    requestedFeatures->push_back(TextureFeaturesFilterType::Inertia);
    requestedFeatures->push_back(TextureFeaturesFilterType::ClusterShade);
    requestedFeatures->push_back(TextureFeaturesFilterType::ClusterProminence);
    requestedFeatures->push_back(TextureFeaturesFilterType::HaralickCorrelation);
#endif

    glcm->SetRequestedFeatures(requestedFeatures);
    // glcm->SetNormalizeOn();
    // glcm->Print( cout );


    glcm->SetInput(inputImage);

    glcm->SetMaskImage(maskImage);
    // glcm->FastCalculationsOn();

    glcm->SetNumberOfBinsPerAxis(NumberOfBins);
    glcm->SetPixelValueMinMax(inputMin, inputMax);

    glcm->FastCalculationsOff();
    glcm->Update();
    cout << "The inputImage features from the ScalarImageToTextureFeaturesFilter:" << endl;
    // outfile<< endl<< endl << "The inputImage features from the ScalarImageToTextureFeaturesFilter:" << endl;


    using GLCMFilterType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InputImageType>;
    using HistogramType = GLCMFilterType::HistogramType;
    using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType>;

    TextureFilterType::OffsetVector::ConstPointer offsets = glcm->GetOffsets();
    for (TextureFilterType::OffsetVector::ConstIterator offset = offsets->Begin(); offset != offsets->End(); offset++) {
        GLCMFilterType::Pointer glcmGenerator = GLCMFilterType::New();
        glcmGenerator->SetOffset(offset->Value());
        glcmGenerator->SetNumberOfBinsPerAxis(NumberOfBins);    // reasonable number of bins
        glcmGenerator->SetPixelValueMinMax(inputMin, inputMax); // for input UCHAR pixel type
        Hist2FeaturesType::Pointer featureCalc = Hist2FeaturesType::New();

        glcmGenerator->SetInput(inputImage);
        glcmGenerator->SetMaskImage(maskImage);
        glcmGenerator->Update();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        cout << offset->Value() << "\t" << featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment) << endl;
    }


    int fsize = glcm->GetFeatureMeans()->Size();
    for (int x = 0; x < fsize; x++) {
        char charactheristic[100];
        decodeCharactheristicGLCM(glcm->GetRequestedFeatures()->GetElement(x), charactheristic);

        cout << "MeanOf" << charactheristic << "= " << glcm->GetFeatureMeans()->GetElement(x) << endl;

        outfile << "MeanOf" << charactheristic << "= " << glcm->GetFeatureMeans()->GetElement(x) << endl;
    }
    for (int x = 0; x < fsize; x++) {
        char charactheristic[100];
        decodeCharactheristicGLCM(glcm->GetRequestedFeatures()->GetElement(x), charactheristic);

        cout << "StandardDeviationOf" << charactheristic << "= "
             << glcm->GetFeatureStandardDeviations()->GetElement(x) << endl;

        outfile << "StandardDeviationOf" << charactheristic << "= " << glcm->GetFeatureStandardDeviations()->GetElement(
            x) << endl;
    }
} // CalcGlcmFeatures

void
CalcGlrmFeatures(InputImageType::Pointer inputImage, InputImageType::Pointer maskImage, InputPixelType inputMin,
  InputPixelType inputMax, double distMax, ostream &outfile)
{
    using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter<InputImageType>;
    RunLengthFilterType::Pointer glrm = RunLengthFilterType::New();

    using RunLengthFeaturesFilterType = RunLengthFilterType::RunLengthFeaturesFilterType;

    RunLengthFilterType::FeatureNameVectorPointer requestedFeatures =
      RunLengthFilterType::FeatureNameVector::New();
#if ((ITK_VERSION_MAJOR == 5) && (ITK_VERSION_MINOR >= 1))
    using RunLengthFeatureEnums = itk::Statistics::HistogramToRunLengthFeaturesFilterEnums::RunLengthFeature;
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::ShortRunEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::LongRunEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::GreyLevelNonuniformity));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::RunLengthNonuniformity));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::LowGreyLevelRunEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::HighGreyLevelRunEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::ShortRunLowGreyLevelEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::ShortRunHighGreyLevelEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::LongRunLowGreyLevelEmphasis));
    requestedFeatures->push_back(static_cast<uint8_t>(RunLengthFeatureEnums::LongRunHighGreyLevelEmphasis));
#else
    requestedFeatures->push_back(RunLengthFeaturesFilterType::ShortRunEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::LongRunEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::GreyLevelNonuniformity);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::RunLengthNonuniformity);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::LowGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::HighGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::ShortRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::ShortRunHighGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::LongRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeaturesFilterType::LongRunHighGreyLevelEmphasis);
#endif
    glrm->SetRequestedFeatures(requestedFeatures);


    glrm->SetInput(inputImage);
    // glrm->SetInput(maskImage);
    glrm->SetMaskImage(maskImage);
    // glrm->FastCalculationsOn();

    glrm->SetNumberOfBinsPerAxis(NumberOfBins);
    glrm->SetPixelValueMinMax(inputMin, (inputMax - inputMin) * 2 + inputMin);
    glrm->SetDistanceValueMinMax(0, distMax * 2 / 3);

    glrm->FastCalculationsOff();
    glrm->Update();
    cout << "The inputImage features from the ScalarImageToRunLengthFeaturesFilter:" << endl;
    // outfile<< endl<< endl << "The inputImage features from the ScalarImageToTextureFeaturesFilter:" << endl;


    int size = glrm->GetFeatureMeans()->Size();
    for (int x = 0; x < size; x++) {
        char charactheristic[100];
        decodeCharactheristicGLRM(glrm->GetRequestedFeatures()->GetElement(x), charactheristic);

        cout << "MeanOf" << charactheristic << "= " << glrm->GetFeatureMeans()->GetElement(x) << endl;

        outfile << "MeanOf" << charactheristic << "= " << glrm->GetFeatureMeans()->GetElement(x) << endl;
    }
    for (int x = 0; x < size; x++) {
        char charactheristic[100];
        decodeCharactheristicGLRM(glrm->GetRequestedFeatures()->GetElement(x), charactheristic);

        cout << "StandardDeviationOf" << charactheristic << "= "
             << glrm->GetFeatureStandardDeviations()->GetElement(x) << endl;

        outfile << "StandardDeviationOf" << charactheristic << "= " << glrm->GetFeatureStandardDeviations()->GetElement(
            x) << endl;
    }
} // CalcGlrmFeatures

MaskImageType::Pointer
StentRemoval(InputImageType::Pointer inputImage, MaskImageType::Pointer maskImage, InputPixelType stentThreshold = 500)
{
    StructuringElementType structuringElement1;

    structuringElement1.SetRadius(1);
    structuringElement1.CreateStructuringElement();

    StructuringElementType structuringElement3;
    structuringElement3.SetRadius(3);
    structuringElement3.CreateStructuringElement();

    StructuringElementType structuringElement5;
    structuringElement5.SetRadius(5);
    structuringElement5.CreateStructuringElement();

    StructuringElementType structuringElement10;
    structuringElement10.SetRadius(10);
    structuringElement10.CreateStructuringElement();

    #if 0
    ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
    erodeFilter->SetInput(maskImage);
    erodeFilter->SetKernel(structuringElement3);
    erodeFilter->SetErodeValue(maskValue);
    erodeFilter->Update();

    maskImage = erodeFilter->GetOutput();
    #endif


    using BinaryThresholdfilterType = itk::BinaryThresholdImageFilter<InputImageType, MaskImageType>;

    BinaryThresholdfilterType::Pointer binaryThresholdfilter = BinaryThresholdfilterType::New();

    binaryThresholdfilter->SetInput(inputImage);
    binaryThresholdfilter->SetInsideValue(1);
    binaryThresholdfilter->SetOutsideValue(0);
    binaryThresholdfilter->SetLowerThreshold(500); // <- Main paraneter

    #if 0
    using LabelOpeningType = itk::BinaryShapeKeepNObjectsImageFilter<MaskImageType>;

    LabelOpeningType::Pointer opening = LabelOpeningType::New();
    opening->SetInput(binaryThresholdfilter->GetOutput());
    opening->SetBackgroundValue(0);
    opening->SetForegroundValue(maskValue);
    opening->SetNumberOfObjects(1);
    // opening->SetReverseOrdering( false );
    // opening->SetAttribute(100); //  extract the largest object  NUMBER_OF_PIXELS
    opening->Update();
    #endif // if 0

    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(binaryThresholdfilter->GetOutput());
    dilateFilter->SetKernel(structuringElement3);
    dilateFilter->SetDilateValue(maskValue);
    dilateFilter->Update();

    using MaskImage2DType = itk::Image<MaskPixelType, Dimension - 1>;

    using SliceBySliceFilterType = itk::SliceBySliceImageFilter<MaskImageType, MaskImageType>;
    SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

    using I2LType = itk::BinaryFillholeImageFilter<MaskImage2DType>;
    I2LType::Pointer reconstruction = I2LType::New();
    reconstruction->SetFullyConnected(true);
    reconstruction->SetForegroundValue(maskValue);

    sliceBySliceFilter->SetInput(dilateFilter->GetOutput());
    sliceBySliceFilter->SetFilter(reconstruction);

    sliceBySliceFilter->Update();


    // remove stent area
    using SubtractImageFilterType = itk::SubtractImageFilter<MaskImageType, MaskImageType>;
    SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
    subtractFilter->SetInput1(maskImage);
    subtractFilter->SetInput2(sliceBySliceFilter->GetOutput());

    using AndImageFilterType = itk::AndImageFilter<MaskImageType>;
    AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
    andFilter->SetInput1(maskImage);
    andFilter->SetInput2(subtractFilter->GetOutput());
    andFilter->Update();
    maskImage = andFilter->GetOutput();

    return maskImage;
} // StentRemoval

int
ExtractLargestAreaSlice(InputImageType::Pointer inputImage, MaskImageType::Pointer maskImage,
  InputImage2DType::Pointer& input2DImage, MaskImage2DType::Pointer& mask2DImage)
{
    MaskImageType::RegionType maskRegion = maskImage->GetLargestPossibleRegion();
    MaskImageType::IndexType start       = maskRegion.GetIndex();
    MaskImageType::SizeType size         = maskRegion.GetSize();

    const unsigned int numberofSlices = size[2];
    size[2] = 0;


    using FilterType = itk::ExtractImageFilter<MaskImageType, MaskImage2DType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetDirectionCollapseToSubmatrix();

    float maxArea         = 0;
    unsigned int maxSlice = 0;
    for (unsigned int i = 0; i < numberofSlices; i++) {
        start[2] = i;

        MaskImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);

        filter->SetExtractionRegion(desiredRegion);
        filter->SetInput(maskImage);
        filter->Update();

        using maskImage2DIteratorType = itk::ImageRegionConstIterator<MaskImage2DType>;
        maskImage2DIteratorType maskImage2Dit(filter->GetOutput(), filter->GetOutput()->GetRequestedRegion());

        float curArea = 0;
        maskImage2Dit.GoToBegin();
        while (!maskImage2Dit.IsAtEnd()) {
            if (maskImage2Dit.Get() == maskValue) ++curArea;
            ++maskImage2Dit;
        }
        if (curArea > maxArea) {
            maxArea  = curArea;
            maxSlice = i;
        }
        cout << i << ": max area " << maxArea << ", current area " << curArea << endl;
    }


    // applying the slice extraction for inputImage
    start[2] = maxSlice;

    MaskImageType::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);

    filter->SetExtractionRegion(desiredRegion);
    filter->SetInput(maskImage);
    filter->Update();
    mask2DImage = filter->GetOutput();


    // applying the slice extraction for inputImage
    {
        // InputImageType::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
        // InputImageType::IndexType start = inputRegion.GetIndex();
        // InputImageType::SizeType size = inputRegion.GetSize();
        //
        // size[2] = 0;
        // start[2] = maxSlice;

        using FilterType = itk::ExtractImageFilter<InputImageType, InputImage2DType>;
        FilterType::Pointer filter = FilterType::New();
        // filter->InPlaceOn();
        filter->SetDirectionCollapseToSubmatrix();

        cout << start << size << endl;

        InputImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);

        filter->SetExtractionRegion(desiredRegion);
        filter->SetInput(inputImage);
        filter->Update();
        input2DImage = filter->GetOutput();
    }

    return 0;
} // ExtractLargestAreaSlice

int
main(int argc, char * argv[])
{
    bool is2DShapeIntensity = false; // 2
    bool isIntensity        = false; // i
    bool isShapeIntensity   = false; // s
    bool isGeometry         = false; // g
    bool isGLCM = false;             // c
    bool isGLRM = false;             // r
    string featureSelect = "2isgcr";


    if (argc < 3) {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage = " << argv[0];
        cerr << " InputImage LabelImage FeatureFile Label={1} FeatureSelect={2isgcr} NumberOfBins={64}" << endl;
        cerr << " 2 - 2DShapeIntensity, i - Intensity, s - Shape Intensity, g - Geometry, c - GLCM, r - GLRM" << endl;
        return EXIT_FAILURE;
    }

    LabelPixelType selectedLabelValue = maskValue;
    if (argc > 4) {
        selectedLabelValue = atoi(argv[4]);
    }

    if (argc > 5) {
        featureSelect = argv[5];
    }

    for (int i = 0; i < featureSelect.length(); i++) {
        switch (featureSelect[i]) {
            case '2':
                is2DShapeIntensity = true; // 2
                break;
            case 'i':
                isIntensity = true; // i
                break;
            case 's':
                isShapeIntensity = true; // s
                break;
            case 'g':
                isGeometry = true; // g
                break;
            case 'c':
                isIntensity = true; // c needs i
                isGLCM      = true; // c
                break;
            case 'r':
                isIntensity = true; // r needs i
                isGLRM      = true; // r
                break;
        }
    }

    string inputImageName  = argv[1];
    string labelImageName  = argv[2];
    string featureFileName = argv[3];

    if (argc > 6) {
        NumberOfBins = atoi(argv[6]);
    }


    cout << "Input Image Name = " << inputImageName << endl;
    cout << "Label Image Name = " << labelImageName << endl;
    cout << "Feature File Name = " << featureFileName << endl;


    /////////To Read CT images ////////////////////
    InputImageType::Pointer inputImage = ReadImageFile<InputImageType>(argv[1]);

    ///////// To Read the Mask images ////////////////////
    cout << "Get the mask from Input" << endl;
    LabelImageType::Pointer labelImage = ReadImageFile<LabelImageType>(argv[2]);

    // To check coordinate of the inputImage Image
    SpacingType inputImageSpacing     = inputImage->GetSpacing();
    OriginType inputImageOrigin       = inputImage->GetOrigin();
    RegionType inputImageRegion       = inputImage->GetLargestPossibleRegion();
    SizeType inputImageSize           = inputImageRegion.GetSize();
    DirectionType inputImageDirection = inputImage->GetDirection();

    // To check coordinate For mask Image
    SpacingType labelImageSpacing     = labelImage->GetSpacing();
    OriginType labelImageOrigin       = labelImage->GetOrigin();
    RegionType labelImageRegion       = labelImage->GetLargestPossibleRegion();
    SizeType labelImageRegionSize     = labelImageRegion.GetSize();
    DirectionType labelImageDirection = labelImage->GetDirection();

    /////////////////////////////

    cout << "Input Image Spacing = " << inputImageSpacing << endl;
    cout << "Input Image Origin = " << inputImageOrigin << endl;
    cout << "Input Image Size = " << inputImageSize << endl;
    cout << "Input Image Direction = " << inputImageDirection << endl << endl << endl;


    cout << "Label Image Spacing = " << labelImageSpacing << endl;
    cout << "Label Image Origin = " << labelImageOrigin << endl;
    cout << "Label Image RegionSize = " << labelImageRegionSize << endl;
    cout << "Label Image Direction = " << labelImageDirection << endl << endl << endl;

    // using TransformType = itk::IdentityTransform<double, Dimension>;

    /*
     * using ResampleImageFilterType = itk::ResampleImageFilter<InputImageType, InputImageType>;
     * ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
     * resampler->SetInput(inputImage);
     * resampler->SetSize(labelImageRegionSize);
     * resampler->SetOutputOrigin(labelImageOrigin);
     * resampler->SetOutputSpacing(labelImageSpacing);
     * resampler->SetOutputDirection(labelImageDirection);
     * //resampler->SetTransform(TransformType::New());
     * resampler->UpdateLargestPossibleRegion();
     * inputImage = resampler->GetOutput();
     *
     * inputImageSpacing = inputImage->GetSpacing();
     * inputImageOrigin  = inputImage->GetOrigin();
     * inputImageRegion  = inputImage->GetLargestPossibleRegion();
     * inputImageSize    = inputImageRegion.GetSize();
     *
     * cout << "Resampled Input Image Spacing = " << inputImageSpacing << endl;
     * cout << "Resampled Input Image Origin = " << inputImageOrigin << endl;
     * cout << "Resampled Input Image Size = " << inputImageSize << endl;
     * cout << "Resampled Input Image Direction = " << inputImageDirection << endl<< endl << endl;
     *
     * WriteImageFile<InputImageType>(inputImage, "test.nrrd");*/

    // Image ROI
    InputImageType::RegionType inputRegion;
    MaskImageType::RegionType maskRegion;

    MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->SetRegions(labelImage->GetRequestedRegion());
    maskImage->CopyInformation(labelImage);
    maskImage->Allocate();

    maskImage = GetMaskImage(labelImage, selectedLabelValue);

    // For Island Removing
    MaskImageType::Pointer maskImageIslandRemoved;
    {
        // using LabelOpeningType = itk::LabelShapeKeepNObjectsImageFilter< MaskImageType >;
        using LabelOpeningType = itk::BinaryShapeKeepNObjectsImageFilter<MaskImageType>;

        LabelOpeningType::Pointer opening = LabelOpeningType::New();
        opening->SetInput(maskImage);
        opening->SetBackgroundValue(0);
        opening->SetFullyConnected(true);
        opening->SetForegroundValue(1);
        opening->SetNumberOfObjects(1);
        // opening->SetReverseOrdering( false );
        opening->SetAttribute(LabelOpeningType::LabelObjectType::NUMBER_OF_PIXELS); //  extract the largest object  NUMBER_OF_PIXELS
        opening->Update();

        maskImageIslandRemoved = opening->GetOutput();
    }
    maskImage = maskImageIslandRemoved;

    maskRegion = GetRoi(maskImage);
    // ExpandRoi(maskImage, maskRegion, 1);
    BoundingCheck(maskImage, inputImage, maskRegion, inputRegion);

    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;

    RegionToIndex(maskRegion, roiStart, roiEnd);

    inputImage = ApplyRoi<InputImageType>(inputImage, inputRegion);
    maskImage  = ApplyRoi<MaskImageType>(maskImage, maskRegion);

    // To check coordinate For mask Image
    SpacingType maskImageSpacing     = maskImage->GetSpacing();
    OriginType maskImageOrigin       = maskImage->GetOrigin();
    RegionType maskImageRegion       = maskImage->GetLargestPossibleRegion();
    SizeType maskImageRegionSize     = maskRegion.GetSize();
    DirectionType maskImageDirection = maskImage->GetDirection();

    cout << "Mask Image Spacing = " << maskImageSpacing << endl;
    cout << "Mask Image Origin = " << maskImageOrigin << endl;
    cout << "Mask Image RegionSize = " << maskImageRegionSize << endl;
    cout << "Mask Image Direction = " << maskImageDirection << endl << endl << endl;

    inputImage->SetSpacing(maskImage->GetSpacing());
    inputImage->SetOrigin(maskImage->GetOrigin());
    inputImage->SetDirection(maskImage->GetDirection());

    InternalImageType::Pointer maskImageFloat;
    using CastToInputFilterType = itk::CastImageFilter<MaskImageType, InputImageType>;
    CastToInputFilterType::Pointer toInputFilter = CastToInputFilterType::New();
    toInputFilter->SetInput(maskImage);
    toInputFilter->Update();
    maskImageFloat = toInputFilter->GetOutput();


    // InternalImageType::Pointer inputImageFloat = InternalImageType::New();
    //    inputImageFloat->SetRegions( inputImage->GetRequestedRegion() );
    //    inputImageFloat->CopyInformation( inputImage );
    //    inputImageFloat->Allocate();

    //    {
    //        // Changing the data type of Input Image (unsigned int)  to float

    //        using inputImageIteratorFloatType = itk::ImageRegionIterator< InternalImageType>;
    //        inputImageIteratorFloatType inputImageFloatIt( inputImageFloat, inputImageFloat->GetRequestedRegion() );

    //        using inputImageIteratorType = itk::ImageRegionConstIterator< InputImageType>;
    //        inputImageIteratorType inputImageit( inputImage, inputImage->GetRequestedRegion() );

    //        inputImageit.GoToBegin();
    //        inputImageFloatIt.GoToBegin();

    //        while ( !inputImageit.IsAtEnd() )
    //        {
    //            inputImageFloatIt.Set(  inputImageit.Get() );
    //            ++inputImageit;
    //            ++inputImageFloatIt;
    //        }
    //    }


    // Output file open
    ofstream outfile((featureFileName).c_str(), ios::out);
    if (!outfile) {
        cerr << "Can't open output file " << featureFileName << endl;
        exit(1);
    }

    #if 0
    using SUVCastFilterType = itk::CastImageFilter<InputImageType, LabelImageType>;
    SUVCastFilterType::Pointer SUVCaster = SUVCastFilterType::New();

    SUVCaster->SetInput(inputImage);
    SUVCaster->Update();
    #endif


    #if 0
    {
        MaskWriterType::Pointer binaryWriter = MaskWriterType::New();

        binaryWriter->SetFileName("test0.dcm");
        binaryWriter->SetInput(maskImage);
        binaryWriter->Update();
    }
    #endif


    // maskImage = StentRemoval(inputImage, maskImage, 500);

    ///////////////////////////////////////////////////////////////////////////
    // Compute 2D shape and intensity feature uisng LabelMap   ////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if (is2DShapeIntensity) {
        // extract the slice having the largest area
        MaskImage2DType::Pointer mask2DImage;
        InputImage2DType::Pointer input2DImage;
        ExtractLargestAreaSlice(inputImage, maskImage, input2DImage, mask2DImage);
        CalcShapeAndIntensityFeatures(input2DImage, mask2DImage, outfile);
    }


    /////////////////////////////////////////////////////////////////
    // Compute Intensity Features for three SUV images using itkLabelStatisticsImageFilter /////////
    ///////////////////////////////////////////////////////////////
    InputPixelType inputMin;
    InputPixelType inputMax;
    if (isIntensity)
        CalcIntensityFeatures(inputImage, maskImage, inputMin, inputMax, outfile);

    /////////////////////////////////////////////////////////////////////
    // Compute Geometry feature using LabelGeometryImageFilter  ////
    /////////////////////////////////////////////////////////////////////
    if (isGeometry)
        CalcGeometryFeatures(inputImage, maskImage, outfile);


    ///////////////////////////////////////////////////////////////////////////
    // Compute shape and intensity feature uisng LabelMap   ////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if (isShapeIntensity)
        CalcShapeAndIntensityFeatures(inputImage, maskImage, outfile);


    ////////////////////////////////////////////////////////////////////
    ///// Compute Feature for the Input images using ScalarImageToTextureFeaturesFilter /////
    ////////////////////////////////////////////////////////////////////
    if (isIntensity && isGLCM)
        CalcGlcmFeatures(inputImage, maskImageFloat, inputMin, inputMax, outfile);


    ////////////////////////////////////////////////////////////////////
    ///// Compute Feature for the Input images using ScalarImageToRunLengthFeaturesFilter /////
    ////////////////////////////////////////////////////////////////////
    if (isIntensity && isGLRM) {
        MaskImageType::PointType pointMin;
        MaskImageType::PointType pointMax;

        maskImage->TransformIndexToPhysicalPoint(roiStart, pointMin);
        maskImage->TransformIndexToPhysicalPoint(roiEnd, pointMax);

        double distMax = pointMin.EuclideanDistanceTo(pointMax);

        CalcGlrmFeatures(inputImage, maskImageFloat, inputMin, inputMax, distMax, outfile);
    }

    #if 0
    binaryWriter->SetFileName("test1.dcm");
    binaryWriter->SetInput(maskImage);
    binaryWriter->Update();
    #endif

    ////////////////////////////////////////////////////////////////////
    ///// Compute Glucosity Feature /////
    ////////////////////////////////////////////////////////////////////
    #if 0
    {
        outfile << "Glucosity = " << (labelObjectS->GetPhysicalSize()) * (labelObjectS->GetMean()) << endl;
    }
    #endif

    outfile.close();

    return EXIT_SUCCESS;
} // main
