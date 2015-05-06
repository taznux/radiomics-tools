
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkCommand.h>


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

// For file operation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// For threshold
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

// For island removing

#include <itkBinaryShapeKeepNObjectsImageFilter.h>


using namespace std;

void decodeCharactheristicGLCM(int x, char name[])
{
    switch (x)
    {
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

void decodeCharactheristicGLRM(int x, char name[])
{
    switch (x)
    {
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
}



int main( int argc, char *argv[] )
{


    if ( argc < 3 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " InputImage MaskImage FeatureFile" << std::endl;
        return EXIT_FAILURE;
    }

    // Definition of Data type
    const    unsigned int Dimension = 3;
    typedef  signed short InputPixelType;
    typedef   float  InternalPixelType;
    typedef signed short OutputPixelType;
    const unsigned int OutputDimension = 2;

    typedef unsigned char MaskPixelType;
    typedef unsigned char LabelPixelType;
    const LabelPixelType labelValue = 1;
    bool fullyConnected =  1;

    typedef itk::Image< InputPixelType, Dimension >  InputImageType;
    typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
    typedef itk::Image< OutputPixelType, OutputDimension > OutputImageType;

    typedef itk::Image<MaskPixelType, Dimension> MaskImageType;


    // Definition of IO
    typedef itk::ImageFileReader< InputImageType > InputImageReaderType;
    typedef itk::ImageFileReader< MaskImageType > MaskImageTypeReaderType;

    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    typedef itk::ImageFileWriter< InternalImageType >  WriterTypeFloat;

    typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;


    // Definition of Properties
    typedef InputImageType::SpacingType    SpacingType;
    typedef InputImageType::PointType      OriginType;
    typedef InputImageType::RegionType     RegionType;
    typedef InputImageType::SizeType       SizeType;
    typedef InputImageType::IndexType      IndexType;


    // Definition of Feature Extraction
    typedef itk::MinimumMaximumImageCalculator< InputImageType > MaxMinFilterType;




    // For Island removing

    //   typedef unsigned long                         LabelType;
    // typedef itk::ShapeLabelObject<LabelType, Dimension> LabelObjectType;
    // typedef itk::LabelMap<LabelObjectType>        LabelMapType;
    //  typedef itk::LabelImageToLabelMapFilter   <MaskImageType, LabelMapType> ConverterType


    /////////To Read CT images ////////////////////
    InputImageReaderType::Pointer InImageReader = InputImageReaderType::New();
    InImageReader->SetFileName(argv[1]);
    InImageReader->Update();

    InputImageType::Pointer InputImage = InImageReader->GetOutput();

    // To check coordinate of the InputImage Image
    SpacingType InputImageSpacing = InputImage->GetSpacing();
    OriginType  InputImageOrigin  = InputImage->GetOrigin();
    RegionType  InputImageRegion  = InputImage->GetLargestPossibleRegion();
    SizeType    InputImageSize    = InputImageRegion.GetSize();

    std::cout << "InputImageSpacing = " << InputImageSpacing << std::endl;
    std::cout << "InputImageOrigin = " << InputImageOrigin << std::endl;
    std::cout << "InputImageSize = " << InputImageSize << std::endl << std::endl << std::endl;


    /////////////////////////////


    MaskImageType::Pointer mask = MaskImageType::New();


    // Using the input mask
    if ( argc > 3 )
    {

        std::cout << "Get the mask from Input" << std::endl;
        ///////// To Read the Mask images ////////////////////
        MaskImageTypeReaderType::Pointer  MaskImageReader = MaskImageTypeReaderType::New();
        
        MaskImageReader->SetFileName(argv[2]);
        MaskImageReader->Update();

        mask = MaskImageReader->GetOutput();

    }


    InputImageType::IndexType index = image->GetLargestPossibleRegion().GetIndex();
    InputImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();

    OutputImageType::IndexType roiStart;
    OutputImageType::IndexType roiEnd;


    roiStart[0] = 0; roiStart[1] = 0; roiStart[2] = 0;
    roiEnd[0] = 0; roiEnd[1] = 0; roiEnd[2] = 0;

    bool foundLabel = false;

    itk::ImageRegionIteratorWithIndex< OutImageType > label(mask, mask->GetBufferedRegion() );

    for (label.GoToBegin(); !label.IsAtEnd();++label)
    {
        OutputImageType::PixelType color = label.Get();
        if (color > 0)
        {
            OutputImageType::IndexType idx = label.GetIndex();
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (!foundLabel)
                {
                    roiStart[i] = idx[i];
                    roiEnd[i] = idx[i];
                }
                else
                {
                    if (idx[i] <= roiStart[i])
                    {
                        roiStart[i] = idx[i];
                    }
                    if (idx[i] >= roiEnd[i])
                    {
                        roiEnd[i] = idx[i];
                    }
                }
            }
            foundLabel = true;
        }
    }

    InputImageType::IndexType istart;
    InputImageType::SizeType isize;

    OutputImageType::IndexType ostart;
    OutputImageType::SizeType osize;

    OutputImageType::IndexType roiCenter;

    for (unsigned n = 0; n < ndims; n++)
    {
        istart[n] = roiStart[n];
        isize[n] = roiEnd[n] - roiStart[n];

        ostart[n] = istart[n];
        osize[n] = isize[n];

        roiCenter[n] = static_cast<OutputPixelType>(round(isize[n] / 2));
    }

    InImageType::RegionType iRegion;
    iRegion.SetSize( isize );
    iRegion.SetIndex( istart );

    typedef itk::RegionOfInterestImageFilter< InImageType, InImageType > iFilterType;
    iFilterType::Pointer fInput = iFilterType::New();
    fInput->SetRegionOfInterest( iRegion );

    fInput->SetInput( InputImage );
    fInput->Update();


    OutImageType::RegionType oRegion;
    oRegion.SetSize(osize);
    oRegion.SetIndex(ostart);

    typedef itk::RegionOfInterestImageFilter< OutImageType, OutImageType > oFilterType;
    oFilterType::Pointer fOutput = oFilterType::New();
    fOutput->SetRegionOfInterest( oRegion );

    fOutput->SetInput( mask );
    fOutput->Update();

    InputImage = fInput->GetOutput();
    mask = fOutput->GetOutput();


    // To check coordinate For mask Image
    SpacingType maskImageSpacing = mask->GetSpacing();
    OriginType  maskImageOrigin  = mask->GetOrigin();
    RegionType  maskImageRegion  = mask->GetLargestPossibleRegion();
    SizeType    maskImageRegionSize = maskImageRegion.GetSize();

    std::cout << "maskImageSpacing = " << maskImageSpacing << std::endl;
    std::cout << "maskImageOrigin = " << maskImageOrigin << std::endl;
    // std::cout << "CTfixedRegion = " << fixedRegion<< std::endl;
    std::cout << "maskImageRegionSize = " << maskImageRegionSize << std::endl << std::endl << std::endl;

    std::string InputImageName = argv[1];
    std::string MaskImageName = argv[2];
    std::string FeatureFileName = argv[3];


    std::cout << "InputImageName = " << InputImageName << std::endl;
    std::cout << "MaskImageName = " << MaskImageName << std::endl;
    std::cout << "FeatureFileName = " << FeatureFileName << std::endl;


    std::ofstream outfile((FeatureFileName).c_str(), ios::out);

    if (!outfile)
    {
        cerr << "Can't open output file " << FeatureFileName << endl;
        exit(1);
    }


    // Producing the binary masks

    typedef itk::BinaryThresholdImageFilter< MaskImageType, MaskImageType>  Thresholdfiltertype;
    Thresholdfiltertype::Pointer      Thresholdfilter =  Thresholdfiltertype::New();


    Thresholdfilter->SetInput( mask  );
    Thresholdfilter->SetOutsideValue(0);
    Thresholdfilter->SetInsideValue(1);
    Thresholdfilter->SetLowerThreshold(1);
    Thresholdfilter->Update();



    // For Island Removing


    // typedef itk::LabelShapeKeepNObjectsImageFilter< MaskImageType > LabelOpeningType;

    // typedef itk::BinaryShapeKeepNObjectsImageFilter<MaskImageType > LabelOpeningType;

    // LabelOpeningType::Pointer opening = LabelOpeningType::New();
    // opening->SetInput( Thresholdfilter->GetOutput() );
    // opening->SetBackgroundValue( 0 );
    // opening->SetForegroundValue( 1 );
    // opening->SetNumberOfObjects( 1 );
    // //opening->SetReverseOrdering( false );
    // opening->SetAttribute( 100); //  extract the largest object  NUMBER_OF_PIXELS
    // opening->Update();


    /////////////////////////////////////////////////////////////////
    // Compute Intensity Features for three SUV images using itkLabelStatisticsImageFilter /////////
    ///////////////////////////////////////////////////////////////


    typedef itk::LabelStatisticsImageFilter< InternalImageType, MaskImageType > LabelStatisticsFilterType;

    LabelStatisticsFilterType::Pointer LabelStatisticsfilter = LabelStatisticsFilterType::New();
    LabelStatisticsfilter->UseHistogramsOn();

    LabelStatisticsfilter->SetInput (InputImageFloat);
    LabelStatisticsfilter->SetLabelInput (mask);

    LabelStatisticsfilter->Update();
    std::cout << "InputImage features from the LabelStatisticsImageFilter:" << std::endl;
    //outfile << std::endl<< "InputImage features from the LabelStatisticsImageFilter:" << std::endl;

    std::cout << "Minimum   = " << LabelStatisticsfilter->GetMinimum( labelValue ) << std::endl;
    std::cout << "Maximum   = " << LabelStatisticsfilter->GetMaximum( labelValue )      << std::endl;
    std::cout << "Median    = " << LabelStatisticsfilter->GetMedian( labelValue )   << std::endl;
    std::cout << "Mean      = " << LabelStatisticsfilter->GetMean( labelValue )     << std::endl;
    std::cout << "Sigma     = " << LabelStatisticsfilter->GetSigma( labelValue )    << std::endl;
    std::cout << "Variance  = " << LabelStatisticsfilter->GetVariance( labelValue ) << std::endl;
    std::cout << "Sum       = " << LabelStatisticsfilter->GetSum( labelValue )      << std::endl;
    //std::cout << "Bounding box = " << LabelStatisticsfilter->GetBoundingBox( labelValue )  << std::endl;
    std::cout << "Region    = " << LabelStatisticsfilter->GetRegion( labelValue )   << std::endl;
    std::cout << "HistogramPointer    = " << LabelStatisticsfilter->GetHistogram( labelValue )   << std::endl;

    outfile << " Minimum= " << LabelStatisticsfilter->GetMinimum( labelValue ) << std::endl;
    outfile << " Maximum= " << LabelStatisticsfilter->GetMaximum( labelValue )      << std::endl;
    //outfile<< "Median    = " << LabelStatisticsfilter->GetMedian( labelValue )   << std::endl;
    outfile << " Mean= " << LabelStatisticsfilter->GetMean( labelValue )     << std::endl;
    outfile << " Sigma= " << LabelStatisticsfilter->GetSigma( labelValue )    << std::endl;
    outfile << " Variance= " << LabelStatisticsfilter->GetVariance( labelValue ) << std::endl;
    outfile << " Sum= " << LabelStatisticsfilter->GetSum( labelValue )      << std::endl;
    //outfile<< "Bounding box = " << LabelStatisticsfilter->GetBoundingBox( labelValue )  << std::endl;
    //outfile<< "Region    = " << LabelStatisticsfilter->GetRegion( labelValue )   << std::endl << std::endl << std::endl;
    std::cout << "HistogramPointer    = " << LabelStatisticsfilter->GetHistogram( labelValue )   << std::endl;


    float InputMin = LabelStatisticsfilter->GetMinimum( labelValue );
    float InputMax = LabelStatisticsfilter->GetMaximum( labelValue );

    /////////////////////////////////////////////////////////////////////
    // Compute feature three SUV Geometry feature using LabelGeometryImageFilter  ////
    /////////////////////////////////////////////////////////////////////


    typedef itk::LabelGeometryImageFilter< MaskImageType, InternalImageType > LabelGeometryType;
    LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();

    // Set up a connected components filter to label the binary objects.
    typedef itk::ConnectedComponentImageFilter< MaskImageType, MaskImageType > ConnectedComponentType;
    ConnectedComponentType::Pointer connectedComponentFilter = ConnectedComponentType::New();
    connectedComponentFilter->SetInput(opening->GetOutput());

    // Relabel the components in order of size.
    typedef itk::RelabelComponentImageFilter< MaskImageType, MaskImageType > RelabelType;
    RelabelType::Pointer relabeler = RelabelType::New();
    relabeler->SetInput( connectedComponentFilter->GetOutput() );
    labelGeometryFilter->SetInput( relabeler->GetOutput() );
    //labelGeometryFilter->SetInput(opening->GetOutput());

    // These generate optional outputs.
    labelGeometryFilter->CalculatePixelIndicesOn();
    labelGeometryFilter->CalculateOrientedBoundingBoxOn();
    labelGeometryFilter->CalculateOrientedLabelRegionsOn();
    labelGeometryFilter->CalculateOrientedIntensityRegionsOn();

    labelGeometryFilter->ReleaseDataFlagOn();

    labelGeometryFilter->SetIntensityInput(InputImageFloat); //
    labelGeometryFilter->Update();
    std::cout << "Fixed SUV image features from the LabelGeometryImageFilter:" << std::endl;
    //outfile << "InputImage features from the LabelGeometryImageFilter:" << std::endl;

    outfile << " Volume= " << labelGeometryFilter->GetVolume(labelValue) << std::endl;
    outfile << " IntegratedIntensity= " << labelGeometryFilter->GetIntegratedIntensity(labelValue) << std::endl;
    //outfile << "\tCentroid: " << labelGeometryFilter->GetCentroid(labelValue) << std::endl;
    //outfile << "\tWeighted Centroid: " << labelGeometryFilter->GetWeightedCentroid(labelValue) << std::endl;
    //outfile << "\tAxes Length: " << labelGeometryFilter->GetAxesLength(labelValue) << std::endl;
    outfile << " MajorAxisLength= " << labelGeometryFilter->GetMajorAxisLength(labelValue) << std::endl;
    outfile << " MinorAxisLength= " << labelGeometryFilter->GetMinorAxisLength(labelValue) << std::endl;
    outfile << " Eccentricity= " << labelGeometryFilter->GetEccentricity(labelValue) << std::endl;
    outfile << " Elongation= " << labelGeometryFilter->GetElongation(labelValue) << std::endl;
    outfile << " Orientation= " << labelGeometryFilter->GetOrientation(labelValue) << std::endl;
    //outfile << "\tBounding box: " << labelGeometryFilter->GetBoundingBox(labelValue) << std::endl;

    outfile << " BoundingBoxVolume= " << labelGeometryFilter->GetBoundingBoxVolume(labelValue) << std::endl;
    //outfile << "\tBounding box size: " << labelGeometryFilter->GetBoundingBoxSize(labelValue) << std::endl;
    outfile << " OrientedBoundingBoxVolume= " << labelGeometryFilter->GetOrientedBoundingBoxVolume(labelValue) << std::endl;
    //outfile << "\tOriented Bounding box size: " << labelGeometryFilter->GetOrientedBoundingBoxSize(labelValue) << std::endl;

    //outfile << "\tEigenvalues: " << labelGeometryFilter->GetEigenvalues(labelValue) << std::endl;
    //outfile << "\tEigenvectors: " << labelGeometryFilter->GetEigenvectors(labelValue) << std::endl;
    //outfile << "\tRotationMatrix: " << labelGeometryFilter->GetRotationMatrix(labelValue) << std::endl<< std::endl << std::endl;
    //outfile << "\n\n";

    std::cout << "Compute shape and intensity feature uisng LabelMap..." << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Compute shape and intensity feature uisng LabelMap   ////////////////////
    ////////////////////////////////////////////////////////////////////////////

    /////////////////////////////

    typedef itk::StatisticsLabelObject< unsigned char, Dimension > LabelObjectTypeS;
    typedef itk::LabelMap<LabelObjectTypeS >            LabelMapTypeS;

    //converting binary image to Statistics label map
    typedef itk::BinaryImageToStatisticsLabelMapFilter< MaskImageType, InternalImageType, LabelMapTypeS > I2LSType;
    I2LSType::Pointer i2ls = I2LSType::New();
    // itk::SimpleFilterWatcher watcher1( i2ls );
    i2ls->SetFullyConnected( fullyConnected );

    i2ls->FullyConnectedOn();
    int inputForegroundValue = 1;   bool computeFeretDiameter =  1;     unsigned int outputBackgroundValue = 0;
    bool computePerimeter =  1;     bool computeHistogram =  1;     unsigned int numberOfBins = 32;
    i2ls->SetInputForegroundValue( inputForegroundValue );
    i2ls->SetOutputBackgroundValue( outputBackgroundValue );
    i2ls->SetComputeFeretDiameter( computeFeretDiameter );
    i2ls->ComputeFeretDiameterOff();
    i2ls->ComputeFeretDiameterOn();
    i2ls->SetComputePerimeter( computePerimeter );
    i2ls->ComputePerimeterOn();
    i2ls->SetComputeHistogram( computeHistogram );
    i2ls->ComputeHistogramOn();
    i2ls->SetNumberOfBins( numberOfBins );


    LabelMapTypeS::Pointer labelMapS;
    LabelObjectTypeS   *labelObjectS;
    i2ls->SetFeatureImage( InputImageFloat );
    i2ls->SetInput( opening->GetOutput());

    i2ls->Update();
    std::cout << "InputImage features from the BinaryImageToStatisticsLabelMapFilter:" << std::endl;
    //outfile << "InputImage features from the BinaryImageToStatisticsLabelMapFilter:" << std::endl;


    labelMapS = i2ls->GetOutput();
    labelObjectS = labelMapS->GetLabelObject(labelValue);


    //    std::cout << "\tBoundingBox: " << labelObjectS->GetBoundingBox() << std::endl;
    std::cout << "\tCenterOfGravity: " << labelObjectS->GetCenterOfGravity() << std::endl;
    std::cout << "\tCentroid: " << labelObjectS->GetCentroid() << std::endl;
    std::cout << "\tElongation: " << labelObjectS->GetElongation() << std::endl;
    //  std::cout << "\tEquivalentEllipsoidDiameter: " << labelObjectS->GetEquivalentEllipsoidDiameter() << std::endl;
    //  std::cout << "\tEquivalentSphericalPerimeter: " << labelObjectS->GetEquivalentSphericalPerimeter() << std::endl;
    //       std::cout << "\tEquivalentSphericalRadius: " << labelObjectS->GetEquivalentSphericalRadius() << std::endl;
    std::cout << "\tFeretDiameter: " << labelObjectS->GetFeretDiameter() << std::endl;
    std::cout << "\tHistogram: " << labelObjectS->GetHistogram() << std::endl;
    //  std::cout << "\tIndex: " << labelObjectS->GetIndex() << std::endl;
    std::cout << "\tKurtosis: " << labelObjectS->GetKurtosis() << std::endl;
    std::cout << "\tFlatness: " << labelObjectS->GetFlatness() << std::endl;
    std::cout << "\tLabel: " << labelObjectS->GetLabel() << std::endl;
    //  std::cout << "\tLine: " << labelObjectS->GetLine() << std::endl;
    //        std::cout << "\tLineContainer: " << labelObjectS->GetLineContainer() << std::endl;
    std::cout << "\tMaximum: " << labelObjectS->GetMaximum() << std::endl;
    std::cout << "\tMaximumIndex: " << labelObjectS->GetMaximumIndex() << std::endl;
    std::cout << "\tMean: " << labelObjectS->GetMean() << std::endl;
    std::cout << "\tMedian: " << labelObjectS->GetMedian() << std::endl;
    std::cout << "\tMinimum: " << labelObjectS->GetMinimum() << std::endl;
    std::cout << "\tMinimumIndex: " << labelObjectS->GetMinimumIndex() << std::endl;
    std::cout << "\tNameOfClass: " << labelObjectS->GetNameOfClass() << std::endl;
    std::cout << "\tNumberOfLines: " << labelObjectS->GetNumberOfLines() << std::endl;
    //        std::cout << "\tNumberOfPixels: " << labelObjectS->GetNumberOfPixels() << std::endl;
    //  std::cout << "\tNumberOfPixelsOnBorder: " << labelObjectS->GetNumberOfPixelsOnBorder() << std::endl;
    std::cout << "\tPerimeter: " << labelObjectS->GetPerimeter() << std::endl;
    //  std::cout << "\tPerimeterOnBorder: " << labelObjectS->GetPerimeterOnBorder() << std::endl;
    //        std::cout << "\tPerimeterOnBorderRatio: " << labelObjectS->GetPerimeterOnBorderRatio() << std::endl;
    std::cout << "\tPhysicalSize: " << labelObjectS->GetPhysicalSize() << std::endl;
    std::cout << "\tPrincipalAxes: " << labelObjectS->GetPrincipalAxes() << std::endl;
    std::cout << "\tPrincipalMoments: " << labelObjectS->GetPrincipalMoments() << std::endl;
    std::cout << "\tReferenceCount: " << labelObjectS->GetReferenceCount() << std::endl;
    std::cout << "\tRegionElongation: " << labelObjectS->GetElongation() << std::endl; // RegionElongation is disappeared
    std::cout << "\tRoundness: " << labelObjectS->GetRoundness() << std::endl;
    std::cout << "\tSizeRegionRatio: " << 0 /*labelObjectS->GetSizeRegionRatio()*/ << std::endl; // This property is no longer provide
    std::cout << "\tSkewness: " << labelObjectS->GetSkewness() << std::endl;
    //  std::cout << "\tStandardDeviation: " << labelObjectS->GetStandardDeviation() << std::endl;
    std::cout << "\tSum: " << labelObjectS->GetSum() << std::endl;
    std::cout << "\tVariance: " << labelObjectS->GetVariance() << std::endl;
    //  std::cout << "\tWeightedElongation: " << labelObjectS->GetWeightedElongation() << std::endl;
    //        std::cout << "\tWeightedFlatness: " << labelObjectS->GetWeightedFlatness() << std::endl;
    //  std::cout << "\tWeightedPrincipalAxes: " << labelObjectS->GetWeightedPrincipalAxes() << std::endl;
    //        std::cout << "\tWeightedPrincipalMoments: " << labelObjectS->GetWeightedPrincipalMoments() << std::endl;
    //    outfile << "\tBoundingBox: " << labelObjectS->GetBoundingBox() << std::endl;
    //outfile << "\tCenterOfGravity: " << labelObjectS->GetCenterOfGravity() << std::endl;
    //outfile << "\tCentroid: " << labelObjectS->GetCentroid() << std::endl;
    outfile << " Elongation= " << labelObjectS->GetElongation() << std::endl;
    //  outfile << "\tEquivalentEllipsoidDiameter: " << labelObjectS->GetEquivalentEllipsoidDiameter() << std::endl;
    //  outfile << "\tEquivalentSphericalPerimeter: " << labelObjectS->GetEquivalentSphericalPerimeter() << std::endl;
    //       outfile<< "\tEquivalentSphericalRadius: " << labelObjectS->GetEquivalentSphericalRadius() << std::endl;
    outfile << " FeretDiameter= " << labelObjectS->GetFeretDiameter() << std::endl;
    //outfile << "\tHistogram: " << labelObjectS->GetHistogram() << std::endl;
    //  outfile << "\tIndex: " << labelObjectS->GetIndex() << std::endl;
    outfile << " Kurtosis= " << labelObjectS->GetKurtosis() << std::endl;
    outfile << " Flatness= " << labelObjectS->GetFlatness() << std::endl;
    //outfile << "\tLabel: " << labelObjectS->GetLabel() << std::endl;
    //  outfile << "\tLine: " << labelObjectS->GetLine() << std::endl;
    //        outfile << "\tLineContainer: " << labelObjectS->GetLineContainer() << std::endl;
    //outfile << "\tMaximum: " << labelObjectS->GetMaximum() << std::endl;
    //outfile << "\tMaximumIndex: " << labelObjectS->GetMaximumIndex() << std::endl;
    //outfile << "\tMean: " << labelObjectS->GetMean() << std::endl;
    outfile << " Median= " << labelObjectS->GetMedian() << std::endl;
    //outfile << "\tMinimum: " << labelObjectS->GetMinimum() << std::endl;
    //outfile << "\tMinimumIndex: " << labelObjectS->GetMinimumIndex() << std::endl;
    //outfile << "\tNameOfClass: " << labelObjectS->GetNameOfClass() << std::endl;
    outfile << " NumberOfLines= " << labelObjectS->GetNumberOfLines() << std::endl;
    //        outfile << "\tNumberOfPixels: " << labelObjectS->GetNumberOfPixels() << std::endl;
    //  outfile << "\tNumberOfPixelsOnBorder: " << labelObjectS->GetNumberOfPixelsOnBorder() << std::endl;
    outfile << " Perimeter= " << labelObjectS->GetPerimeter() << std::endl;
    //  outfile << "\tPerimeterOnBorder: " << labelObjectS->GetPerimeterOnBorder() << std::endl;
    //        outfile << "\tPerimeterOnBorderRatio: " << labelObjectS->GetPerimeterOnBorderRatio() << std::endl;
    outfile << " PhysicalSize= " << labelObjectS->GetPhysicalSize() << std::endl;
    //outfile << "\tPrincipalAxes: " << labelObjectS->GetPrincipalAxes() << std::endl;
    //  outfile << "\tPrincipalMoments: " << labelObjectS->GetPrincipalMoments() << std::endl;
    //  outfile << "\tReferenceCount: " << labelObjectS->GetReferenceCount() << std::endl;
    outfile << " RegionElongation= " << labelObjectS->GetElongation() << std::endl; // deleted
    outfile << " Roundness= " << labelObjectS->GetRoundness() << std::endl;
    outfile << " SizeRegionRatio= " << 0 /*labelObjectS->GetSizeRegionRatio()*/ << std::endl; // deleted
    outfile << " Skewness= " << labelObjectS->GetSkewness() << std::endl;
    //  outfile << "\tStandardDeviation: " << labelObjectS->GetStandardDeviation() << std::endl;
    //outfile << "\tSum: " << labelObjectS->GetSum() << std::endl;

    //outfile << "\tVariance: " << labelObjectS->GetVariance() << std::endl<< std::endl<< std::endl;
    //  outfile << "\tWeightedElongation: " << labelObjectS->GetWeightedElongation() << std::endl;
    //        outfile << "\tWeightedFlatness: " << labelObjectS->GetWeightedFlatness() << std::endl;
    //  outfile << "\tWeightedPrincipalAxes: " << labelObjectS->GetWeightedPrincipalAxes() << std::endl;
    //        outfile << "\tWeightedPrincipalMoments: " << labelObjectS->GetWeightedPrincipalMoments() << std::endl;




    ////////////////////////////////////////////////////////////////////
    ///// Compute Feature for the Input images using ScalarImageToTextureFeaturesFilter /////
    ////////////////////////////////////////////////////////////////////



    typedef itk::CastImageFilter<InputImageType, MaskImageType > SUVCastFilterType;
    SUVCastFilterType::Pointer SUVCaster   = SUVCastFilterType::New();

    typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<MaskImageType> TextureFilterType;
    TextureFilterType::Pointer glcm = TextureFilterType::New();

    typedef TextureFilterType::TextureFeaturesFilterType   TextureFeaturesFilterType;

    TextureFilterType::FeatureNameVectorPointer requestedFeatures =
        TextureFilterType::FeatureNameVector::New();

    requestedFeatures->push_back(TextureFeaturesFilterType::Energy);
    requestedFeatures->push_back(TextureFeaturesFilterType::Entropy);
    requestedFeatures->push_back(TextureFeaturesFilterType::Correlation);
    requestedFeatures->push_back(TextureFeaturesFilterType::InverseDifferenceMoment);
    requestedFeatures->push_back(TextureFeaturesFilterType::Inertia);
    requestedFeatures->push_back(TextureFeaturesFilterType::ClusterShade);
    requestedFeatures->push_back(TextureFeaturesFilterType::ClusterProminence);
    requestedFeatures->push_back(TextureFeaturesFilterType::HaralickCorrelation);

    glcm->SetRequestedFeatures(requestedFeatures);
    //glcm->SetNormalizeOn();
    //glcm->Print( std::cout );



    SUVCaster->SetInput( InputImage );
    SUVCaster->Update();


    // Normalize image by setting min and max
    /*typedef itk::RescaleIntensityImageFilter< InputImageType, MaskImageType> RescaleFilterType;
     RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
     rescaleFilter->SetInput(InputImage);
     rescaleFilter->SetOutputMinimum(0);
     rescaleFilter->SetOutputMaximum(63);
    */





    InternalImageType::Pointer RevisedImage   = InternalImageType::New();
    RevisedImage->SetRegions( InputImage->GetRequestedRegion() );
    RevisedImage->CopyInformation( InputImage );
    RevisedImage->Allocate();

    typedef itk::ImageRegionIterator< InternalImageType>       RevisedImageIteratorFloatType;
    RevisedImageIteratorFloatType RevisedImageIt( RevisedImage, RevisedImage->GetRequestedRegion() );


    InputImageFloatIt.GoToBegin();
    RevisedImageIt.GoToBegin();

    while ( !InputImageFloatIt.IsAtEnd() )
    {
        RevisedImageIt.Set( InputImageFloatIt.Get() );
        if (InputImageFloatIt.Get() > InputMax)
        {

            // std::cout << "The current intensity is "<< RevisedImageIt.Get() << std::endl;
            // std::cout << "The maximum is "<< InputMax << std::endl;
            RevisedImageIt.Set(InputMax);
            //std::cout << "The revised intensity is "<< RevisedImageIt.Get() << std::endl;

            //int jmax;
            //scanf("%d", &jmax);


        }
        if (InputImageFloatIt.Get() < InputMin)
        {

            //std::cout << "The current intensity is "<< RevisedImageIt.Get() << std::endl;
            //std::cout << "The minimum is "<< InputMin << std::endl;
            RevisedImageIt.Set(InputMin);
            //std::cout << "The revised intensity is "<< RevisedImageIt.Get() << std::endl;

            //int jmin;
            //scanf("%d", &jmin);
        }
        ++InputImageFloatIt;
        ++RevisedImageIt;
    }





    typedef itk::RescaleIntensityImageFilter< InternalImageType, MaskImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(RevisedImage);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(63);





    // glcm->SetInput(SUVCaster->GetOutput());

    glcm->SetInput(rescaleFilter->GetOutput());
    // glcm->SetInput(opening->GetOutput());
    glcm->SetMaskImage(opening->GetOutput());
    // glcm->FastCalculationsOn();
    glcm->SetNumberOfBinsPerAxis(64);

    glcm->FastCalculationsOff();
    glcm->Update();
    std::cout << "The InputImage features from the ScalarImageToTextureFeaturesFilter:" << std::endl;
    // outfile<< std::endl<< std::endl << "The InputImage features from the ScalarImageToTextureFeaturesFilter:" << std::endl;


    int size = glcm->GetFeatureMeans()->Size();
    for (int x = 0; x < size; x++)
    {
        char charactheristic[100];
        decodeCharactheristicGLCM(glcm->GetRequestedFeatures()->GetElement(x), charactheristic);

        std::cout << " Mean of " << charactheristic << "= " <<  glcm->GetFeatureMeans()->GetElement(x) << std::endl;
        //  std::cout << " StandardDeviation of"<< charactheristic << "= " <<  glcm->GetFeatureStandardDeviations()->GetElement(x) << std::endl;

        outfile << " Mean of" << charactheristic << "= " << glcm->GetFeatureMeans()->GetElement(x) << std::endl;
        //  outfile << " StandardDeviationof"<< charactheristic << "= " << glcm->GetFeatureStandardDeviations()->GetElement(x) << std::endl;
    }


    {
        typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter<MaskImageType> RunLengthFilterType;
        RunLengthFilterType::Pointer glrm = RunLengthFilterType::New();

        typedef RunLengthFilterType::RunLengthFeaturesFilterType   RunLengthFeaturesFilterType;

        RunLengthFilterType::FeatureNameVectorPointer requestedFeatures =
            RunLengthFilterType::FeatureNameVector::New();

        requestedFeatures->push_back( RunLengthFeaturesFilterType::ShortRunEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::LongRunEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::GreyLevelNonuniformity );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::RunLengthNonuniformity );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::LowGreyLevelRunEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::HighGreyLevelRunEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::ShortRunLowGreyLevelEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::ShortRunHighGreyLevelEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::LongRunLowGreyLevelEmphasis );
        requestedFeatures->push_back( RunLengthFeaturesFilterType::LongRunHighGreyLevelEmphasis );

        glrm->SetRequestedFeatures(requestedFeatures);


        // glrm->SetInput(SUVCaster->GetOutput());

        glrm->SetInput(rescaleFilter->GetOutput());
        // glrm->SetInput(opening->GetOutput());
        glrm->SetMaskImage(opening->GetOutput());
        // glrm->FastCalculationsOn();
        glrm->SetNumberOfBinsPerAxis(64);

        glrm->FastCalculationsOff();
        glrm->Update();
        std::cout << "The InputImage features from the ScalarImageToRunLengthFeaturesFilter:" << std::endl;
        // outfile<< std::endl<< std::endl << "The InputImage features from the ScalarImageToTextureFeaturesFilter:" << std::endl;


        int size = glrm->GetFeatureMeans()->Size();
        for (int x = 0; x < size; x++)
        {
            char charactheristic[100];
            decodeCharactheristicGLRM(glrm->GetRequestedFeatures()->GetElement(x), charactheristic);

            std::cout << " Mean of " << charactheristic << "= " <<  glrm->GetFeatureMeans()->GetElement(x) << std::endl;
            std::cout << " StandardDeviation of"<< charactheristic << "= " <<  glrm->GetFeatureStandardDeviations()->GetElement(x) << std::endl;

            outfile << " Mean of" << charactheristic << "= " << glrm->GetFeatureMeans()->GetElement(x) << std::endl;
            outfile << " StandardDeviationof"<< charactheristic << "= " << glrm->GetFeatureStandardDeviations()->GetElement(x) << std::endl;
        }
    }


    ////////////////////////////////////////////////////////////////////
    ///// Compute Glucosity Feature /////
    ////////////////////////////////////////////////////////////////////





    outfile << " Glucosity= " << (labelObjectS->GetPhysicalSize())*(labelObjectS->GetMean())  << std::endl;

    outfile.close();
    return EXIT_SUCCESS;
}
