#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkImageFileWriter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkOtsuMultipleThresholdsCalculator.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"

#include "itkConnectedComponentImageFilter.h"

#include "itkGrowCutSegmentationImageFilter.h"
#include "itkFastGrowCutSegmentationImageFilter.h"


using namespace std;

// To print the progress
class ShowProgressObject
{
public:
    ShowProgressObject( itk::ProcessObject *o )
    {
        m_Process = o;
    }
    void ShowProgress()
    {
        std::cout << "\rProgress: "
                  << static_cast<unsigned int>( 100.0 * m_Process->GetProgress() ) << "%" << endl;
    }
    itk::ProcessObject::Pointer m_Process;
};

int main( int argc, char *argv[] )
{
    if ( argc < 3 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << "InputImageFile,SeedImageFile,OutputImageFile" << std::endl;
        return EXIT_FAILURE;
    }

    // Definition of Data type
    typedef signed short InputPixelType;
    typedef unsigned char OutputPixelType;
    typedef float InternalPixelType;


    typedef itk::Image<InputPixelType, 3> InImageType;
    typedef itk::Image<OutputPixelType, 3> OutImageType;
    typedef itk::Image<OutputPixelType, 2> OutDCMImageType;
    typedef itk::Image<InternalPixelType, 3> WeightImageType;

    // Definition of IO
    typedef itk::ImageSeriesReader< InImageType > InImageReaderType;
    typedef itk::ImageSeriesReader< OutImageType > LabImageReaderType;
    typedef itk::ImageSeriesWriter< OutImageType, OutDCMImageType > OutImageWriterType;
    typedef itk::ImageFileWriter< OutImageType> SingleOutImageWriterType;


    typedef itk::GDCMImageIO    ImageIOType;
    ImageIOType::Pointer gdcmIO = ImageIOType::New();

    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    typedef std::vector< std::string > FileNamesContainer;
    NamesGeneratorType::Pointer ImageNameGenerator = NamesGeneratorType::New();


    double objectSize = 20;
    double priorSegmentStrength = 1.0;
    double contrastNoiseRatio = 0.003;

    unsigned char foreground = 255;
    unsigned char background = 0;


    /////////To Read CT images ////////////////////
    InImageReaderType::Pointer InImageReader = InImageReaderType::New();
    InImageReader->SetImageIO( gdcmIO );
    FileNamesContainer InImageNames;
    ImageNameGenerator->SetInputDirectory( argv[1] );
    InImageNames = ImageNameGenerator->GetInputFileNames();
    InImageReader->SetFileNames(InImageNames);
    InImageReader->Update();

    /////////To Read Seed images ////////////////////
    LabImageReaderType::Pointer LabImageReader = LabImageReaderType::New();
    LabImageReader->SetImageIO( gdcmIO );
    FileNamesContainer LabImageNames;
    ImageNameGenerator->SetInputDirectory( argv[2] );
    LabImageNames = ImageNameGenerator->GetInputFileNames();
    LabImageReader->SetFileNames(LabImageNames);
    LabImageReader->Update();


    InImageType::Pointer image = InImageReader->GetOutput();
    OutImageType::Pointer labelImage = LabImageReader->GetOutput();
    OutImageType::Pointer prevSegmentedImage = LabImageReader->GetOutput(); // inPtr3


    OutImageType::Pointer outputImageROI = OutImageType::New();
    OutImageType::Pointer outputImage = OutImageType::New();
    WeightImageType::Pointer weightImage = WeightImageType::New();


    weightImage->CopyInformation(image);
    weightImage->SetBufferedRegion( image->GetBufferedRegion() );
    weightImage->Allocate();
    weightImage->FillBuffer( 0 );



    if (contrastNoiseRatio > 1.0)
    {
        contrastNoiseRatio /= 100.0;
    }

    if (priorSegmentStrength > 1.0)
    {
        priorSegmentStrength /= 100.0;
    }

    InImageType::IndexType index = image->GetLargestPossibleRegion().GetIndex();
    InImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();

    OutImageType::IndexType roiStart;
    OutImageType::IndexType roiEnd;


    roiStart[0] = 0; roiStart[1] = 0; roiStart[2] = 0;
    roiEnd[0] = 0; roiEnd[1] = 0; roiEnd[2] = 0;

    unsigned int ndims = image->GetImageDimension();

    bool foundLabel = false;

    itk::ImageRegionIterator< WeightImageType > weight(weightImage, weightImage->GetBufferedRegion() );
    itk::ImageRegionIteratorWithIndex< OutImageType > label(labelImage, labelImage->GetBufferedRegion() );

    for (weight.GoToBegin(), label.GoToBegin(); !weight.IsAtEnd();
            ++weight, ++label)
    {
        OutImageType::PixelType color = label.Get();
        if (color == 0)
        {
            weight.Set(0.0);
        }
        else
        {
            weight.Set( contrastNoiseRatio );

            OutImageType::IndexType idx = label.GetIndex();
            for (unsigned i = 0; i < ndims; i++)
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


    std::cout << " objectSize (radius) " << objectSize << std::endl;

    OutImageType::PixelType radius = static_cast< OutImageType::PixelType> (objectSize);

    for (unsigned i = 0; i < ndims; i++)
    {
        int diff = static_cast< int > (roiStart[i] - radius);
        if (diff >= index[i])
        {
            roiStart[i] -= radius;
        }
        else
        {
            roiStart[i] = index[i];
        }
        roiEnd[i] = (static_cast<unsigned int>(roiEnd[i] + radius) < size[i]) ?
                    (roiEnd[i] + radius) : size[i] - 1;

        std::cout << " roi[ " << roiStart[i] << " " << roiEnd[i] << "] " << std::endl;
    }


    typedef itk::GrowCutSegmentationImageFilter<InImageType, OutImageType> GrowCutFilterType;
    GrowCutFilterType::Pointer filter = GrowCutFilterType::New();

    //
    ShowProgressObject progressWatch(filter);
    itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
    command = itk::SimpleMemberCommand<ShowProgressObject>::New();
    command->SetCallbackFunction(&progressWatch,
                                 &ShowProgressObject::ShowProgress);

    filter->AddObserver(itk::ProgressEvent(), command);

    InImageType::IndexType istart;
    InImageType::SizeType isize;

    OutImageType::IndexType ostart;
    OutImageType::SizeType osize;

    WeightImageType::IndexType wstart;
    WeightImageType::SizeType wsize;

    OutImageType::IndexType roiCenter;

    for (unsigned n = 0; n < ndims; n++)
    {
        istart[n] = roiStart[n];
        isize[n] = roiEnd[n] - roiStart[n];

        ostart[n] = istart[n];
        osize[n] = isize[n];

        roiCenter[n] = static_cast<OutputPixelType>(round(isize[n] / 2));

        wstart[n] = istart[n];
        wsize[n] = isize[n];
    }

    std::cout << " istart " << istart << " isize " << isize << " " << std::endl;
    std::cout << " roiCenter " << roiCenter << std::endl;

    InImageType::RegionType iRegion;
    iRegion.SetSize( isize );
    iRegion.SetIndex( istart );

    typedef itk::RegionOfInterestImageFilter< InImageType, InImageType > iFilterType;
    iFilterType::Pointer fInput = iFilterType::New();
    fInput->SetRegionOfInterest( iRegion );

    fInput->SetInput( image );
    fInput->Update();

    OutImageType::RegionType oRegion;
    oRegion.SetSize(osize);
    oRegion.SetIndex(ostart);

    typedef itk::RegionOfInterestImageFilter< OutImageType, OutImageType > oFilterType;
    oFilterType::Pointer fOutput = oFilterType::New();
    fOutput->SetRegionOfInterest( oRegion );

    fOutput->SetInput( labelImage );
    fOutput->Update();

    WeightImageType::RegionType wRegion;
    wRegion.SetSize(wsize);
    wRegion.SetIndex(wstart);

    typedef itk::RegionOfInterestImageFilter< WeightImageType, WeightImageType > wFilterType;
    wFilterType::Pointer fWeight = wFilterType::New();
    fWeight->SetRegionOfInterest( wRegion );

    fWeight->SetInput( weightImage );
    fWeight->Update();

    InImageType::Pointer inImage = InImageType::New();
    inImage = fInput->GetOutput();

    OutImageType::Pointer labImage = OutImageType::New();
    labImage = fOutput->GetOutput();

    WeightImageType::Pointer wtImage = WeightImageType::New();
    wtImage = fWeight->GetOutput();

    typedef itk::Statistics::ScalarImageToHistogramGenerator <InImageType > ScalarImageToHistogramGeneratorType;
    typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
    typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

    typedef itk::BinaryThresholdImageFilter <InImageType, OutImageType >  ThresholdFilterType;
    typedef itk::ConnectedComponentImageFilter <OutImageType, OutImageType > ConnectedComponentImageFilterType;

    {
        ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
        CalculatorType::Pointer calculator = CalculatorType::New();
        ThresholdFilterType::Pointer filter = ThresholdFilterType::New();
        ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();

        scalarImageToHistogramGenerator->SetNumberOfBins( 128 );
        calculator->SetNumberOfThresholds(1);

        const OutputPixelType outsideValue = 0;
        const OutputPixelType insideValue = 255;

        filter->SetInput(inImage);
        connected->SetInput(filter->GetOutput());

        filter->SetOutsideValue( outsideValue );
        filter->SetInsideValue(  insideValue  );

        //Connect Pipeline
        scalarImageToHistogramGenerator->SetInput( inImage);
        calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );

        scalarImageToHistogramGenerator->Compute();
        try
        {
            calculator->Compute();
        }
        catch ( itk::ExceptionObject &excp )
        {
            std::cerr << "Exception thrown " << excp << std::endl;
        }

        const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();
        CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

        InputPixelType lowerThreshold = static_cast<InputPixelType>(*itNum);
        cout << "Threshold: " << lowerThreshold << endl;
        filter->SetLowerThreshold( lowerThreshold );

        connected->Update();

        labImage = connected->GetOutput();
    }

    //foreground = labImage->GetPixel(roiCenter);
    //cout << roiCenter << "Foreground: " << static_cast<int>(foreground) << endl;
    wtImage->SetPixel(roiCenter, 1.0);
    labImage->SetPixel(roiCenter, foreground);


    // Definition of Morphological operation sturcture elements
    typedef itk::BinaryBallStructuringElement< OutputPixelType, 3 > StructuringElementType;
    typedef itk::BinaryDilateImageFilter<OutImageType, OutImageType, StructuringElementType> DilateFilterType;
    typedef itk::BinaryErodeImageFilter<OutImageType, OutImageType, StructuringElementType> ErodeFilterType;

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

    // Automated weight generation, Object and Background
    DilateFilterType::Pointer dilate3Filter = DilateFilterType::New();

    dilate3Filter->SetKernel(structuringElement3);
    dilate3Filter->SetDilateValue(foreground);
    dilate3Filter->SetInput(labImage);

    dilate3Filter->Update();


    OutImageType::Pointer objImage = OutImageType::New();
    objImage = dilate3Filter->GetOutput();

    DilateFilterType::Pointer dilate1Filter = DilateFilterType::New();
    dilate1Filter->SetKernel(structuringElement1);
    dilate1Filter->SetDilateValue(background);
    dilate1Filter->SetInput(labImage);

    dilate1Filter->Update();

    OutImageType::Pointer bgImage = OutImageType::New();
    bgImage = dilate1Filter->GetOutput();


    {
        itk::ImageRegionIterator< WeightImageType > weight(wtImage, wtImage->GetBufferedRegion() );
        itk::ImageRegionIterator< OutImageType > label(labImage, labImage->GetBufferedRegion() );
        itk::ImageRegionConstIterator< OutImageType > olabel(objImage, objImage->GetBufferedRegion() );
        itk::ImageRegionConstIterator< OutImageType > blabel(bgImage, bgImage->GetBufferedRegion() );


        for (weight.GoToBegin(), olabel.GoToBegin(), blabel.GoToBegin(), label.GoToBegin();
                !weight.IsAtEnd();
                ++weight, ++olabel, ++blabel, ++label)
        {
            OutImageType::PixelType ocolor = olabel.Get();
            OutImageType::PixelType bcolor = blabel.Get();
            if (ocolor == foreground) // object
            {
                weight.Set( priorSegmentStrength );
                label.Set (foreground);
                cout << "o";
            }
            else if (bcolor == background) // background
            {
                if (weight.Get() == 0)
                {
                    weight.Set( priorSegmentStrength );
                    cout << ".";
                }
                else
                    cout << "-";
                label.Set (background);
            }
            else // grey area
            {
                label.Set (0);
                cout << "|";
            }
        }
        cout << endl;
    }


    // Segmentation fitler parameters
    filter->SetInput( inImage );
    filter->SetLabelImage( labImage );

    filter->SetStrengthImage( wtImage );

    filter->SetSeedStrength( contrastNoiseRatio );
    filter->SetObjectRadius((unsigned int)objectSize);

    {
        ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
        ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
        DilateFilterType::Pointer dilateFilter = DilateFilterType::New();

        erodeFilter->SetKernel(structuringElement3);
        dilateFilter->SetKernel(structuringElement3);

        erodeFilter->SetErodeValue(foreground);


        erodeFilter->SetInput(filter->GetOutput());
        connected->SetInput(erodeFilter->GetOutput());
        dilateFilter->SetInput(connected->GetOutput());

        connected->Update();

        foreground = connected->GetOutput()->GetPixel(roiCenter);

        dilateFilter->SetDilateValue(foreground);

        dilateFilter->Update();

        outputImageROI = dilateFilter->GetOutput();
    }
    cout << roiCenter << "Foreground: " << static_cast<int>(foreground) << endl;

    // typedef itk::BinaryMorphologicalClosingImageFilter<OutImageType, OutImageType, StructuringElementType> ClosingFilterType;
    // ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();

    // closingFilter->SetKernel(structuringElement3);
    // closingFilter->SetForegroundValue(1);
    // closingFilter->SetInput(filter->GetOutput());

    // closingFilter->Update();

    // //    filter->Update();


    std::cout << "Done running filter " << std::endl;

    // allocate outputImage first
    outputImage->CopyInformation( labelImage );
    outputImage->SetBufferedRegion( labelImage->GetBufferedRegion() );
    outputImage->Allocate();
    outputImage->FillBuffer(0);

    itk::ImageRegionIterator< OutImageType > filterOut(outputImageROI, outputImageROI->GetBufferedRegion());
    itk::ImageRegionIterator< OutImageType > out(outputImage, oRegion);

    for (filterOut.GoToBegin(), out.GoToBegin(); !filterOut.IsAtEnd(); ++filterOut, ++out)
    {
        if (filterOut.Get() == foreground)
            out.Set(1);
    }


    //std::cout << image << std::endl;
    //std::cout << outputImage << std::endl;


    std::cout << "writeImage_spacing = " << outputImage->GetSpacing() << std::endl ;
    std::cout << "writeImage_origin  = "  << outputImage->GetOrigin() << std::endl ;
    std::cout << "writeImage_LargestPossibleRegion = " << outputImage->GetLargestPossibleRegion() << std::endl ;


    ////////To write Output images /////////////////////
    ImageNameGenerator->SetOutputDirectory( argv[3] );
    FileNamesContainer OutputImageNames;
    OutputImageNames = ImageNameGenerator->GetOutputFileNames();

    OutImageWriterType::Pointer OutImageWriter = OutImageWriterType::New();
    OutImageWriter->SetInput(outputImage);
    OutImageWriter->SetImageIO(gdcmIO);
    OutImageWriter->SetFileNames(OutputImageNames);
    OutImageWriter->SetMetaDataDictionaryArray(InImageReader->GetMetaDataDictionaryArray());

    SingleOutImageWriterType::Pointer OutWriter = SingleOutImageWriterType::New();
    OutWriter->SetInput(outputImage);
    OutWriter->SetFileName("../output.mha");
    OutWriter->Update();

    try
    {
        OutImageWriter->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
        std::cerr << "Exception thrown while writing the series " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }


    std::cout << "Saved the tumor mask image" << std::endl;


    return EXIT_SUCCESS;
}
