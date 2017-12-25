#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConstantPadImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"

#include "itkAreaOpeningImageFilter.h"

#include "itkBinaryFillholeImageFilter.h"
#include "itkSliceBySliceImageFilter.h"

#include "itkOtsuMultipleThresholdsImageFilter.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkCylinderSpatialObject.h"

#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"

#include "itkImageDuplicator.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#include "itkGrowCutSegmentationImageFilter.h"

using namespace std;

// Definition of Data type
typedef float InputPixelType;
typedef unsigned char OutputPixelType;
typedef float InternalPixelType;


typedef itk::Image<InputPixelType, 3> InImageType;
typedef itk::Image<OutputPixelType, 3> OutImageType;
typedef itk::Image<OutputPixelType, 2> OutImage2DType;
typedef itk::Image<InternalPixelType, 3> WeightImageType;

// Definition of IO
typedef itk::ImageFileReader< InImageType > InImageReaderType;
typedef itk::ImageFileWriter< InImageType > InImageWriterType;
typedef itk::ImageFileWriter< WeightImageType > WeightImageWriterType;
typedef itk::ImageFileWriter< OutImageType> OutImageWriterType;

typedef itk::ImageDuplicator< OutImageType > DuplicatorType;

typedef itk::RegionOfInterestImageFilter< InImageType, InImageType > RoiFilterType;


typedef itk::BinaryBallStructuringElement< OutputPixelType, 3 > StructuringElementType;
typedef itk::BinaryBallStructuringElement< OutputPixelType, 2 > StructuringElement2DType;
typedef itk::BinaryDilateImageFilter<OutImageType, OutImageType, StructuringElementType> DilateFilterType;
typedef itk::BinaryErodeImageFilter<OutImageType, OutImageType, StructuringElementType> ErodeFilterType;

typedef itk::BinaryMorphologicalClosingImageFilter<OutImage2DType, OutImage2DType, StructuringElement2DType> ClosingFilter2DType;
typedef itk::BinaryMorphologicalOpeningImageFilter<OutImage2DType, OutImage2DType, StructuringElement2DType> OpeningFilter2DType;


typedef itk::SliceBySliceImageFilter< OutImageType, OutImageType> SliceBySliceFilterType;
typedef itk::BinaryFillholeImageFilter<OutImage2DType > FillholeFilter2DType;

typedef itk::AreaOpeningImageFilter <OutImageType, OutImageType> AreaOpeningImageFilterType;
typedef itk::ConnectedComponentImageFilter <OutImageType, OutImageType> ConnectedComponentImageFilterType;

typedef itk::CurvatureAnisotropicDiffusionImageFilter<InImageType, WeightImageType > DiffusionFilterType;

typedef itk::BinaryThresholdImageFilter <WeightImageType, OutImageType >  ThresholdFilterType;
typedef itk::ConnectedComponentImageFilter <OutImageType, OutImageType > ConnectedComponentImageFilterType;

typedef itk::GrowCutSegmentationImageFilter<WeightImageType, OutImageType> GrowCutFilterType;

typedef itk::SubtractImageFilter< OutImageType, OutImageType, OutImageType > SubFilterType;

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
        cout << "\rProgress: "
                  << static_cast<unsigned int>( 100.0 * m_Process->GetProgress() ) << "%";
		if(m_Process->GetProgress() == 1)
		{
			cout << endl;
		}
    }
    itk::ProcessObject::Pointer m_Process;
};


int main( int argc, char *argv[] )
{
    if ( argc < 8 )
    {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage: " << argv[0];
        cerr << " InputImageFile SeedPoint[x y z] NoduleSize[long short] OutputImageFile" << endl;
        return EXIT_FAILURE;
    }

    unsigned int x = atoi( argv[2] );
    unsigned int y = atoi( argv[3] );
    int sliceNumber = atoi( argv[4] );
    double longD = atof( argv[5] );
    double shortD = atof( argv[6] );

    if (longD <= 10) shortD = longD * 0.8;
    if (shortD < 4) shortD = 4;

    double objectSize = longD;
    double priorSegmentStrength = 1.0;
    double contrastNoiseRatio = 0.003;

    double shortR = shortD / 2;
    double longR = longD / 2;

    unsigned char foreground = 255;
    unsigned char background = 0;


    /////////To Read CT images ////////////////////
    InImageReaderType::Pointer InImageReader = InImageReaderType::New();
    InImageReader->SetFileName(argv[1]);
    InImageReader->Update();

    InImageType::Pointer oimage = InImageReader->GetOutput();
    InImageType::Pointer image = InImageReader->GetOutput();

    if (contrastNoiseRatio > 1.0)
    {
        contrastNoiseRatio /= 100.0;
    }

    if (priorSegmentStrength > 1.0)
    {
        priorSegmentStrength /= 100.0;
    }

    // set ROI
    InImageType::IndexType index = image->GetLargestPossibleRegion().GetIndex();
    InImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    InImageType::PointType origin  = image->GetOrigin();
    InImageType::SpacingType spacing = image->GetSpacing();
    InImageType::SizeType roiSize;

    roiSize[0] = static_cast<unsigned int>((longD + 20) / spacing[0]);
    roiSize[1] = static_cast<unsigned int>((longD + 20) / spacing[1]);
    roiSize[2] = static_cast<unsigned int>((longD + 20) / spacing[2]);

    if(size[0] < roiSize[0] || size[1] < roiSize[1] || size[2] < roiSize[2])
    {
      typedef itk::ConstantPadImageFilter <InImageType, InImageType> ConstantPadImageFilterType;

      cout << "Padding " << endl;
      cout << "Input Image Spacing = " << spacing << endl;
      cout << "Input Image Origin = " << origin << endl;
      cout << "Input Image Size = " << size << endl<< endl << endl;

      InImageType::PointType diff;
      diff[0] = roiSize[0]-size[0];
      diff[1] = roiSize[1]-size[1];
      diff[2] = roiSize[2]-size[2];
      diff[0] = diff[0]<0?0:diff[0];
      diff[1] = diff[1]<0?0:diff[1];
      diff[2] = diff[2]<0?0:diff[2];

      InImageType::SizeType lowerExtendRegion;
      lowerExtendRegion[0] = diff[0];
      lowerExtendRegion[1] = diff[1];
      lowerExtendRegion[2] = diff[2];

      InImageType::SizeType upperExtendRegion;
      upperExtendRegion[0] = diff[0];
      upperExtendRegion[1] = diff[1];
      upperExtendRegion[2] = diff[2];

      cout << x << "," << y << "," << sliceNumber << endl;

      x += diff[0];
      y += diff[1];
      sliceNumber += diff[2];

      cout << x << "," << y << "," << sliceNumber << endl;

      InImageType::PixelType constantPixel = 0;

      ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
      padFilter->SetInput(image);
      padFilter->SetPadLowerBound(lowerExtendRegion);
      padFilter->SetPadUpperBound(upperExtendRegion);
      padFilter->SetConstant(constantPixel);
      padFilter->Update();

      image = padFilter->GetOutput();
      image->SetOrigin(origin-diff);
    }

    index = image->GetLargestPossibleRegion().GetIndex();
    size = image->GetLargestPossibleRegion().GetSize();
    origin  = image->GetOrigin();
    spacing = image->GetSpacing();

    cout << "Input Image Spacing = " << spacing << endl;
    cout << "Input Image Origin = " << origin << endl;
    cout << "Input Image Size = " << size << endl<< endl << endl;

    OutImageType::IndexType roiStart;
    OutImageType::IndexType roiEnd;

    OutImageType::IndexType center;

    center[0] = x;
    center[1] = y;

    if (sliceNumber > 0)
        center[2] = size[2] - sliceNumber - 1;
    else
        center[2] = -sliceNumber;

    cout << center << endl;

    roiStart[0] = center[0] - roiSize[0] / 2 + 1;
    roiStart[1] = center[1] - roiSize[1] / 2 + 1;
    roiStart[2] = center[2] - roiSize[2] / 2 + 1;
    roiEnd[0] = roiStart[0] + roiSize[0];
    roiEnd[1] = roiStart[1] + roiSize[1];
    roiEnd[2] = roiStart[2] + roiSize[2];

    cout << roiSize << endl;
    cout << roiStart << endl;
    cout << roiEnd << endl;


    // Definition of Morphological operation sturcture elements
    StructuringElementType::RadiusType elementRadius;

    StructuringElement2DType structuringElement2D1;
    structuringElement2D1.SetRadius(1);
    structuringElement2D1.CreateStructuringElement();

    StructuringElementType structuringElement1;
    structuringElement1.SetRadius(1);
    structuringElement1.CreateStructuringElement();

    elementRadius.Fill(static_cast<unsigned int>(round(3/spacing[0])));
    elementRadius[2] = static_cast<unsigned int>(round(3/spacing[2]));

    StructuringElementType structuringElement3;
    structuringElement3.SetRadius(elementRadius);
    structuringElement3.CreateStructuringElement();

    elementRadius.Fill(static_cast<unsigned int>(round(5/spacing[0])));
    elementRadius[2] = static_cast<unsigned int>(round(5/spacing[2]));

    StructuringElementType structuringElement5;
    structuringElement5.SetRadius(elementRadius);
    structuringElement5.CreateStructuringElement();

    elementRadius.Fill(static_cast<unsigned int>(round(10/spacing[0])));
    elementRadius[2] = static_cast<unsigned int>(round(10/spacing[2]));

    StructuringElementType structuringElement10;
    structuringElement10.SetRadius(elementRadius);
    structuringElement10.CreateStructuringElement();


    unsigned int ndims = image->GetImageDimension();

    cout << " objectSize (radius) " << objectSize << endl;
    OutImageType::PixelType radius = static_cast< OutImageType::PixelType> (objectSize);

    /////////////////////////////////////////////////////////////////

    InImageType::IndexType istart;
    InImageType::SizeType isize;

    OutImageType::IndexType ostart;
    OutImageType::SizeType osize;

    WeightImageType::IndexType wstart;
    WeightImageType::SizeType wsize;

    OutImageType::IndexType roiCenter;

    for (unsigned n = 0; n < ndims; n++)
    {
        istart[n] = roiStart[n] + origin[n];
        isize[n] = roiEnd[n] - roiStart[n];

        ostart[n] = istart[n];
        osize[n] = isize[n];

        roiCenter[n] = static_cast<OutputPixelType>(round(isize[n] / 2));

        wstart[n] = istart[n];
        wsize[n] = isize[n];
    }

    cout << " roiStart " << roiStart << " roiEnd " << roiEnd << " " << endl;
    cout << " istart " << istart << " isize " << isize << " " << endl;
    cout << " roiCenter " << roiCenter << endl;



    InImageType::RegionType iRegion;
    iRegion.SetSize( isize );
    iRegion.SetIndex( istart );

    InImageType::RegionType oRegion;
    oRegion.SetSize( osize );
    oRegion.SetIndex( ostart );

    InImageType::RegionType wRegion;
    wRegion.SetSize( wsize );
    wRegion.SetIndex( wstart );


    RoiFilterType::Pointer fInput = RoiFilterType::New();
    //cout << image ;
    //cout << iRegion;
    fInput->SetRegionOfInterest( iRegion );
    fInput->SetInput( image );
    fInput->Update();


    InImageType::Pointer inImage = InImageType::New();
    inImage = fInput->GetOutput();
    origin[0] = istart[0];
    origin[1] = istart[1];
    origin[2] = istart[2];
    inImage->SetOrigin(origin);

    //cout << inImage ;

    OutImageType::Pointer labImage = OutImageType::New();
    WeightImageType::Pointer wtImage = WeightImageType::New();

    labImage->CopyInformation(inImage);
    labImage->SetBufferedRegion( inImage->GetBufferedRegion() );
    labImage->Allocate();
    labImage->FillBuffer( 0 );

    wtImage->CopyInformation(inImage);
    wtImage->SetBufferedRegion( inImage->GetBufferedRegion() );
    wtImage->Allocate();
    wtImage->FillBuffer( 0 );

    /////////////////////////////////////////////////////////////////

    DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();

    // Diffusion filter parameters
    diffusionFilter->SetInput(inImage);

    diffusionFilter->SetNumberOfIterations( 5 );
    diffusionFilter->SetTimeStep( (spacing[0] < spacing[2]) ? (spacing[0] / 16) : (spacing[2] / 16));
    diffusionFilter->SetConductanceParameter( 2 );
    diffusionFilter->UseImageSpacingOn();
    diffusionFilter->Update();

    // {
    //     InImageWriterType::Pointer writer = InImageWriterType::New();
    //     writer->SetUseCompression(true);
    //     writer->SetInput(inImage);
    //     writer->SetFileName(string(argv[7]) + "in.nrrd");
    //     writer->Update();
    // }
    //
    // {
    //     WeightImageWriterType::Pointer writer = WeightImageWriterType::New();
    //     writer->SetUseCompression(true);
    //     writer->SetInput(diffusionFilter->GetOutput());
    //     writer->SetFileName(string(argv[7]) + "filt.nrrd");
    //     writer->Update();
    // }

    /////////////////////////////////////////////////////////////////
    {
        InImageType::IndexType start;
        InImageType::SizeType size;
        InImageType::RegionType centerRegion;

        size.Fill(4);
        start[0] = roiCenter[0] - 1;
        start[1] = roiCenter[1] - 1;
        start[2] = roiCenter[2] - 1;
        centerRegion.SetSize(size);
        centerRegion.SetIndex(start);

        cout << centerRegion << endl;

        //itk::ImageRegionIterator< InImageType > centerPixels(image, centerRegion);
        //itk::ImageRegionIterator< InImageType > centerPixels(inImage, centerRegion);
        itk::ImageRegionIterator< WeightImageType > centerPixels(diffusionFilter->GetOutput(), centerRegion);
        double Tc = 0.5;
        double pTc = 0;
        double varTc = 0;
        double stdTc = 0.4;
        double numPixels = 4 * 4 * 4;
        int iter = 0;

        while(numPixels > 4 && iter++ < 50 && stdTc > 0.01)
        {
            double maxPix = 0;
            double sumPix = 0;
            double sumPix2 = 0;
            double Tc2 = 0;

            numPixels = 4 * 4 * 4;

            centerPixels.GoToBegin();
            while (!centerPixels.IsAtEnd())
            {
                double pixel = centerPixels.Get();
                if(pixel > maxPix)
                    maxPix = pixel;
                if (pixel > Tc-2*stdTc && pixel < Tc+2*stdTc)
                {
                    sumPix = sumPix + pixel;
                    sumPix2 = sumPix2 + pixel*pixel;
                    cout << pixel << " ";
                }
                else
                {
                    numPixels--;
                    cout << "(" << pixel << ") ";
                }

                ++centerPixels;
            }
            pTc = Tc;
            Tc = sumPix / numPixels;
            Tc2 = sumPix2 / numPixels;
            varTc = Tc2 - Tc*Tc;
            stdTc = sqrt(varTc);

            if(stdTc > 0.2 || (pTc == Tc && stdTc > 0.1))
            {
                Tc = maxPix;
                stdTc = stdTc/2;
                cout << endl << "Reset initial: Tc="<< Tc <<" stdTc=" << stdTc << endl;
                continue;
            }

            cout << endl;
            cout << "--------" << endl;
            cout << numPixels << "     => " << Tc << "+-" << stdTc << endl;
        }

        if(stdTc > 0.2) stdTc = 0.2;
        cout << "Threshold: " << Tc - 3*stdTc << endl;

        ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
        ThresholdFilterType::Pointer thresholdFilter_bg = ThresholdFilterType::New();
        //thresholdFilter->SetInput(inImage);
        //thresholdFilter_bg->SetInput(inImage);
        thresholdFilter->SetInput(diffusionFilter->GetOutput());
        thresholdFilter_bg->SetInput(diffusionFilter->GetOutput());

        thresholdFilter->SetLowerThreshold(Tc-3*stdTc);
        thresholdFilter->SetUpperThreshold(Tc+3*stdTc);

        //thresholdFilter_bg->SetLowerThreshold(Tc-stdTc*2);
        //thresholdFilter_bg->SetUpperThreshold(Tc+200);

        thresholdFilter_bg->SetLowerThreshold(Tc-3*stdTc);
        thresholdFilter_bg->SetUpperThreshold(Tc+3*stdTc);

        thresholdFilter->SetInsideValue(1);
        thresholdFilter->SetOutsideValue(0);
        thresholdFilter->Update();

        thresholdFilter_bg->SetInsideValue(0);
        thresholdFilter_bg->SetOutsideValue(1);
        thresholdFilter_bg->Update();

        OutImageType::Pointer thresholdImage = thresholdFilter->GetOutput();
        {
            SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

            FillholeFilter2DType::Pointer fillholeFilter2D = FillholeFilter2DType::New();
            fillholeFilter2D->SetFullyConnected(true);
            fillholeFilter2D->SetForegroundValue(1);

            sliceBySliceFilter->SetInput(thresholdImage);
            sliceBySliceFilter->SetFilter(fillholeFilter2D);

            sliceBySliceFilter->Update();

            thresholdImage = sliceBySliceFilter->GetOutput();
        }

        {
            SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

            OpeningFilter2DType::Pointer openingFilter2D = OpeningFilter2DType::New();
            openingFilter2D->SetKernel(structuringElement2D1);
            openingFilter2D->SetForegroundValue(1);

            sliceBySliceFilter->SetInput(thresholdImage);
            sliceBySliceFilter->SetFilter(openingFilter2D);

            sliceBySliceFilter->Update();

            thresholdImage = sliceBySliceFilter->GetOutput();
        }


        {
    		AreaOpeningImageFilterType::Pointer openingFilter = AreaOpeningImageFilterType::New();
	    	openingFilter->SetInput(thresholdImage);
		    openingFilter->SetLambda(3);
            //openingFilter->FullyConnectedOn();
    		openingFilter->Update();
            cout <<  "Opening" << endl;
            thresholdImage = openingFilter->GetOutput();
	    	//thresholdImage->SetPixel(roiCenter,1);
        }

        cout <<  "Thresholding" << endl;


        /////////////////////////////////////////////////////////////////

		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
		connected->SetInput(thresholdImage);
		connected->Update();


        /////////////////////////////////////////////////////////////////
        double target = 0;
        {
            itk::ImageRegionIterator< OutImageType > centerPixels(connected->GetOutput(), centerRegion);
            centerPixels.GoToBegin();
            while (!centerPixels.IsAtEnd())
            {
                double label = centerPixels.Get();
                if(label > 0)
                {
                    target = label;
                }

                ++centerPixels;
            }
        }
        cout << target << endl;



        OutImageType::Pointer targetImage = OutImageType::New();
        int numSeeds = 0;

        targetImage->CopyInformation(inImage);
        targetImage->SetBufferedRegion( inImage->GetBufferedRegion() );
        targetImage->Allocate();
        targetImage->FillBuffer( 0 );
        {
            itk::ImageRegionIterator< OutImageType > label(connected->GetOutput(), connected->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator< OutImageType > targetLabel(targetImage, targetImage->GetBufferedRegion());
            label.GoToBegin();
            targetLabel.GoToBegin();
            while (!label.IsAtEnd())
            {
                if(label.Get() == target)
                {
                    targetLabel.Set(1);
                    ++numSeeds;
                }
                else targetLabel.Set(0);
                ++label;
                ++targetLabel;
            }
        }

        cout << numSeeds << endl;

        {
            SubFilterType::Pointer subFilter = SubFilterType::New();
            subFilter->SetInput1(targetImage);
            subFilter->SetInput2(thresholdFilter_bg->GetOutput());
            subFilter->Update();

            targetImage = subFilter->GetOutput();

            itk::ImageRegionIterator< OutImageType > targetLabel(targetImage, targetImage->GetBufferedRegion());
            targetLabel.GoToBegin();
            while (!targetLabel.IsAtEnd())
            {
                if(targetLabel.Get() != 1) targetLabel.Set(0);
                ++targetLabel;
            }

     		AreaOpeningImageFilterType::Pointer openingFilter = AreaOpeningImageFilterType::New();
	    	openingFilter->SetInput(targetImage);
		    openingFilter->SetLambda(5);
    		openingFilter->Update();
            targetImage = openingFilter->GetOutput();
        }


        ErodeFilterType::Pointer erodeFilterB = ErodeFilterType::New();
		erodeFilterB->SetInput(thresholdFilter_bg->GetOutput());
		erodeFilterB->SetKernel(structuringElement1);
        erodeFilterB->SetErodeValue(1);
		erodeFilterB->Update();


		DilateFilterType::Pointer dilateFilter5 = DilateFilterType::New();
		dilateFilter5->SetInput(targetImage);
		dilateFilter5->SetKernel(structuringElement5);
        dilateFilter5->SetDilateValue(1);
		dilateFilter5->Update();

		DilateFilterType::Pointer dilateFilter10 = DilateFilterType::New();
		dilateFilter10->SetInput(targetImage);
		dilateFilter10->SetKernel(structuringElement10);
        dilateFilter10->SetDilateValue(1);
		dilateFilter10->Update();
        {
            ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    		connected->SetInput(dilateFilter10->GetOutput());
	    	connected->Update();

            OutImageType::Pointer image = connected->GetOutput();

            itk::ImageRegionIterator< OutImageType > centerLabels(image, centerRegion);

            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(image);
            duplicator->Update();
            OutImageType::Pointer tempImage = duplicator->GetOutput();


            double target = connected->GetOutput()->GetPixel(roiCenter);

            OutImageType::SizeType roiSize;
            roiSize[0] = static_cast<unsigned int>((longD + 6) / spacing[0]);
            roiSize[1] = static_cast<unsigned int>((longD + 6) / spacing[1]);
            roiSize[2] = static_cast<unsigned int>((longD + 6) / spacing[2]);

            OutImageType::IndexType roiStart;
            OutImageType::IndexType roiEnd;

            roiStart[0] = roiCenter[0] - roiSize[0] / 2;
            roiStart[1] = roiCenter[1] - roiSize[1] / 2;
            roiStart[2] = roiCenter[2] - roiSize[2] / 2;
            roiEnd[0] = roiStart[0] + roiSize[0];
            roiEnd[1] = roiStart[1] + roiSize[1];
            roiEnd[2] = roiStart[2] + roiSize[2];

            cout << roiSize << endl;
            cout << roiStart << endl;
            cout << roiEnd << endl;


            itk::ImageRegionIterator< OutImageType > tempLabel(tempImage, tempImage->GetBufferedRegion());
            itk::ImageRegionIteratorWithIndex< OutImageType > label(dilateFilter10->GetOutput(), image->GetBufferedRegion());
            label.GoToBegin();
            tempLabel.GoToBegin();
            while (!label.IsAtEnd())
            {
                OutImageType::IndexType idx = label.GetIndex();

                if (idx[0] >= roiStart[0] && idx[1] >= roiStart[1] && idx[2] >= roiStart[2] &&
                    idx[0] < roiEnd[0] && idx[1] < roiEnd[1] && idx[2] < roiEnd[2])
                {
                    if(tempLabel.Get() == target) label.Set(1);
                }
                else label.Set(0);
                ++label;
                ++tempLabel;
            }
        }


        // generate seed areas
        {
            //itk::ImageRegionIterator< OutImageType > bg(erodeFilterB->GetOutput(), erodeFilterB->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator< OutImageType > bg(thresholdFilter_bg->GetOutput(), erodeFilterB->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator< OutImageType > gray(dilateFilter5->GetOutput(), dilateFilter5->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator< OutImageType > truebg(dilateFilter10->GetOutput(), dilateFilter10->GetOutput()->GetBufferedRegion());

            itk::ImageRegionIterator< OutImageType > targetLabel(targetImage, targetImage->GetBufferedRegion());
            itk::ImageRegionIterator< OutImageType > newLabel(labImage, labImage->GetBufferedRegion());
            itk::ImageRegionIterator< WeightImageType > weight(wtImage, wtImage->GetBufferedRegion());

            targetLabel.GoToBegin();
            newLabel.GoToBegin();
            weight.GoToBegin();

            bg.GoToBegin();
            gray.GoToBegin();
            truebg.GoToBegin();
            while (!targetLabel.IsAtEnd())
            {
                weight.Set(0);
                if(gray.Get() == 0)
                {
                    weight.Set(contrastNoiseRatio);
                    newLabel.Set(2);
                }

                if(targetLabel.Get() == 1) // taget
                {
                    weight.Set(priorSegmentStrength);
                    newLabel.Set(1);
                }


                if(bg.Get() == 1) // background
                {
                    weight.Set(contrastNoiseRatio);
                    //weight.Set(priorSegmentStrength/2);
                    newLabel.Set(2);
                }

                if(truebg.Get() == 0) // definetely background, opposite value (0 is bg)
                {
                    weight.Set(priorSegmentStrength);
                    newLabel.Set(2);
                }

                ++targetLabel;
                ++newLabel;
                ++weight;
                ++bg;
                ++gray;
                ++truebg;
            }
        }

        // {
        //     OutImageWriterType::Pointer writer = OutImageWriterType::New();
        //     writer->SetUseCompression(true);
        //     writer->SetInput(labImage);
        //     writer->SetFileName(string(argv[7]) + "Tc-label.nrrd");
        //     writer->Update();
        // }

    }

    cout << "label image generation completed" << endl;

    /////////////////////////////////////////////////////////////////

    // GrowCut Segmentation
    GrowCutFilterType::Pointer filter = GrowCutFilterType::New();

    // Segmentation filter parameters
    filter->SetInput( diffusionFilter->GetOutput() );
    filter->SetLabelImage( labImage );

    filter->SetStrengthImage( wtImage );

    filter->SetSeedStrength( contrastNoiseRatio );
    filter->SetObjectRadius((unsigned int)objectSize);

    //
    ShowProgressObject progressWatch(filter);
    itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
    command = itk::SimpleMemberCommand<ShowProgressObject>::New();
    command->SetCallbackFunction(&progressWatch,
                                 &ShowProgressObject::ShowProgress);
    filter->AddObserver(itk::ProgressEvent(), command);

    filter->Update();


    OutImageType::Pointer outputImageROI = OutImageType::New();

    outputImageROI = filter->GetOutput();

    cout << "Done running segmentation filter " << endl;

    /////////////////////////////////////////////////////////////////

    // Post processing
    {
        itk::ImageRegionIterator< OutImageType > label(outputImageROI, outputImageROI->GetBufferedRegion());
        label.GoToBegin();
        while (!label.IsAtEnd())
        {
            if(label.Get() != 1) label.Set(0);
            ++label;
        }

        SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

        OpeningFilter2DType::Pointer openingFilter2D = OpeningFilter2DType::New();
        openingFilter2D->SetKernel(structuringElement2D1);
        openingFilter2D->SetForegroundValue(1);

        sliceBySliceFilter->SetInput(outputImageROI);
        sliceBySliceFilter->SetFilter(openingFilter2D);

        sliceBySliceFilter->Update();

        outputImageROI = sliceBySliceFilter->GetOutput();

        AreaOpeningImageFilterType::Pointer openingFilter = AreaOpeningImageFilterType::New();
        openingFilter->SetInput(outputImageROI);
        openingFilter->SetLambda(3);
        openingFilter->FullyConnectedOn();
        openingFilter->Update();
        cout <<  "Opening" << endl;
        outputImageROI = openingFilter->GetOutput();
    }

    // {
    //     OutImageWriterType::Pointer writer = OutImageWriterType::New();
    //     writer->SetUseCompression(true);
    //     writer->SetInput(outputImageROI);
    //     writer->SetFileName(string(argv[7]) + "mask-label.nrrd");
    //     writer->Update();
    // }


    /////////////////////////////////////////////////////////////////
    OutImageType::Pointer outputImage = OutImageType::New();
    // allocate outputImage first
    outputImage->CopyInformation( oimage );
    outputImage->SetSpacing(spacing);
    outputImage->SetBufferedRegion( oimage->GetBufferedRegion() );
    outputImage->Allocate();
    outputImage->FillBuffer(0);

    //cout << outputImageROI;
    oRegion = outputImage->GetLargestPossibleRegion();
    oRegion.SetIndex(0, oRegion.GetIndex(0) - origin[0]);
    oRegion.SetIndex(1, oRegion.GetIndex(1) - origin[1]);
    oRegion.SetIndex(2, oRegion.GetIndex(2) - origin[2]);

    itk::ImageRegionIterator< OutImageType > filterOut(outputImageROI, oRegion);
    itk::ImageRegionIterator< OutImageType > out(outputImage, outputImage->GetLargestPossibleRegion());

    for (filterOut.GoToBegin(), out.GoToBegin(); !filterOut.IsAtEnd(); ++filterOut, ++out)
    {
        if (filterOut.Get() == 1)
            out.Set(1);
        //out.Set(filterOut.Get());
    }

    cout << "writeImage_spacing = " << outputImage->GetSpacing() << endl ;
    cout << "writeImage_origin  = "  << outputImage->GetOrigin() << endl ;
    cout << "writeImage_LargestPossibleRegion = " << outputImage->GetLargestPossibleRegion() << endl ;


    ////////To write Output images /////////////////////
    OutImageWriterType::Pointer OutImageWriter = OutImageWriterType::New();
    OutImageWriter->SetUseCompression(true);
    OutImageWriter->SetInput(outputImage);
    OutImageWriter->SetFileName( argv[7] );

    try
    {
        OutImageWriter->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
        cerr << "Exception thrown while writing the series " << endl;
        cerr << excp << endl;
        return EXIT_FAILURE;
    }


    cout << "Saved the tumor mask image" << endl;


    return EXIT_SUCCESS;
}
