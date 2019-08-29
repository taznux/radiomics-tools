#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageSpatialObject.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryClosingByReconstructionImageFilter.h"
#include "itkBinaryOpeningByReconstructionImageFilter.h"

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

#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"

#include "itkImageDuplicator.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkIsotropicResamplerImageFilter.h"

#include "itkSatoVesselnessSigmoidFeatureGenerator.h"

#include "itkGrowCutSegmentationImageFilter.h"

using namespace std;

// Definition of Data type
using InputPixelType = signed short;
using OutputPixelType = unsigned char;
using InternalPixelType = float;


using InImageType = itk::Image<InputPixelType, 3>;
using OutImageType = itk::Image<OutputPixelType, 3>;
using OutImage2DType = itk::Image<OutputPixelType, 2>;
using WeightImageType = itk::Image<InternalPixelType, 3>;


// Definition of IO
using InImageReaderType = itk::ImageFileReader<InImageType>;
using InImageWriterType = itk::ImageFileWriter<InImageType>;
using WeightImageWriterType = itk::ImageFileWriter<WeightImageType>;
using OutImageWriterType = itk::ImageFileWriter<OutImageType>;

using DuplicatorType = itk::ImageDuplicator<OutImageType>;

using RoiFilterType = itk::RegionOfInterestImageFilter<InImageType, InImageType>;


using StructuringElementType = itk::BinaryBallStructuringElement<OutputPixelType, 3>;
using StructuringElement2DType = itk::BinaryBallStructuringElement<OutputPixelType, 2>;
using DilateFilterType = itk::BinaryDilateImageFilter<OutImageType, OutImageType, StructuringElementType>;
using ErodeFilterType = itk::BinaryErodeImageFilter<OutImageType, OutImageType, StructuringElementType>;

using ClosingFilter2DType = itk::BinaryClosingByReconstructionImageFilter<OutImage2DType, StructuringElement2DType>;
using OpeningFilter2DType = itk::BinaryOpeningByReconstructionImageFilter<OutImage2DType, StructuringElement2DType>;


using SliceBySliceFilterType = itk::SliceBySliceImageFilter<OutImageType, OutImageType>;
using FillholeFilter2DType = itk::BinaryFillholeImageFilter<OutImage2DType>;

using AreaOpeningImageFilterType = itk::AreaOpeningImageFilter<OutImageType, OutImageType>;
using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<OutImageType, OutImageType>;

using DiffusionFilterType = itk::CurvatureAnisotropicDiffusionImageFilter<WeightImageType, WeightImageType>;
using IsotropicResamplerType = itk::IsotropicResamplerImageFilter<InImageType, WeightImageType>;
using VesselnessGeneratorType = itk::SatoVesselnessSigmoidFeatureGenerator<3>;
using SpatialObjectType = VesselnessGeneratorType::SpatialObjectType;
using InputImageSpatialObjectType = VesselnessGeneratorType::InputImageSpatialObjectType;
using OutputImageSpatialObjectType = itk::ImageSpatialObject<VesselnessGeneratorType::Dimension, InternalPixelType>;


using ThresholdFilterType = itk::BinaryThresholdImageFilter<WeightImageType, OutImageType>;
using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<OutImageType, OutImageType>;

using GrowCutFilterType = itk::GrowCutSegmentationImageFilter<WeightImageType, OutImageType>;

using SubFilterType = itk::SubtractImageFilter<OutImageType, OutImageType, OutImageType>;

// To print the progress
class ShowProgressObject
{
public:
    ShowProgressObject(itk::ProcessObject * o)
    {
        m_Process = o;
    }

    void
    ShowProgress()
    {
        cout << "\rProgress: "
             << static_cast<unsigned int>(100.0 * m_Process->GetProgress()) << "%";
        if (m_Process->GetProgress() == 1) {
            cout << endl;
        }
    }

    itk::ProcessObject::Pointer m_Process;
};


int
main(int argc, char * argv[])
{
    if (argc < 8) {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage: " << argv[0];
        cerr << " InputImageFile SeedPoint[x y z] NoduleSize[long short] OutputImageFile" << endl;
        return EXIT_FAILURE;
    }

    const unsigned int x  = atoi(argv[2]);
    const unsigned int y  = atoi(argv[3]);
    const int sliceNumber = atoi(argv[4]);
    double longD  = atof(argv[5]);
    double shortD = atof(argv[6]);

    double objectSize = longD;
    double priorSegmentStrength = 0.003;
    double contrastNoiseRatio   = 0.8;

    double shortR = shortD / 2;
    double longR  = longD / 2;

    unsigned char foreground = 255;
    unsigned char background = 0;


    /////////To Read CT images ////////////////////
    InImageReaderType::Pointer InImageReader = InImageReaderType::New();
    InImageReader->SetFileName(argv[1]);
    InImageReader->Update();


    InImageType::Pointer image = InImageReader->GetOutput();


    if (contrastNoiseRatio > 1.0) {
        contrastNoiseRatio /= 100.0;
    }

    if (priorSegmentStrength > 1.0) {
        priorSegmentStrength /= 100.0;
    }

    // set ROI
    InImageType::IndexType index     = image->GetLargestPossibleRegion().GetIndex();
    InImageType::SizeType size       = image->GetLargestPossibleRegion().GetSize();
    InImageType::SpacingType spacing = image->GetSpacing();
    InImageType::SizeType roiSize;
    InImageType::PointType center_phy;

    roiSize[0] = static_cast<unsigned int>((longD + 10) / spacing[0]);
    roiSize[1] = static_cast<unsigned int>((longD + 10) / spacing[1]);
    roiSize[2] = static_cast<unsigned int>((longD + 10) / spacing[2]);

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
    image->TransformIndexToPhysicalPoint(center, center_phy);
    cout << center_phy << endl;


    roiStart[0] = center[0] - roiSize[0] / 2;
    roiStart[1] = center[1] - roiSize[1] / 2;
    roiStart[2] = center[2] - roiSize[2] / 2;
    roiEnd[0]   = roiStart[0] + roiSize[0];
    roiEnd[1]   = roiStart[1] + roiSize[1];
    roiEnd[2]   = roiStart[2] + roiSize[2];

    cout << roiSize << endl;
    cout << roiStart << endl;
    cout << roiEnd << endl;


    unsigned int ndims = image->GetImageDimension();

    cout << " objectSize (radius) " << objectSize << endl;
    OutImageType::PixelType radius = static_cast<OutImageType::PixelType>(objectSize);

    /////////////////////////////////////////////////////////////////

    InImageType::IndexType istart;
    InImageType::SizeType isize;

    for (unsigned n = 0; n < ndims; n++) {
        roiStart[n] = roiStart[n] < 0 ? 0 : roiStart[n];
        istart[n]   = roiStart[n];
        isize[n]    = roiEnd[n] - roiStart[n];
        isize[n]    = (static_cast<unsigned int>(isize[n] + istart[n]) <= size[n]) ?
          isize[n] : (size[n] - istart[n] - 1);
    }

    cout << " roiStart " << roiStart << " roiEnd " << roiEnd << " " << endl;
    cout << " istart " << istart << " isize " << isize << " " << endl;


    InImageType::RegionType iRegion;
    iRegion.SetIndex(istart);
    iRegion.SetSize(isize);

    RoiFilterType::Pointer fInput = RoiFilterType::New();
    fInput->SetRegionOfInterest(iRegion);
    fInput->SetInput(image);
    fInput->Update();

    InImageType::Pointer inImage = InImageType::New();
    inImage = fInput->GetOutput();


    {
        InImageWriterType::Pointer writer = InImageWriterType::New();
        writer->SetUseCompression(true);
        writer->SetInput(inImage);
        writer->SetFileName(string(argv[7]) + "in.nrrd");
        writer->Update();
    }


    /////////////////////////////////////////////////////////////////

    WeightImageType::SpacingType outSpacing = image->GetSpacing();
    WeightImageType::Pointer inImageIso;
    if (outSpacing[2] != outSpacing[0]) {
		if(outSpacing[2] > outSpacing[0])
			outSpacing[2] = outSpacing[0];
		else {
			outSpacing[0] = outSpacing[2];
			outSpacing[1] = outSpacing[2];
		}

        cout << "Out spacing " << outSpacing << endl;

        IsotropicResamplerType::Pointer isotropicResampler = IsotropicResamplerType::New();

        isotropicResampler->SetInput(inImage);
        isotropicResampler->SetOutputSpacing(outSpacing);
        isotropicResampler->Update();

        inImageIso = isotropicResampler->GetOutput();
    } else {
        using CastImageFilterType = itk::CastImageFilter<InImageType, WeightImageType>;
        CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
        castImageFilter->SetInput(inImage);
        castImageFilter->Update();

        inImageIso = castImageFilter->GetOutput();
    }
    iRegion = inImageIso->GetLargestPossibleRegion();


    {
        WeightImageWriterType::Pointer writer = WeightImageWriterType::New();

        writer->SetUseCompression(true);
        writer->SetInput(inImageIso);
        writer->SetFileName(string(argv[7]) + "in-iso.nrrd");
        writer->Update();
    }

    // Definition of Morphological operation sturcture elements
    StructuringElementType::RadiusType elementRadius;

    StructuringElement2DType structuringElement2D1;
    structuringElement2D1.SetRadius(1);
    structuringElement2D1.CreateStructuringElement();

    StructuringElementType structuringElement1;
    structuringElement1.SetRadius(1);
    structuringElement1.CreateStructuringElement();

    /*elementRadius.Fill(static_cast<unsigned int>(round(3 / outSpacing[0])));

    StructuringElementType structuringElement3;
    structuringElement3.SetRadius(elementRadius);
    structuringElement3.CreateStructuringElement();


	elementRadius.Fill(static_cast<unsigned int>(round(5 / outSpacing[0])));

	StructuringElementType structuringElement5;
	structuringElement5.SetRadius(elementRadius);
	structuringElement5.CreateStructuringElement();*/

    elementRadius.Fill(static_cast<unsigned int>(round(10 / outSpacing[0])));

    StructuringElementType structuringElement10;
    structuringElement10.SetRadius(elementRadius);
    structuringElement10.CreateStructuringElement();

    elementRadius.Fill(static_cast<unsigned int>(round(15 / outSpacing[0])));

    StructuringElementType structuringElement15;
    structuringElement15.SetRadius(elementRadius);
    structuringElement15.CreateStructuringElement();


    OutImageType::IndexType roiCenter;
    inImageIso->TransformPhysicalPointToIndex(center_phy, roiCenter);
    for (unsigned n = 0; n < ndims; n++) {
        roiStart[n] = iRegion.GetIndex(n);
        roiEnd[n]   = roiStart[n] + iRegion.GetSize(n) - 1;
    }

    cout << " roiStart " << roiStart << " roiEnd " << roiEnd << " " << endl;
    cout << " roiCenter " << roiCenter << endl;


    OutImageType::Pointer labImage   = OutImageType::New();
    WeightImageType::Pointer wtImage = WeightImageType::New();

    labImage->CopyInformation(inImageIso);
    labImage->SetBufferedRegion(inImageIso->GetBufferedRegion());
    labImage->Allocate();
    labImage->FillBuffer(0);

    wtImage->CopyInformation(inImageIso);
    wtImage->SetBufferedRegion(inImageIso->GetBufferedRegion());
    wtImage->Allocate();
    wtImage->FillBuffer(0);


	DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
	// Diffusion filter parameters
	diffusionFilter->SetInput(inImageIso);
	diffusionFilter->SetNumberOfIterations(5);
	diffusionFilter->SetTimeStep(outSpacing[0] / 16);
	diffusionFilter->SetConductanceParameter(2);
	diffusionFilter->UseImageSpacingOn();
	diffusionFilter->Update();

    {
		WeightImageWriterType::Pointer writer = WeightImageWriterType::New();
		writer->SetUseCompression(true);
		writer->SetInput(diffusionFilter->GetOutput());
		writer->SetFileName(string(argv[7]) + "filt.nrrd");
		writer->Update();
    }
	inImageIso = diffusionFilter->GetOutput();


    DilateFilterType::Pointer dilateFilter10 = DilateFilterType::New();
    DilateFilterType::Pointer dilateFilter15 = DilateFilterType::New();

    /////////////////////////////////////////////////////////////////
    {
        InImageType::IndexType start;
        InImageType::SizeType size;
        InImageType::RegionType centerRegion;

        size.Fill(round(max(shortD/ outSpacing[0], 5.0)));
        start[0] = roiCenter[0] - round(max(shortD / 2 / outSpacing[0],2.0));
        start[1] = roiCenter[1] - round(max(shortD / 2 / outSpacing[1],2.0));
        start[2] = roiCenter[2] - round(max(shortD / 2 / outSpacing[2],2.0));
        centerRegion.SetSize(size);
        centerRegion.SetIndex(start);

        cout << centerRegion << endl;

        // itk::ImageRegionIterator< InImageType > centerPixels(image, centerRegion);
        // itk::ImageRegionIterator< InImageType > centerPixels(inImage, centerRegion);
        // itk::ImageRegionIterator< WeightImageType > centerPixels(diffusionFilter->GetOutput(), centerRegion);
        itk::ImageRegionIterator<WeightImageType> centerPixels(inImageIso, centerRegion);
		const double th_low = -800;
		const double th_high = 200;
        double Tc        = -300;
        double stdTc     = 250;
		double pTc = 0;
		double pStdTc = 0;
		double varTc = 0;
        double numPixels = 125;
        int iter         = 0;
		int numMaxPix = 0;

        while (	numPixels > 30 &&
				numMaxPix < 100 &&
				iter++ < 20 && 
				stdTc > 200 && 
				(pTc != Tc || pStdTc != stdTc)) {
            double sumPix  = 0;
            double sumPix2 = 0;
            double Tc2     = 0;

            numPixels = 0;
			numMaxPix = 0;

            centerPixels.GoToBegin();
            while (!centerPixels.IsAtEnd()) {
                double pixel = centerPixels.Get();
                if (pixel > th_high){
					numMaxPix++;
				}
                if ((pixel > th_low  && pixel < th_high) && (pixel > Tc - 1 * stdTc && pixel < Tc + 2 * stdTc)) {
                    sumPix  = sumPix + pixel;
                    sumPix2 = sumPix2 + pixel * pixel;
					numPixels++;
                    //cout << pixel << " ";
                } else {
                    //cout << "(" << pixel << ") ";
                }

                ++centerPixels;
            }
            pTc    = Tc;
            pStdTc = stdTc;
            if(numPixels>10){
              Tc     = sumPix / numPixels;
              Tc2    = sumPix2 / numPixels;
              varTc  = Tc2 - Tc * Tc;
              stdTc  = sqrt(varTc);
            }
            if(stdTc == 0) stdTc = 50;

            if (Tc < th_low || Tc > th_high || stdTc > 400  || stdTc < 50) {
                if (Tc > th_high) Tc = th_high-250;
				if (Tc < th_low) Tc = th_low+250;
                stdTc = 100;
                cout << endl << "Reset initial: Tc=" << Tc << " stdTc=" << stdTc << endl;
                continue;
            }

            cout << endl;
            cout << "--------" << endl;
            cout << numPixels << "     => " << Tc << "+-" << stdTc << endl;
        }
        cout << "Threshold: " << Tc - 2 * stdTc << endl;

        ThresholdFilterType::Pointer thresholdFilter    = ThresholdFilterType::New();
        ThresholdFilterType::Pointer thresholdFilter_bg = ThresholdFilterType::New();
        // thresholdFilter->SetInput(inImage);
        // thresholdFilter_bg->SetInput(inImage);
        thresholdFilter->SetInput(inImageIso);
        thresholdFilter_bg->SetInput(inImageIso);

        thresholdFilter->SetLowerThreshold(Tc - 2 * stdTc);
        thresholdFilter->SetUpperThreshold(5000);

        thresholdFilter_bg->SetLowerThreshold(th_low);
        thresholdFilter_bg->SetUpperThreshold(5000);

        // thresholdFilter_bg->SetLowerThreshold(Tc - 200);
        // thresholdFilter_bg->SetUpperThreshold(Tc + 400);

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

        /*{
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
            // openingFilter->FullyConnectedOn();
            openingFilter->Update();
            cout << "Opening" << endl;
            thresholdImage = openingFilter->GetOutput();
        }*/
        cout << "Thresholding" << endl;

        // Fixed seed
		{
        	InImageType::IndexType start;
	        InImageType::SizeType size;
	        InImageType::RegionType centerRegion;

        	size.Fill(3);
	        start[0] = roiCenter[0] - 1;
	        start[1] = roiCenter[1] - 1;
	        start[2] = roiCenter[2] - 1;
	        centerRegion.SetSize(size);
	        centerRegion.SetIndex(start);
	        
            itk::ImageRegionIterator<OutImageType> centerPixels(thresholdImage, centerRegion);
			
            centerPixels.GoToBegin();
            while (!centerPixels.IsAtEnd()) {
				centerPixels.Set(1);

                ++centerPixels;
            }
        }

		/////////////////////////////////////////////////////////////////
		cout << "Connected Component" << endl;
		double target = 0;
		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
		try
		{
			AreaOpeningImageFilterType::Pointer openingFilter = AreaOpeningImageFilterType::New(); // delete small particles
			openingFilter->SetInput(thresholdImage);
			openingFilter->SetLambda(3);
			openingFilter->FullyConnectedOn();
			openingFilter->Update();
			cout << "Opening" << endl;
			
			connected->SetInput(openingFilter->GetOutput());
			connected->Update();


			/////////////////////////////////////////////////////////////////
			cout << "Select Center Component" << endl;
			unsigned int count[255] = { 0, };
			{
				itk::ImageRegionIterator<OutImageType> centerPixels(connected->GetOutput(), centerRegion);
				itk::ImageRegionIterator<OutImageType> centerPixels1(thresholdFilter_bg->GetOutput(), centerRegion);

				centerPixels.GoToBegin();
				centerPixels1.GoToBegin();
				while (!centerPixels.IsAtEnd()) {
					unsigned char label = centerPixels.Get();
					if (label > 0) {
						++count[label];
						centerPixels1.Set(0);
					}

					++centerPixels;
					++centerPixels1;
				}
			}
			unsigned int max_count = 0;
			for (unsigned c = 1; c < 255; c++) {
				if (count[c] > max_count) {
					max_count = count[c];
					target = c;
					cout << c << ", " << max_count << endl;
				}
			}
			cout << target << endl;

		}
		catch (itk::ExceptionObject &excp)
		{
			cerr << "Exception thrown while writing the series " << endl;
			cerr << excp << endl;
		}

        OutImageType::Pointer targetImage = OutImageType::New();
        int numSeeds = 0;

        targetImage->CopyInformation(inImageIso);
        targetImage->SetBufferedRegion(inImageIso->GetBufferedRegion());
        targetImage->Allocate();
        targetImage->FillBuffer(0);
        {
            itk::ImageRegionIterator<OutImageType> label(connected->GetOutput(), centerRegion);
            itk::ImageRegionIterator<OutImageType> targetLabel(targetImage, centerRegion);
            label.GoToBegin();
            targetLabel.GoToBegin();
            while (!label.IsAtEnd()) {
                if (label.Get() == target) {
                    targetLabel.Set(1);
                    ++numSeeds;
                } else {
					targetLabel.Set(0); 
				}
                ++label;
                ++targetLabel;
            }
        }

        cout << numSeeds << endl;



        /*{
            SubFilterType::Pointer subFilter = SubFilterType::New();
            subFilter->SetInput1(targetImage);
            subFilter->SetInput2(thresholdFilter_bg->GetOutput());
            subFilter->Update();

            targetImage = subFilter->GetOutput();

            itk::ImageRegionIterator<OutImageType> targetLabel(targetImage, targetImage->GetBufferedRegion());
            targetLabel.GoToBegin();
            while (!targetLabel.IsAtEnd()) {
                if (targetLabel.Get() != 1) targetLabel.Set(0);
                ++targetLabel;
            }

            AreaOpeningImageFilterType::Pointer openingFilter = AreaOpeningImageFilterType::New();
            openingFilter->SetInput(targetImage);
            openingFilter->SetLambda(5);
            openingFilter->Update();
            targetImage = openingFilter->GetOutput();
        }*/


        ErodeFilterType::Pointer erodeFilterB = ErodeFilterType::New();
        erodeFilterB->SetInput(thresholdFilter_bg->GetOutput());
        erodeFilterB->SetKernel(structuringElement1);
        erodeFilterB->SetErodeValue(1);
        erodeFilterB->Update();


        dilateFilter10->SetInput(targetImage);
        dilateFilter10->SetKernel(structuringElement10);
        dilateFilter10->SetDilateValue(1);
        dilateFilter10->Update();

        dilateFilter15->SetInput(targetImage);
        dilateFilter15->SetKernel(structuringElement15);
        dilateFilter15->SetDilateValue(1);
        dilateFilter15->Update();
        {
            ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
            connected->SetInput(dilateFilter15->GetOutput());
            connected->Update();

            OutImageType::Pointer image = connected->GetOutput();

            itk::ImageRegionIterator<OutImageType> centerLabels(image, centerRegion);

            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(image);
            duplicator->Update();
            OutImageType::Pointer tempImage = duplicator->GetOutput();


            double target = connected->GetOutput()->GetPixel(roiCenter);

            OutImageType::SizeType roiSize;
            roiSize[0] = static_cast<unsigned int>(longD / outSpacing[0]);
            roiSize[1] = static_cast<unsigned int>(longD / outSpacing[1]);
            roiSize[2] = static_cast<unsigned int>(longD / outSpacing[2]);

            OutImageType::IndexType roiStart;
            OutImageType::IndexType roiEnd;

            roiStart[0] = roiCenter[0] - roiSize[0] / 2;
            roiStart[1] = roiCenter[1] - roiSize[1] / 2;
            roiStart[2] = roiCenter[2] - roiSize[2] / 2;
            roiEnd[0]   = roiStart[0] + roiSize[0];
            roiEnd[1]   = roiStart[1] + roiSize[1];
            roiEnd[2]   = roiStart[2] + roiSize[2];

            cout << roiSize << endl;
            cout << roiStart << endl;
            cout << roiEnd << endl;


            itk::ImageRegionIterator<OutImageType> tempLabel(tempImage, tempImage->GetBufferedRegion());
            itk::ImageRegionIteratorWithIndex<OutImageType> label(dilateFilter15->GetOutput(),
              image->GetBufferedRegion());
            label.GoToBegin();
            tempLabel.GoToBegin();
            while (!label.IsAtEnd()) {
                OutImageType::IndexType idx = label.GetIndex();

                if (idx[0] >= roiStart[0] && idx[1] >= roiStart[1] && idx[2] >= roiStart[2] &&
                  idx[0] < roiEnd[0] && idx[1] < roiEnd[1] && idx[2] < roiEnd[2])
                {
                    if (tempLabel.Get() == target) label.Set(1);
                } else { 
					label.Set(0); 
				}
                ++label;
                ++tempLabel;
            }
        }


        // generate seed areas
        {
            // itk::ImageRegionIterator< OutImageType > bg(erodeFilterB->GetOutput(), erodeFilterB->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<OutImageType> bg(thresholdFilter_bg->GetOutput(),
              thresholdFilter_bg->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<OutImageType> gray(dilateFilter10->GetOutput(),
              dilateFilter10->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<OutImageType> truebg(dilateFilter15->GetOutput(),
              dilateFilter15->GetOutput()->GetBufferedRegion());

            itk::ImageRegionIterator<OutImageType> targetLabel(targetImage, targetImage->GetBufferedRegion());
            itk::ImageRegionIterator<OutImageType> newLabel(labImage, labImage->GetBufferedRegion());
            itk::ImageRegionIterator<WeightImageType> weight(wtImage, wtImage->GetBufferedRegion());

            targetLabel.GoToBegin();
            newLabel.GoToBegin();
            weight.GoToBegin();

            bg.GoToBegin();
            gray.GoToBegin();
            truebg.GoToBegin();
            while (!targetLabel.IsAtEnd()) {
                weight.Set(0.0);

                if (targetLabel.Get() == 1) { // taget
                    weight.Set(contrastNoiseRatio);
                    newLabel.Set(1);
                }

				if (bg.Get() == 1) { // background
                    weight.Set(contrastNoiseRatio);
                    newLabel.Set(2);
                }
				else if (gray.Get() == 0) {
					weight.Set(contrastNoiseRatio);
					newLabel.Set(3);
				}

                if (truebg.Get() == 0) { // definetely background, opposite value (0 is bg)
                    //weight.Set(contrastNoiseRatio);
                    weight.Set(priorSegmentStrength);
                    newLabel.Set(4);
                }

                ++targetLabel;
                ++newLabel;
                ++weight;
                ++bg;
                ++gray;
                ++truebg;
            }
        }

        {
            OutImageWriterType::Pointer writer = OutImageWriterType::New();
            writer->SetUseCompression(true);
            writer->SetInput(labImage);
            writer->SetFileName(string(argv[7]) + "Tc-label.nrrd");
            writer->Update();
        }
    }

    cout << "label image generation completed" << endl;


    /////////////////////////////////////////////////////////////////

    // GrowCut Segmentation
    GrowCutFilterType::Pointer filter = GrowCutFilterType::New();

    // Segmentation filter parameters
    filter->SetInput(inImageIso);
    filter->SetLabelImage(labImage);

    filter->SetStrengthImage(wtImage);

    filter->SetSeedStrength(contrastNoiseRatio);
    filter->SetObjectRadius((unsigned int) objectSize);

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
        itk::ImageRegionIterator<OutImageType> label(outputImageROI, outputImageROI->GetBufferedRegion());
        itk::ImageRegionIterator<OutImageType> truebg(dilateFilter15->GetOutput(), dilateFilter15->GetOutput()->GetBufferedRegion());

        label.GoToBegin();
        truebg.GoToBegin();
        while (!label.IsAtEnd()) {
            if (label.Get() != 1) label.Set(0);
            if (truebg.Get() == 0) label.Set(0);
            ++label;
            ++truebg;
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
        cout << "Opening" << endl;
        outputImageROI = openingFilter->GetOutput();
    }

    /*{
     *  OutImageWriterType::Pointer writer = OutImageWriterType::New();
     *  writer->SetUseCompression(true);
     *  writer->SetInput(outputImageROI);
     *  writer->SetFileName(string(argv[7]) + "mask-label.nrrd");
     *  writer->Update();
     * }*/


    /////////////////////////////////////////////////////////////////

    /*OutImageType::Pointer outputImage = OutImageType::New();
     * // allocate outputImage first
     * outputImage->CopyInformation(image);
     * outputImage->SetSpacing(spacing);
     * outputImage->SetBufferedRegion(image->GetBufferedRegion());
     * outputImage->Allocate();
     * outputImage->FillBuffer(0);
     *
     * itk::ImageRegionIterator< OutImageType > filterOut(outputImageROI, outputImageROI->GetBufferedRegion());
     * itk::ImageRegionIterator< OutImageType > out(outputImage, oRegion);
     *
     * for (filterOut.GoToBegin(), out.GoToBegin(); !filterOut.IsAtEnd(); ++filterOut, ++out)
     * {
     *  if (filterOut.Get() == 1)
     *      out.Set(1);
     *  //out.Set(filterOut.Get());
     * }
     *
     * cout << "writeImage_spacing = " << outputImage->GetSpacing() << endl;
     * cout << "writeImage_origin  = " << outputImage->GetOrigin() << endl;
     * cout << "writeImage_LargestPossibleRegion = " << outputImage->GetLargestPossibleRegion() << endl;*/


    ////////To write Output images /////////////////////
    OutImageWriterType::Pointer OutImageWriter = OutImageWriterType::New();
    OutImageWriter->SetUseCompression(true);
    OutImageWriter->SetInput(outputImageROI);
    OutImageWriter->SetFileName(argv[7]);

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

    // Vesselness feature generate
    {
        VesselnessGeneratorType::Pointer vesselnessFeatureGenerator = VesselnessGeneratorType::New();
        InputImageSpatialObjectType::Pointer inputSpatialObject     = InputImageSpatialObjectType::New();

        using CastImageFilterType = itk::CastImageFilter<WeightImageType, InImageType>;
        CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
        castImageFilter->SetInput(inImageIso);

        inputSpatialObject->SetImage(castImageFilter->GetOutput());

        vesselnessFeatureGenerator->SetInput(inputSpatialObject);
        vesselnessFeatureGenerator->SetSigma(1.0);
        vesselnessFeatureGenerator->SetAlpha1(0.1);
        vesselnessFeatureGenerator->SetAlpha2(2.0);
        vesselnessFeatureGenerator->SetSigmoidAlpha(-10.0);
        vesselnessFeatureGenerator->SetSigmoidBeta(40.0);
        vesselnessFeatureGenerator->SetUseVesselEnhancingDiffusion(true);

        vesselnessFeatureGenerator->Update();

        cout << "vesselness feature generation completed" << endl;

        typename SpatialObjectType::Pointer vesselness =
          const_cast<SpatialObjectType *>(vesselnessFeatureGenerator->GetFeature());
        typename OutputImageSpatialObjectType::Pointer outputObject =
          dynamic_cast<OutputImageSpatialObjectType *>( vesselness.GetPointer() );
        typename WeightImageType::Pointer vesselImage =
          const_cast<WeightImageType *>(outputObject->GetImage());
        vesselImage->DisconnectPipeline();

    itk::ImageRegionIterator<WeightImageType> vesselnessit(vesselImage, vesselImage->GetBufferedRegion());
    vesselnessit.GoToBegin();
    while (!vesselnessit.IsAtEnd()) {
      vesselnessit.Set(1 - vesselnessit.Get());

      ++vesselnessit;
    }

        {
            WeightImageWriterType::Pointer writer = WeightImageWriterType::New();
            writer->SetUseCompression(true);
            writer->SetInput(vesselImage);
            writer->SetFileName(string(argv[7]) + "-vessel.nrrd");
            writer->Update();
        }
    }

    return EXIT_SUCCESS;
} // main
