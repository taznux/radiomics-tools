#if defined(_MSC_VER)
# pragma warning(disable : 4786)
#endif

#include <itkCastImageFilter.h>
#include <itkCommand.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>

#include <itkExtractImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>

// Operations on Images

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <vnl/vnl_math.h>

// For Feature
#include <itkRescaleIntensityImageFilter.h>

#include <itkMinimumMaximumImageCalculator.h>
#include <itkStatisticsImageFilter.h>

#include <itkScalarImageToRunLengthFeaturesFilter.h>
#include <itkScalarImageToTextureFeaturesFilter.h>

#include <itkConnectedComponentImageFilter.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

// For Label Object Representation
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkBinaryImageToStatisticsLabelMapFilter.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMap.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelObject.h>
#include <itkShapeLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkStatisticsLabelMapFilter.h>
#include <itkStatisticsLabelObject.h>

// For island removing
#include <itkBinaryShapeKeepNObjectsImageFilter.h>

#include <itkAndImageFilter.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkSliceBySliceImageFilter.h>
#include <itkSubtractImageFilter.h>

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkFlatStructuringElement.h"

#include <itkMaskedSTAPLEImageFilter.h>
#include <itkSTAPLEImageFilter.h>

// For file operation
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// For threshold
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>

using namespace std;

#include "ITKUtils.h"

int
main(int argc, char * argv[])
{
    if (argc < 4) {
        cerr << "Missing Parameters " << endl;
        cerr << "Usage = " << argv[0];
        cerr << " OutputImage OutputText Label InputLabelImage1 InputLabelImage2 "
            "InputLabelImage3 ..."
             << endl;
        return EXIT_FAILURE;
    }

    int numInputs = argc - 4;

    cout << "The number of input imeages = " << numInputs << endl;

    string outputImageName   = argv[1];
    string outputTxtFileName = argv[2];

    // Output Text file open
    ofstream outfile((outputTxtFileName).c_str(), ios::out);
    if (!outfile) {
        cerr << "Can't open output file " << outputTxtFileName << endl;
        exit(1);
    }

    LabelPixelType selectedLabelValue = atoi(argv[3]);

    vector<string> inputLabelImageNames;
    cout << "Output Image Name = " << outputImageName << endl;
    for (int i = 0; i < numInputs; i++) {
        inputLabelImageNames.push_back(argv[i + 4]);
        cout << "Input Label Image " << i << " Name = " << inputLabelImageNames[i]
             << endl;
    }

    ///////// To Read the Mask image ////////////////////
    cout << "Get the mask from Input" << endl;
    vector<MaskImageType::Pointer> inputMaskImages;
    for (int i = 0; i < numInputs; i++) {
        MaskImageType::Pointer maskImage;
        if (selectedLabelValue > 0) {
            maskImage =
              GetMaskImage(ReadImageFile<LabelImageType>(inputLabelImageNames[i]),
                selectedLabelValue);
        } else {
            maskImage =
              GetMaskImage(ReadImageFile<LabelImageType>(inputLabelImageNames[i]));
        }

        inputMaskImages.push_back(maskImage);
    }

    // To check coordinate of the input label image
    SpacingType inputImageSpacing = inputMaskImages[0]->GetSpacing();
    OriginType inputImageOrigin   = inputMaskImages[0]->GetOrigin();
    RegionType inputImageRegion   = inputMaskImages[0]->GetLargestPossibleRegion();
    SizeType inputImageSize       = inputImageRegion.GetSize();

    /////////////////////////////

    cout << "Input Image Spacing = " << inputImageSpacing << endl;
    cout << "Input Image Origin = " << inputImageOrigin << endl;
    cout << "Input Image Size = " << inputImageSize << endl << endl << endl;

    // union maskImages
    MaskImageType::Pointer unionMaskImage = unionImages(inputMaskImages);

    RegionType roiRegion;
    IndexType roiStart;
    IndexType roiEnd;

    roiRegion = GetRoi(unionMaskImage);
    ExpandRoi(unionMaskImage, roiRegion, 4);
    RegionToIndex(roiRegion, roiStart, roiEnd);

    unionMaskImage = ApplyRoi<MaskImageType>(unionMaskImage, roiRegion);

    InternalImageType::Pointer outputImage;
    typedef itk::STAPLEImageFilter<MaskImageType, InternalImageType>
      StapleFilterType;

    StapleFilterType::Pointer sFilter = StapleFilterType::New();
    for (int i = 0; i < inputMaskImages.size(); i++) {
        inputMaskImages[i] = ApplyRoi<MaskImageType>(inputMaskImages[i], roiRegion);
        sFilter->SetInput(i, inputMaskImages[i]);
    }
    // sFilter->SetMaskImage(unionMaskImage);
    // sFilter->SetMaskValue(maskValue);
    sFilter->SetForegroundValue(maskValue);
    sFilter->SetCoordinateTolerance(0.001);
    sFilter->SetDirectionTolerance(0.001);

    ShowProgressObject progressWatch(sFilter);
    itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
    command = itk::SimpleMemberCommand<ShowProgressObject>::New();
    command->SetCallbackFunction(&progressWatch,
      &ShowProgressObject::ShowProgress);
    sFilter->AddObserver(itk::ProgressEvent(), command);

    sFilter->Update();
    outputImage = sFilter->GetOutput();

    #if 0
    priority_queue<float> probs;
    itk::ImageRegionIterator<InternalImageType> output(outputImage, outputImage->GetBufferedRegion());
    for (output.GoToBegin(); !output.IsAtEnd(); ++output) {
        InternalImageType::PixelType prob = output.Get();
        probs.push(prob);
    }
    int volumeG      = 0;
    int volumeU      = 0;
    int error        = 1000000000;
    float threshold  = 0;
    float pthreshold = 0;
    for (int count = 0; probs.size() > 0; count++) {
        float new_threshold = probs.top();
        if (new_threshold >= 0.5) {
            volumeG++;
        } else if (pthreshold != new_threshold)   {
            int new_error = abs(count - 2 * volumeG);
            if (new_error < error) {
                error     = new_error;
                threshold = new_threshold;
                volumeU   = count;
                cout << "Error " << error << ", T " << threshold << ", Volume" << count << endl;
            }
            pthreshold = new_threshold;
        }
        probs.pop();
    }


    int volumeU1 = 0;
    itk::ImageRegionIterator<MaskImageType> mask(unionMaskImage, unionMaskImage->GetBufferedRegion());
    for (mask.GoToBegin(); !mask.IsAtEnd(); ++mask) {
        MaskImageType::PixelType m = mask.Get();
        if (m > 0)
            volumeU1++;
    }
    cout << "Threshold is " << threshold << endl;
    cout << "G Volume is " << volumeG << endl;
    cout << "U Volume is " << volumeU << endl;
    cout << "U1 Volume is " << volumeU1 << endl;
    #endif // if 0
    #if 0
    typedef itk::Image<MaskPixelType, Dimension - 1> MaskImage2DType;
    typedef itk::FlatStructuringElement<Dimension - 1> StructuringElement2DType;
    typedef itk::BinaryDilateImageFilter<MaskImage2DType, MaskImage2DType, StructuringElement2DType> DilateFilter2DType;
    typedef itk::BinaryErodeImageFilter<MaskImage2DType, MaskImage2DType, StructuringElement2DType> ErodeFilter2DType;
    typedef itk::SliceBySliceImageFilter<MaskImageType, MaskImageType> SliceBySliceFilterType;

    StructuringElement2DType::RadiusType elementRadius2D;
    elementRadius2D.Fill(1);

    StructuringElement2DType structuringElement2D = StructuringElement2DType::Ball(elementRadius2D, true);

    int volumeU = 0;
    int volumeG = 0;
    itk::ImageRegionIterator<InternalImageType> output(outputImage, outputImage->GetBufferedRegion());
    for (output.GoToBegin(); !output.IsAtEnd(); ++output) {
        InternalImageType::PixelType prob = output.Get();
        if (prob > 0.5) {
            volumeG++;
        }
    }
    while (volumeU < volumeG * 2) {
        itk::ImageRegionIterator<MaskImageType> mask(unionMaskImage, unionMaskImage->GetBufferedRegion());

        for (mask.GoToBegin(); !mask.IsAtEnd(); ++mask) {
            MaskImageType::PixelType m = mask.Get();
            if (m > 0)
                volumeU++;
        }

        cout << "G Volume is " << volumeG << endl;
        cout << "U Volume is " << volumeU << endl;

        if (volumeU >= volumeG * 2) break;

        SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();

        DilateFilter2DType::Pointer dilate2DFilter = DilateFilter2DType::New();
        dilate2DFilter->SetKernel(structuringElement2D);
        dilate2DFilter->SetDilateValue(maskValue);

        sliceBySliceFilter->SetInput(unionMaskImage);
        sliceBySliceFilter->SetFilter(dilate2DFilter);
        sliceBySliceFilter->Update();
        unionMaskImage = sliceBySliceFilter->GetOutput();
    }
    #endif // if 0
    double avg_p = 0.0;
    double avg_q = 0.0;
    double avg_j = 0.0;
    // Print out the specificities
    cout << "Number of elapsed iterations = " << sFilter->GetElapsedIterations()
         << endl;

    //  cout.precision(5);
    //  cout.setf(ios::fixed, ios::floatfield);

    cout << "File \t\t"
         << "\tSensitivity(p) "
         << "\tSpecificity(q)"
         << "\tJaccard Index(j)" << endl;
    cout << "-----\t\t"
         << "\t-------------- "
         << "\t--------------"
         << "\t--------------" << endl;
    outfile << "File "
            << "\tSensitivity(p) "
            << "\tSpecificity(q)"
            << "\tJaccard Index(j)" << endl;
    itk::ImageRegionIterator<InternalImageType> output(
        outputImage, outputImage->GetBufferedRegion());
    itk::ImageRegionIterator<MaskImageType> unionMask(
        unionMaskImage, unionMaskImage->GetBufferedRegion());

    double R  = 0;
    double G  = 0;
    double cG = 0;
    double U  = 0;

    for (output.GoToBegin(), unionMask.GoToBegin();
      !output.IsAtEnd() && !unionMask.IsAtEnd(); ++output, ++unionMask)
    {
        InternalImageType::PixelType prob = output.Get();
        char u = unionMask.Get() > 0;
        char g = prob >= 0.5;
        G  += g;
        cG += (!g) && u;
        U  += u;
        R  += 1;
    }

    for (int i = 0; i < numInputs; i++) {
        double sensitivity = 0;
        double specificity = 0;
        double jaccard     = 0;
        double GD   = 0;
        double GuD  = 0;
        double D    = 0;
        double cGcD = 0;

        itk::ImageRegionIterator<MaskImageType> individualMask(
            inputMaskImages[i], inputMaskImages[i]->GetBufferedRegion());
        for (output.GoToBegin(), unionMask.GoToBegin(), individualMask.GoToBegin();
          !output.IsAtEnd() && !unionMask.IsAtEnd() && !individualMask.IsAtEnd();
          ++output, ++unionMask, ++individualMask)
        {
            InternalImageType::PixelType prob = output.Get();
            char u = unionMask.Get() > 0;
            char d = individualMask.Get() > 0;
            char g = prob >= 0.5;
            GD   += (g && d);
            GuD  += (g || d);
            cGcD += (!g && !d) && u;
            D    += d;
        }
        double padding = 2 * G - U;
        if (padding < 0)
            padding = 0;
        double G2 = U + padding;
        cout << "R = " << R << ", "
             << "G = " << G << ", "
             << "!G = " << cG << ", "
             << "D = " << D << ", "
             << "U = " << U << ", "
             << "G&D = " << GD << ", "
             << "G|D = " << GuD << ", "
             << "!G&!D = " << cGcD << ", "
             << "padding = " << padding << endl;
        sensitivity = sFilter->GetSensitivity(i);
        specificity = (sFilter->GetSpecificity(i) * (R - G) - R + G2) / (G2 - G);
        jaccard     = GD / GuD;

        string name = inputLabelImageNames[i];
        #ifdef WIN32
        size_t found = name.rfind("\\");
        #else
        size_t found = name.rfind("/");
        #endif

        if (found > 0) {
            name = inputLabelImageNames[i].substr(
                found + 1, inputLabelImageNames[i].length() - (found + 1)
                - 11); // post 11 -> -label.nrrd
        } else                                                                  {
            name = inputLabelImageNames[i].substr(
                0, inputLabelImageNames[i].length() - 11); // post 11 -> -label.nrrd
        }
        // string name = inputLabelImageNames[i];
        // avg_q += sFilter->GetSpecificity(i);
        // avg_p += sFilter->GetSensitivity(i);
        avg_p += sensitivity;
        avg_q += specificity;
        avg_j += jaccard;
        cout << i << ": " << name << "\t" << sFilter->GetSensitivity(i) << "\t"
             << sFilter->GetSpecificity(i) << endl;
        cout << i << ": " << name << "\t" << sensitivity << "\t" << specificity
             << "\t" << jaccard << endl;
        // outfile << name << "\t"
        //      << sFilter->GetSensitivity(i) << "\t"
        //      << sFilter->GetSpecificity(i) << endl;
        outfile << name << "\t" << sensitivity << "\t" << specificity << "\t"
                << jaccard << endl;
    }
    avg_p /= static_cast<double>(numInputs);
    avg_q /= static_cast<double>(numInputs);
    avg_j /= static_cast<double>(numInputs);

    cout << "Mean:\t" << avg_p << "\t" << avg_q << "\t" << avg_j << endl << endl;
    outfile << "Mean:\t" << avg_p << "\t" << avg_q << "\t" << avg_j << endl
            << endl;

    cout << "writeImage_spacing = " << outputImage->GetSpacing() << endl;
    cout << "writeImage_origin  = " << outputImage->GetOrigin() << endl;
    cout << "writeImage_LargestPossibleRegion = "
         << outputImage->GetLargestPossibleRegion() << endl;

    ////////To write Output images /////////////////////
    try {
        WriteImageFile<InternalImageType>(outputImage, outputImageName);
    } catch (itk::ExceptionObject &excp) {
        cerr << "Exception thrown while writing the series " << endl;
        cerr << excp << endl;
        return EXIT_FAILURE;
    }

    cout << "Saved the normal tissue and tumor boundary image" << endl;

    outfile.close();

    return EXIT_SUCCESS;
} // main
