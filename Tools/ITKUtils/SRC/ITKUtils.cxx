#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "ITKUtils.h"

MaskImageType::Pointer GetMaskImage(LabelImageType::Pointer labelImage, LabelPixelType selectedLabelValue)
{
    MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->SetRegions( labelImage->GetRequestedRegion() );
    maskImage->CopyInformation( labelImage );
    maskImage->Allocate();

    itk::ImageRegionIterator< LabelImageType > label(labelImage, labelImage->GetBufferedRegion() );
    itk::ImageRegionIterator< MaskImageType > mask(maskImage, maskImage->GetBufferedRegion() );

    for (label.GoToBegin(), mask.GoToBegin(); !label.IsAtEnd(); ++label, ++mask)
    {
        LabelImageType::PixelType color = label.Get();
        if (color == selectedLabelValue)
        {
            mask.Set(maskValue);
        }
        else
        {
            mask.Set(0);
        }
    }

    return maskImage;
}

MaskImageType::RegionType RoiIndexToRegion(MaskImageType::IndexType roiStart, MaskImageType::IndexType roiEnd)
{
    MaskImageType::RegionType outputReigion;
    MaskImageType::IndexType outputStart;
    MaskImageType::SizeType outputSize;

    for (unsigned n = 0; n < Dimension; n++)
    {
        outputStart[n] = roiStart[n];
        outputSize[n] = roiEnd[n] - roiStart[n] + 1;
    }

    outputReigion.SetIndex(outputStart);
    outputReigion.SetSize(outputSize);

    return outputReigion;
}

void RegionToIndex(MaskImageType::RegionType outputReigion, MaskImageType::IndexType &roiStart, MaskImageType::IndexType &roiEnd)
{
    MaskImageType::IndexType outputStart = outputReigion.GetIndex();
    MaskImageType::SizeType outputSize = outputReigion.GetSize();

    for (unsigned n = 0; n < Dimension; n++)
    {
        roiStart[n] = outputStart[n];
        roiEnd[n] = outputStart[n] + outputSize[n] - 1;
    }
}

MaskImageType::RegionType ReadROI(string roiName)
{
    MaskImageType::RegionType roiRegion;
    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;
    MaskImageType::SizeType roiSize;

    ifstream roiFile;
    char buf[256];
    roiFile.open(roiName.c_str(), ios::in);

    roiFile >> roiStart[0];roiFile >> roiStart[1];roiFile >> roiStart[2];
    roiFile.getline(buf,256);
    cout << roiStart << endl;

    roiFile >> roiEnd[0];roiFile >> roiEnd[1];roiFile >> roiEnd[2];
    roiFile.getline(buf,256);
    cout << roiEnd << endl;

    roiFile >> roiSize[0];roiFile >> roiSize[1];roiFile >> roiSize[2];
    roiFile.getline(buf,256);
    cout << roiSize << endl;

    roiRegion.SetIndex(roiStart);
    roiRegion.SetSize(roiSize);

    roiFile.close();

    return roiRegion;
}

void WriteROI(MaskImageType::RegionType inputRegion, string roiName)
{
    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;
    RegionToIndex(inputRegion,roiStart,roiEnd);
    ofstream roiFile;
    
    roiFile.open(roiName.c_str(), ios::out);
    for (unsigned i = 0; i < Dimension; i++)
        roiFile << roiStart[i] << " ";
    roiFile << endl;
    for (unsigned i = 0; i < Dimension; i++)
        roiFile << roiEnd[i] << " ";
    roiFile << endl;
    for (unsigned i = 0; i < Dimension; i++)
        roiFile << inputRegion.GetSize()[i] << " ";
    roiFile << endl;
    roiFile.close();
}

MaskImageType::RegionType GetRoi(MaskImageType::Pointer maskImage)
{
    MaskImageType::RegionType roiRegion;
    MaskImageType::IndexType roiStart;
    MaskImageType::IndexType roiEnd;

    roiStart[0] = 0; roiStart[1] = 0; roiStart[2] = 0;
    roiEnd[0] = 0; roiEnd[1] = 0; roiEnd[2] = 0;

    bool foundMask = false;

    itk::ImageRegionIteratorWithIndex< MaskImageType > mask(maskImage, maskImage->GetBufferedRegion() );
    for (mask.GoToBegin(); !mask.IsAtEnd(); ++mask)
    {
        MaskImageType::PixelType color = mask.Get();
        if (color > 0)
        {
            MaskImageType::IndexType idx = mask.GetIndex();
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (!foundMask)
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
            foundMask = true;
        }
    }

    cout << "roiStart " << roiStart << ", roiEnd " << roiEnd << "" << endl;

    roiRegion = RoiIndexToRegion(roiStart, roiEnd);

    return roiRegion;
}

void ExpandRoi(MaskImageType::Pointer maskImage, MaskImageType::RegionType &maskRegion,const double objectSize)
{
    RegionType maskImageRegion = maskImage->GetLargestPossibleRegion();
    SizeType maskImageSize = maskImageRegion.GetSize();

    MaskImageType::IndexType maskStart = maskRegion.GetIndex();
    MaskImageType::SizeType maskSize = maskRegion.GetSize();

    MaskImageType::IndexType index = maskImage->GetBufferedRegion().GetIndex();
    MaskImageType::SizeType size = maskImage->GetBufferedRegion().GetSize();
    MaskImageType::SpacingType spacing = maskImage->GetSpacing();


    for (unsigned i = 0; i < Dimension; i++)
    {
        MaskImageType::PixelType radius = static_cast< MaskImageType::PixelType> (objectSize/spacing[i]);
        int diff = static_cast< int > (maskStart[i] - radius);
        if (diff >= index[i])
        {
            maskStart[i] -= radius;
            maskSize[i] += 2*radius;
        }
        else
        {
            maskStart[i] = index[i];
        }
        maskStart[i] = (maskStart[i] < 0) ? 0 : maskStart[i];
        maskSize[i] = (static_cast<unsigned int>(maskStart[i] + maskSize[i] - 1) < maskImageSize[i]) ?
                      maskSize[i] : (maskImageSize[i] - maskStart[i]);
    }

    cout << "maskStart " << maskStart << ", maskEnd " << maskStart + maskSize << "" << endl;

    maskRegion.SetIndex(maskStart);
    maskRegion.SetSize(maskSize);
}

void BoundingCheck(MaskImageType::Pointer maskImage, InputImageType::Pointer inputImage, MaskImageType::RegionType &maskRegion, InputImageType::RegionType &inputRegion)
{
    RegionType inputImageRegion = inputImage->GetLargestPossibleRegion();
    SizeType inputImageSize = inputImageRegion.GetSize();
    
    MaskImageType::IndexType maskStart = maskRegion.GetIndex();
    MaskImageType::SizeType maskSize = maskRegion.GetSize();
    MaskImageType::IndexType maskCenter;

    InputImageType::IndexType istart;
    InputImageType::SizeType isize;

    MaskImageType::IndexType ostart;
    MaskImageType::SizeType osize;

    MaskImageType::PointType ptRoiStart;
    maskImage->TransformIndexToPhysicalPoint(maskStart, ptRoiStart);
    inputImage->TransformPhysicalPointToIndex(ptRoiStart, istart);

    for (unsigned n = 0; n < Dimension; n++)
    {
        isize[n] = maskSize[n];
        ostart[n] = maskStart[n];

        // Input image
        if (istart[n] < 0)
        {
            ostart[n] = ostart[n] - istart[n];
            isize[n] = isize[n] + istart[n];
            istart[n] = 0;
        }
        isize[n] = (static_cast<unsigned int>(isize[n] + istart[n]) <= inputImageSize[n]) ?
                   isize[n] : (inputImageSize[n] - istart[n] - 1);
        // Mask
        osize[n] = isize[n];
    
        maskCenter[n] = static_cast<OutputPixelType>(round(isize[n] / 2));
    }

    cout << "maskStart " << maskStart << ", maskEnd " << maskStart + maskSize << "" << endl;
    cout << "istart " << istart << ", isize " << isize << "" << endl;
    cout << "ostart " << ostart << ", osize " << osize << "" << endl;
    cout << "maskCenter " << maskCenter << endl;

    inputRegion.SetIndex(istart);
    inputRegion.SetSize(isize);

    maskRegion.SetIndex(ostart);
    maskRegion.SetSize(osize);
}

MaskImageType::Pointer unionImages(vector<MaskImageType::Pointer> inputMaskImages)
{
    MaskImageType::Pointer outputMaskImage = MaskImageType::New();
    outputMaskImage->SetRegions( inputMaskImages[0]->GetLargestPossibleRegion() );
    outputMaskImage->CopyInformation( inputMaskImages[0] );
    outputMaskImage->Allocate();
    outputMaskImage->FillBuffer(0);

    int volume = 0;
    for (int i = 0; i < inputMaskImages.size(); i++)
    {   
        itk::ImageRegionIterator< MaskImageType > in(inputMaskImages[i], inputMaskImages[i]->GetBufferedRegion() );
        itk::ImageRegionIterator< MaskImageType > out(outputMaskImage, outputMaskImage->GetBufferedRegion() );

        for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd(); ++in, ++out)
        {
            MaskImageType::PixelType mask = in.Get();
            if (mask > 0)
            {
                out.Set(maskValue);
            }
        }
    }
 
    return outputMaskImage;
}
