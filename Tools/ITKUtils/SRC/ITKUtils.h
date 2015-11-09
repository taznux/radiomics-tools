#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include <itkCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRegionOfInterestImageFilter.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

#include "itkMaskedSTAPLEImageFilter.h"

#include <string>
using namespace std;

// Definition of Data type
typedef short InputPixelType;
typedef float InternalPixelType;
typedef short OutputPixelType;

typedef short LabelPixelType;
typedef short MaskPixelType;

const unsigned int Dimension = 3;
const unsigned int OutputDimension = 2;
const MaskPixelType maskValue = 1;

typedef itk::Image < InputPixelType, Dimension >  InputImageType;
typedef itk::Image < InternalPixelType, Dimension > InternalImageType;
typedef itk::Image < OutputPixelType, OutputDimension > OutputImageType;
typedef itk::Image < InputPixelType, Dimension - 1 > InputImage2DType;

typedef itk::Image < LabelPixelType, Dimension> LabelImageType;
typedef itk::Image < MaskPixelType, Dimension> MaskImageType;
typedef itk::Image < MaskPixelType, Dimension - 1 > MaskImage2DType;


// Definition of Properties
typedef InputImageType::SpacingType    SpacingType;
typedef InputImageType::PointType      OriginType;
typedef InputImageType::RegionType     RegionType;
typedef InputImageType::SizeType       SizeType;
typedef InputImageType::IndexType      IndexType;


template <typename TImageType>
typename TImageType::Pointer ReadImageFile(string fileName)
{
    typedef itk::ImageFileReader<TImageType> ImageReaderType;

    typename ImageReaderType::Pointer ImageReader = ImageReaderType::New();
    ImageReader->SetFileName(fileName);
    ImageReader->Update();

    itk::SmartPointer<TImageType> Image = ImageReader->GetOutput();

    return Image;
}


template <typename TImageType>
void WriteImageFile(typename TImageType::Pointer outputImage, string fileName)
{
    ////////To write Output images /////////////////////
    typedef itk::ImageFileWriter<TImageType> ImageWriterType;

    typename ImageWriterType::Pointer OutImageWriter = ImageWriterType::New();
    OutImageWriter->SetUseCompression(true);
    OutImageWriter->SetInput(outputImage);
    OutImageWriter->SetFileName( fileName );

    try
    {
        OutImageWriter->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
        std::cerr << "Exception thrown while writing the series " << std::endl;
        std::cerr << excp << std::endl;
        //return EXIT_FAILURE;
    }
}


template <typename TImageType>
typename TImageType::Pointer ApplyRoi(typename TImageType::Pointer inputImage, itk::ImageRegion<Dimension> inputRegion)
{
    typename TImageType::Pointer roiImage;

    typedef itk::RegionOfInterestImageFilter<TImageType, TImageType> RoiFilterType;
    typename RoiFilterType::Pointer RoiFilter = RoiFilterType::New();
    RoiFilter->SetRegionOfInterest(inputRegion);

    RoiFilter->SetInput( inputImage );
    RoiFilter->Update();

    roiImage = RoiFilter->GetOutput();

    return roiImage;
 }

MaskImageType::Pointer GetMaskImage(LabelImageType::Pointer, LabelPixelType);

MaskImageType::RegionType RoiIndexToRegion(MaskImageType::IndexType, MaskImageType::IndexType);
void RegionToIndex(MaskImageType::RegionType, MaskImageType::IndexType &, MaskImageType::IndexType &);

MaskImageType::RegionType GetRoi(MaskImageType::Pointer);
void ExpandRoi(MaskImageType::Pointer, MaskImageType::RegionType &, const double objectSize = 10);
void BoundingCheck(MaskImageType::Pointer, InputImageType::Pointer, MaskImageType::RegionType &, InputImageType::RegionType &);

MaskImageType::RegionType ReadROI(string);
void WriteROI(MaskImageType::RegionType, string);


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
                  << static_cast<unsigned int>( 100.0 * m_Process->GetProgress() ) << "%";
    }
    itk::ProcessObject::Pointer m_Process;
};
