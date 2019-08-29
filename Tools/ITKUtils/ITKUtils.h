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

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>


// std
#include <iostream>
#include <fstream>
#include <limits>


#include <csignal>

#include <string>
using namespace std;

// Definition of Data type
using InputPixelType = float;
using InternalPixelType = float;
using OutputPixelType = float;

using LabelPixelType = short;
using MaskPixelType = short;

constexpr unsigned int Dimension = 3;
constexpr unsigned int OutputDimension = 2;
constexpr MaskPixelType maskValue = 1;

using InputImageType = itk::Image < InputPixelType, Dimension >;
using InternalImageType = itk::Image < InternalPixelType, Dimension >;
using OutputImageType = itk::Image < OutputPixelType, OutputDimension >;
using InputImage2DType = itk::Image < InputPixelType, Dimension - 1 >;

using LabelImageType = itk::Image < LabelPixelType, Dimension>;
using MaskImageType = itk::Image < MaskPixelType, Dimension>;
using MaskImage2DType = itk::Image < MaskPixelType, Dimension - 1 >;


// Definition of Properties
using SpacingType = InputImageType::SpacingType;
using OriginType = InputImageType::PointType;
using RegionType = InputImageType::RegionType;
using SizeType = InputImageType::SizeType;
using IndexType = InputImageType::IndexType;
using DirectionType = InputImageType::DirectionType;


//Definition of Morphological operation sturcture elements
using StructuringElementType = itk::BinaryBallStructuringElement< MaskPixelType, Dimension >;
using DilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;
using ErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;


template <typename TImageType>
typename TImageType::Pointer ReadImageFile(string fileName)
{
  using ImageReaderType = itk::ImageFileReader<TImageType>;

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
  using ImageWriterType = itk::ImageFileWriter<TImageType>;

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

  using RoiFilterType = itk::RegionOfInterestImageFilter<TImageType, TImageType>;
  typename RoiFilterType::Pointer RoiFilter = RoiFilterType::New();
  RoiFilter->SetRegionOfInterest(inputRegion);

  RoiFilter->SetInput( inputImage );
  RoiFilter->Update();

  roiImage = RoiFilter->GetOutput();

  return roiImage;
}

MaskImageType::Pointer unionImages(vector<MaskImageType::Pointer> inputMaskImages);

namespace FGC
{

/*
* createImage
*/

template<typename itkImageType>
typename itkImageType::Pointer createImage(typename itkImageType::SizeType size, int iniValue)
{

  typename itkImageType::Pointer img = itkImageType::New();
  typename itkImageType::SpacingType spacing;
  typename itkImageType::IndexType start;
  start.Fill(0);
  spacing.Fill(1);
  spacing[itkImageType::ImageDimension - 1] = 2;

  typename itkImageType::RegionType region;
  img->SetSpacing(spacing);
  region.SetSize(size);
  region.SetIndex(start);
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(static_cast<typename itkImageType::PixelType>(iniValue));

  const typename itkImageType::SpacingType &sp = img->GetSpacing();
  std::cout << "Spacing = ";
  std::cout << sp[0] << ", " << sp[1] << std::endl;


  return img;

}


template<typename itkImageType>
void FindITKImageROI(typename itkImageType::Pointer &im, std::vector<long> &imROI)
{

  typename itkImageType::IndexType roiStart;
  typename itkImageType::IndexType roiEnd;
  typename itkImageType::IndexType start;
  typename itkImageType::SizeType size;

  size = im->GetLargestPossibleRegion().GetSize();
  start = im->GetLargestPossibleRegion().GetIndex();


  roiStart[0] = 0; roiStart[1] = 0; roiStart[2] = 0;
  roiEnd[0] = 0; roiEnd[1] = 0; roiEnd[2] = 0;

  unsigned int ndims = im->GetImageDimension();

  bool foundLabel = false;
  itk::ImageRegionIteratorWithIndex< itkImageType > label(im, im->GetBufferedRegion() );
  for (label.GoToBegin(); !label.IsAtEnd(); ++label)
  {
    if (label.Get() != 0)
    {
      typename itkImageType::IndexType idx = label.GetIndex();
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

  int radius = 17;
  for (unsigned i = 0; i < ndims; i++)
  {
    int diff = static_cast< int > (roiStart[i] - radius);
    if (diff >= start[i])
    {
      roiStart[i] -= radius;
    }
    else
    {
      roiStart[i] = start[i];
    }
    roiEnd[i] = (static_cast<unsigned int>(roiEnd[i] + radius) < size[i]) ? (roiEnd[i] + radius) : size[i] - 1;

  }

  // copy ROI to vector
  imROI.resize(6);
  for (unsigned i = 0; i < 3; i++)
  {
    imROI[i] = roiStart[i];
    imROI[i + 3] = roiEnd[i];
  }
}

template<typename PixelType>
void ExtractITKImageROI(typename itk::Image<PixelType, 3>::Pointer  im, const std::vector<long> &imROI, \
                        std::vector<PixelType> &imROIVec)
{

  // Copy itk image ROI to vector
  using ImageType = itk::Image<PixelType, 3>;
  typename ImageType::IndexType index;
  long i, j, k, kk, DIMXYZ;

  DIMXYZ = (imROI[3] - imROI[0]) * (imROI[4] - imROI[1]) * (imROI[5] - imROI[2]);
  imROIVec.clear();
  imROIVec.resize(DIMXYZ);
  kk = 0;
  for (k = imROI[2]; k < imROI[5]; k++)
    for (j = imROI[1]; j < imROI[4]; j++)
      for (i = imROI[0]; i < imROI[3]; i++)
      {
        index[0] = i; index[1] = j; index[2] = k;
        imROIVec[kk++] = im->GetPixel(index);
      }
}

template<typename PixelType>
void UpdateITKImageROI(const std::vector<PixelType> &imROIVec, const std::vector<long> &imROI,  \
                       typename itk::Image<PixelType, 3>::Pointer im)
{

  using ImageType = itk::Image<PixelType, 3>;
  typename ImageType::IndexType index;
  long i, j, k, kk;

  // Set non-ROI as zeros
  im->FillBuffer(0);
  kk = 0;
  for (k = imROI[2]; k < imROI[5]; k++)
    for (j = imROI[1]; j < imROI[4]; j++)
      for (i = imROI[0]; i < imROI[3]; i++)
      {
        index[0] = i; index[1] = j; index[2] = k;
        im->SetPixel(index, imROIVec[kk++]);
      }
}




/**
 * readImage
 */
template< typename TImageType >
typename TImageType::Pointer readImage(const char *fileName)
{
  using ImageReaderType = itk::ImageFileReader< TImageType >;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(fileName);

  typename TImageType::Pointer image;

  try
  {
    reader->Update();
    image = reader->GetOutput();
  }
  catch ( itk::ExceptionObject &err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    raise(SIGABRT);
  }

  return image;
}

/**
 * writeImage
 */
template< typename TImageType > void writeImage(typename TImageType::Pointer img, const char *fileName)
{
  using WriterType = itk::ImageFileWriter< TImageType >;

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( fileName );
  writer->SetInput(img);
  writer->UseCompressionOn();

  try
  {
    writer->Update();
  }
  catch ( itk::ExceptionObject &err)
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    raise(SIGABRT);
  }
}

template<typename VType>
void WriteVectorIntoFile(const char *fnSave, const std::vector<VType> &vec)
{

  try
  {
    std::ofstream outfile(fnSave);
    outfile.write((const char *)&vec[0], vec.size()*sizeof(VType));
  }
  catch ( itk::ExceptionObject &err)
  {
    std::cout << "Fail to write to file " << fnSave << std::endl;
    std::cout << err << std::endl;
    raise(SIGABRT);
  }
}

template<typename VType>
void LoadVectorIntoFile(const char *fnLoad, std::vector<VType> &vec, const long VEC_SIZE)
{

  try
  {
    std::ifstream infile(fnLoad);
    vec.resize(VEC_SIZE);
    infile.read((char *)&vec[0], vec.size()*sizeof(VType));
  }
  catch ( itk::ExceptionObject &err)
  {
    std::cout << "Fail to load file " << fnLoad << std::endl;
    std::cout << err << std::endl;
    raise(SIGABRT);
  }
}
}// douher

MaskImageType::Pointer GetMaskImage(LabelImageType::Pointer);
MaskImageType::Pointer GetMaskImage(LabelImageType::Pointer, LabelPixelType);

MaskImageType::RegionType RoiIndexToRegion(MaskImageType::IndexType, MaskImageType::IndexType);
void RegionToIndex(MaskImageType::RegionType, MaskImageType::IndexType &, MaskImageType::IndexType &);

MaskImageType::RegionType GetRoi(MaskImageType::Pointer);
void ExpandRoi(MaskImageType::Pointer, MaskImageType::RegionType &, const double objectSize = 10);
void BoundingCheck(MaskImageType::Pointer, InputImageType::Pointer, MaskImageType::RegionType &, InputImageType::RegionType &);
void BoundingCheck(MaskImageType::Pointer, LabelImageType::Pointer, MaskImageType::RegionType &, InputImageType::RegionType &);

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
