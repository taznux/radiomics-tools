#include "ITKUtils.h"

#include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"

#include "itkImageDuplicator.h"

using namespace std;

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " InputLabelImageFile SeedPoint[x y z] HoleSize OutputPrefix" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int x = atoi( argv[2] );
  const unsigned int y = atoi( argv[3] );
  const int sliceNumber = atoi( argv[4] );
  const double holeD = atof( argv[5] );
  const double holeR = holeD / 2;
  const string outputPrefix = argv[6];


  /////////To Read Label images ////////////////////
  LabelImageType::Pointer labelImage = ReadImageFile<MaskImageType>(argv[1]);
  LabelImageType::Pointer sphereImage = LabelImageType::New();
  LabelImageType::Pointer holedImage;

  sphereImage->CopyInformation(labelImage);
  sphereImage->SetBufferedRegion(labelImage->GetBufferedRegion());
  sphereImage->Allocate();
  sphereImage->FillBuffer( 0 );

  typedef itk::ImageDuplicator< LabelImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(labelImage);
  duplicator->Update();
  holedImage = duplicator->GetOutput();;

  // set ROI
  LabelImageType::PointType origin = labelImage->GetOrigin();
  LabelImageType::IndexType index = labelImage->GetLargestPossibleRegion().GetIndex();
  LabelImageType::SizeType size = labelImage->GetLargestPossibleRegion().GetSize();
  LabelImageType::SpacingType spacing = labelImage->GetSpacing();
  LabelImageType::SizeType roiSize;
  LabelImageType::IndexType roiStart;
  LabelImageType::IndexType roiEnd;
  LabelImageType::IndexType center;

  roiSize[0] = static_cast<unsigned int>((holeD + spacing[2]*2) / spacing[0]);
  roiSize[1] = static_cast<unsigned int>((holeD + spacing[2]*2) / spacing[1]);
  roiSize[2] = static_cast<unsigned int>((holeD + spacing[2]*2) / spacing[2]);

  center[0] = x;
  center[1] = y;

  if (sliceNumber > 0)
    center[2] = size[2] - sliceNumber - 1;
  else
    center[2] = -sliceNumber;

  std::cout << center << std::endl;

  roiStart[0] = center[0] - roiSize[0] / 2;
  roiStart[1] = center[1] - roiSize[1] / 2;
  roiStart[2] = center[2] - roiSize[2] / 2;
  roiEnd[0] = roiStart[0] + roiSize[0];
  roiEnd[1] = roiStart[1] + roiSize[1];
  roiEnd[2] = roiStart[2] + roiSize[2];

  std::cout << roiSize << std::endl;
  std::cout << roiStart << std::endl;
  std::cout << roiEnd << std::endl;

  // make a sphere
  {
    typedef itk::EllipseSpatialObject<3> EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, LabelImageType> SpatialObjectToImageFilterType;

    typedef EllipseType::TransformType TransformType;

    SpatialObjectToImageFilterType::Pointer sphereFilter = SpatialObjectToImageFilterType::New();

    TransformType::Pointer transform = TransformType::New();

    transform->SetIdentity();

    TransformType::OutputVectorType  translation;
    translation[0] = roiSize[0] * spacing[0] / 2.0;
    translation[1] = roiSize[1] * spacing[1] / 2.0;
    translation[2] = roiSize[2] * spacing[2] / 2.0;
    transform->Translate( translation, false );

    std::cout << translation << std::endl;

    EllipseType::ArrayType radiusArray;
    EllipseType::Pointer sphere = EllipseType::New();
    radiusArray[0] = holeR;
    radiusArray[1] = holeR;
    radiusArray[2] = holeR;

    sphere->SetRadius(radiusArray);
    sphere->SetDefaultInsideValue(maskValue);
    sphere->SetDefaultOutsideValue(0);
    sphere->SetObjectToParentTransform(transform);

    sphereFilter->SetSize(roiSize);
    sphereFilter->SetSpacing(spacing);
    sphereFilter->SetInput(sphere);
    sphereFilter->SetUseObjectValue(true);
    sphereFilter->SetOutsideValue(0);

    try
    {
      sphereFilter->Update();

      itk::ImageRegionIteratorWithIndex<LabelImageType> sphere(sphereImage, sphereImage->GetBufferedRegion() );
      itk::ImageRegionIterator<LabelImageType> holed(holedImage, holedImage->GetBufferedRegion() );
      itk::ImageRegionConstIterator<LabelImageType> sphereROI(sphereFilter->GetOutput(), sphereFilter->GetOutput()->GetBufferedRegion());

      for (sphere.GoToBegin(), holed.GoToBegin(), sphereROI.GoToBegin(); !sphere.IsAtEnd();
           ++sphere,++holed)
      {
        LabelImageType::IndexType idx = sphere.GetIndex();

        if (idx[0] >= roiStart[0] && idx[1] >= roiStart[1] && idx[2] >= roiStart[2] &&
            idx[0] < roiEnd[0] && idx[1] < roiEnd[1] && idx[2] < roiEnd[2])
        {
          LabelImageType::PixelType color = sphereROI.Get();

          ++sphereROI;

          if (color > 0)
          {
            sphere.Set(maskValue); // sphere
            holed.Set(0); // hole
          }
        }
      }
    }
    catch ( itk::ExceptionObject &excp )
    {
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  }


  try
  {
    WriteImageFile<LabelImageType>(sphereImage, outputPrefix + "_sphere_" + argv[5] + "-label.nrrd");
    WriteImageFile<LabelImageType>(holedImage, outputPrefix + "_holed_" + argv[5] + "-label.nrrd");
  }
  catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << "Saved the label image" << std::endl;


  return EXIT_SUCCESS;
}
