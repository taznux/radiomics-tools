#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDenseFrequencyContainer2.h>
#include "itkRegionOfInterestImageFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkOffset.h"
#include "itkImageBase.h"
#include "itkIndex.h"
#include <cmath>

//definitions of used types
using InternalImageType = itk::Image<float, 3>;
using NeighborhoodType = itk::Neighborhood<float, 3>;
using Image2runLengthType = itk::Statistics::ScalarImageToRunLengthMatrixFilter<InternalImageType>;
using HistogramType = Image2runLengthType::HistogramType;
using Hist2FeaturesType = itk::Statistics::HistogramToRunLengthFeaturesFilter<HistogramType>;
using HistogramType = Image2runLengthType::HistogramType;
using OffsetType = InternalImageType::OffsetType;
using SpacingType = InternalImageType::SpacingType;

//calculate features for one offset
float * calcTextureFeatureImage (OffsetType offset,
    InternalImageType::Pointer inputImage, int numBins,
	float minVal, float maxVal, InternalImageType::Pointer maskImage)

{
    
	// Initialize an array of floats to store global features
	static float featuresV[10];

	// create an rlm generator
    Image2runLengthType::Pointer rlGenerator = Image2runLengthType::New();

	// define offset (hard-coded for now). should assign from the input offset.
	OffsetType list_offset;
	list_offset[0] = 1;
	list_offset[1] = 0;
	list_offset[2] = 0;

	// set offset
	rlGenerator->SetOffset(offset);
    rlGenerator->SetNumberOfBinsPerAxis(numBins); //reasonable number of bins
    rlGenerator->SetPixelValueMinMax(minVal, maxVal); //for input UCHAR pixel type
    
	/*
	InternalImageType::PointType pointMin;
    InternalImageType::PointType pointMax;
	itk::Index<3> minIndex;
	itk::Index<3> maxIndex;
	minIndex.SetElement(0,0);
	minIndex.SetElement(1,0);
	minIndex.SetElement(2,0);
	maxIndex.SetElement(0,64);
	maxIndex.SetElement(1,64);
	maxIndex.SetElement(2,64);	

    maskImage->TransformIndexToPhysicalPoint(minIndex, pointMin);
    maskImage->TransformIndexToPhysicalPoint(maxIndex, pointMax);

    double distMax = pointMin.EuclideanDistanceTo(pointMax);

	std::cout << "----- dist ------ ";
	std::cout << distMax << std::endl;

	rlGenerator -> SetDistanceValueMinMax (0,distMax);
	*/

	// create pointer to histogram object
	Hist2FeaturesType::Pointer featureCalc = Hist2FeaturesType::New();	

	// assign image and mask to the rlm generator
	rlGenerator -> SetInput(inputImage);
	rlGenerator -> SetMaskImage(maskImage);
	
	featureCalc -> SetInput (rlGenerator->GetOutput());
    featureCalc -> Update();

	// set correct max distance for this direction
	double numLevels = rlGenerator -> GetMax();	
	double voxelLength = 0;
	SpacingType spacing = inputImage -> GetSpacing();
	if (abs(offset[0]) == 1){
		voxelLength = spacing[0];
	}
	if (abs(offset[1]) == 1){
		voxelLength = pow(voxelLength,2) + pow(spacing[1],2);
		voxelLength = sqrt(voxelLength);
	}
	if (abs(offset[2]) == 1){
		voxelLength = pow(voxelLength,2) + pow(spacing[2],2);
		voxelLength = sqrt(voxelLength);
	}
	double dMax = voxelLength * numLevels;
	std::cout << "Maximum distance for: ";
	std::cout << offset << std::endl;
	std::cout << " is: ";
	std::cout << dMax << std::endl;
	rlGenerator -> SetDistanceValueMinMax (0,dMax);
	featureCalc -> Update();

	featuresV[0] = featureCalc->GetFeature(Hist2FeaturesType::ShortRunEmphasis); 
	featuresV[1] = featureCalc->GetFeature(Hist2FeaturesType::LongRunEmphasis); 
	featuresV[2] = featureCalc->GetFeature(Hist2FeaturesType::GreyLevelNonuniformity); 
	featuresV[3] = featureCalc->GetFeature(Hist2FeaturesType::RunLengthNonuniformity); 
	featuresV[4] = featureCalc->GetFeature(Hist2FeaturesType::LowGreyLevelRunEmphasis); 
	featuresV[5] = featureCalc->GetFeature(Hist2FeaturesType::HighGreyLevelRunEmphasis); 
	featuresV[6] = featureCalc->GetFeature(Hist2FeaturesType::ShortRunLowGreyLevelEmphasis); 
	featuresV[7] = featureCalc->GetFeature(Hist2FeaturesType::ShortRunHighGreyLevelEmphasis); 
	featuresV[8] = featureCalc->GetFeature(Hist2FeaturesType::LongRunLowGreyLevelEmphasis); 
	featuresV[9] = featureCalc->GetFeature(Hist2FeaturesType::LongRunHighGreyLevelEmphasis); 
    
	return featuresV;

}

int main(int argc, char*argv[])
{
  if(argc < 5)
    {
    std::cerr << "Usage: " << argv[0] << " Required image.mha" << std::endl;
    return EXIT_FAILURE;
    }
  
  std::string fileName = argv[1];
  std::string numBinsStr = argv[2];
  int numBins = atoi(numBinsStr.c_str());
  std::string minValStr = argv[3];
  float minVal = atof(minValStr.c_str());
  std::string maxValStr = argv[4];
  float maxVal = atof(maxValStr.c_str());
  std::string maskFileName = argv[5];

  // read image
  using ReaderType = itk::ImageFileReader<InternalImageType>;
  ReaderType::Pointer reader=ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();
  InternalImageType::Pointer image=reader->GetOutput();

  // read mask
  ReaderType::Pointer maskReader=ReaderType::New();
  maskReader->SetFileName(maskFileName);
  maskReader->Update();
  InternalImageType::Pointer mask=maskReader->GetOutput();

  // neighborhood to get offset vectors
  NeighborhoodType neighborhood;
  neighborhood.SetRadius(1);
  unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
  
  OffsetType offset;

  using WriterType = itk::ImageFileWriter<InternalImageType>;
  WriterType::Pointer writer=WriterType::New();

  for ( unsigned int d = 0; d < centerIndex; d++ )
  {
      offset = neighborhood.GetOffset(d);

	  float * featuresV;

	  featuresV = calcTextureFeatureImage(offset, image, numBins, minVal, maxVal, mask);

	  for ( int i = 0; i < 10; i++ )
	  {
		  std::cout << "*(featuresV + " << i << ") : ";
		  std::cout << *(featuresV + i) << std::endl;
      }
	
  }

  return EXIT_SUCCESS;
}
