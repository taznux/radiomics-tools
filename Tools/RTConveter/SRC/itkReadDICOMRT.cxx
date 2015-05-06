/*=========================================================================
  Program: DICOM-RT Structure file extractor
  Module: itkDICOMRT.cxx
  Author: Jason Dowling  
  Contact: Jason.Dowling@csiro.au - please email with any bugs or suggestions :)
  Modified by:
  Language: C++
  This ITK4+ version created Feb 2013  Original ITKv3 with local gdcm Sept 2009 ( http://www.insight-journal.org/browse/publication/701 )
  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
  Mathieu Malaterre: http://mathieumalaterre.com/ 

  Requires:  
  1. Recent ITK (>= v4.1.0 ) 
  
  Test data provided is a pelvic MRI with manual contours of bladder, prostate, rectum, etc.
 
  Update History:
  March 2013 Jason Dowling 
	1. mergeImages() now correctly handles holes within structures.  Thanks to Sebastian Steger for the code fix.
 
=========================================================================*/

#include <string>
#include <itkImage.h>
#include <itkPoint.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageIteratorWithIndex.h>
#include <itkPointSetToImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkPolygonSpatialObject.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "itkGDCMImageIO.h" 
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include <gdcmTypes.h>  
#include <gdcmReader.h>
#include <gdcmSmartPointer.h>
#include <gdcmAttribute.h>
#include <gdcmSequenceOfItems.h>
#include <gdcmFileMetaInformation.h>


const unsigned int Dimension = 3;
typedef int PixelType;
typedef itk::Point< double, Dimension > PointType;

typedef itk::PolygonSpatialObject<2> PolygonType;
typedef PolygonType::Pointer PolygonPointer;
typedef itk::SpatialObjectPoint<2> PolygonPointType;

typedef itk::Image< PixelType, Dimension >   ImageType;
typedef itk::Image< PixelType, 2 >   ImageSliceType;

typedef itk::GroupSpatialObject<2> GroupType;
typedef itk::SpatialObjectToImageFilter<GroupType, ImageSliceType> SpatialObjectToImageFilterType;

void FillMedicalImageInformation(const gdcm::Reader &reader);

ImageType::Pointer readFile(std::string fileName){

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName);
  try {
    reader->Update();
  } catch (itk::ExceptionObject & err) {
    std::cerr << "ExceptionObject caught:" << std::endl;
    std::cerr << err << std::endl;
    exit(1);
  }
  ImageType::Pointer image = reader->GetOutput();
  return image;
}

void writeFile(ImageType::Pointer image, std::string outputFileName){
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetUseCompression(true);
  writer->SetFileName(outputFileName);
  try {
    writer->Update();
  } catch (itk::ExceptionObject & err) {
    std::cerr << "ExceptionObject caught:" << std::endl;
    std::cerr << err << std::endl;
    exit(1);
  }
  std::cout << "File written: " << outputFileName <<  std::endl;
}

void write2DFile(ImageSliceType::Pointer image, std::string outputFileName){
  typedef itk::ImageFileWriter< ImageSliceType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName(outputFileName);
  try {
    writer->Update();
  } catch (itk::ExceptionObject & err) {
    std::cerr << "ExceptionObject caught:" << std::endl;
    std::cerr << err << std::endl;
    exit(1);
  }
  std::cout << "File written: " << outputFileName <<  std::endl;
}

void resetImage(ImageType::Pointer image)
{
  itk::ImageRegionIterator<ImageType> imageIt(image, image->GetLargestPossibleRegion()) ;
  imageIt.GoToBegin();
  while(!imageIt.IsAtEnd())  
  {
    imageIt.Set(0);
    ++imageIt;
  }
}
void reset2DImage(ImageSliceType::Pointer imageSlice)
{
  itk::ImageRegionIterator<ImageSliceType> imageIt(imageSlice, imageSlice->GetLargestPossibleRegion()) ;
  imageIt.GoToBegin();
  while(!imageIt.IsAtEnd())  
  {
    imageIt.Set(0);
    ++imageIt;
  }
}

void mergeImages(ImageSliceType::Pointer tempSlice, ImageType::Pointer finalImage, int iRequiredSlice )
{
  PixelType pixelValue =0;
  ImageType::IndexType pixelIndex;
  ImageSliceType::IndexType sliceIndex;
  int iX = finalImage->GetLargestPossibleRegion().GetSize()[0]; 
  int iY = finalImage->GetLargestPossibleRegion().GetSize()[1];

  
  if (iRequiredSlice>0)
  {  
  for (int i=0;i<iX;i++)
    for (int j=0;j<iY;j++)
    {
      pixelIndex[0] = i;
      pixelIndex[1] = j;
      sliceIndex[0] =i;
      sliceIndex[1] = j;

      pixelValue = tempSlice->GetPixel(sliceIndex);
      pixelIndex[2] = iRequiredSlice;    
      
       //Disable hole filling (if required please uncomment the next line (and comment the following line)).  
       //if (pixelValue != 0)  finalImage->SetPixel(pixelIndex, pixelValue  );
      finalImage->SetPixel(pixelIndex, finalImage->GetPixel(pixelIndex) ^ (pixelValue != 0));
    
    }
  }
}

//remove empty spaces from contour names
void trim (std::string &str)
{
    std::string temp;
    for (unsigned int i = 0; i < str.length(); i++)
        if (str[i] != ' ') temp += str[i];
    str = temp;
}

namespace gdcm { class Reader; }

PointType point;
//ITK spatial object to image filter segmentation definitions
PolygonType::PointListType pointList ;   typedef itk::ImageFileWriter< ImageType > WriterType;
 
PolygonType::Pointer polygon; 
PolygonPointType p;
SpatialObjectToImageFilterType::Pointer imageFilter =SpatialObjectToImageFilterType::New();
GroupType::Pointer group = GroupType::New();
std::string templateFilename;
float dcmOrigin =0;
ImageSliceType::Pointer temp2Dimage = ImageSliceType::New();

void insertRegion(GroupType * group, ImageType::Pointer finalImage, int iSlice )
{
  std::cout << "Inserting region with " << pointList.size() << " points into slice: " << iSlice << std::endl;  
  reset2DImage(temp2Dimage);
 
  //need to create a 2D slice here, put the polygon on it, and insert it back into the 3D volume...  
  group->AddSpatialObject(polygon); //add a new polygon group
 
  try
  {
     polygon->SetPoints(pointList);  //so copy them to a polygon object
     imageFilter->SetInput(group);
     imageFilter->SetSize(temp2Dimage->GetLargestPossibleRegion().GetSize());
     imageFilter->Update();
     temp2Dimage=imageFilter->GetOutput();
   }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "Problem setting polygon->SetPoints for this region (non-planar)" << std::endl;
    std::cerr << err << std::endl;

  }

  //merge new polygon from temp image into the contour image
  mergeImages(temp2Dimage, finalImage, iSlice);
   //remove the polygon and clean up pointlist
  group->RemoveSpatialObject(polygon);
   pointList.clear();
 }

int main(int argc, char * argv[])
{

  polygon=PolygonType::New();
  GroupType::Pointer group = GroupType::New();
  ImageType::IndexType pixelIndex;
  int iCurrentSlice = 0;
  int iPointsOutsideBoundary = 0;
  if( argc < 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  DicomDirectory inputDICOMRTImageFile OutputFilePrefix" << std::endl; 
    std::cerr << "eg. " << argv[0] << " E027 E027/RS.M1356941_.dcm Output_E027" << std::endl; 
    return 1;
  }

  //check the RS file is available before conversion
  const char *filename = argv[2];
  gdcm::Reader RTreader;
  RTreader.SetFileName( filename );
  if( !RTreader.Read() ) 
  {
    std::cout << "Problem reading file: " << filename << std::endl;
    return 0;
  }	

  //Step 1. is to read in the dicom image to retreive origin, spacing, etc.
  //the following code is based on the itk example: DicomSeriesReadImageWrite2.cxx
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  
  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );

  nameGenerator->SetDirectory( argv[1] );
  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << argv[1] << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;
    typedef std::vector< std::string >    SeriesIdContainer;
    
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }

    std::string seriesIdentifier;
      seriesIdentifier = seriesUID.begin()->c_str();


    std::cout << "Now reading series: " << seriesIdentifier << std::endl;

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
    reader->SetFileNames( fileNames );

     try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    std::string imageFilename = argv[3];
    imageFilename += "nrrd";
    
    writer->SetFileName( imageFilename);
    writer->SetInput( reader->GetOutput() );
    std::cout  << "Writing the image as " << imageFilename << std::endl;
    try
      {
      writer->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }
 
  //Step 2. Process the RS file  
  templateFilename = argv[2];
  ImageType::Pointer image = reader->GetOutput() ;
  resetImage( image); 
  dcmOrigin= image->GetOrigin()[2];
  ImageType::PointType origin;
  origin[0]=image->GetOrigin()[0];
  origin[1]=image->GetOrigin()[1];
  origin[2]=dcmOrigin;
  
  //we need to create a temporary 2D slice as well...
  ImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
  typedef itk::ExtractImageFilter< ImageType, ImageSliceType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  ImageType::SizeType size = inputRegion.GetSize();
  size[2] = 0;
  ImageType::IndexType start = inputRegion.GetIndex();
  start[2] =0;
  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );
  filter->SetDirectionCollapseToIdentity(); //22.02.2013
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( image );
  filter->Update();
  temp2Dimage = filter->GetOutput();

//  const gdcm::FileMetaInformation &h = RTreader.GetFile().GetHeader();
  const gdcm::DataSet& ds = RTreader.GetFile().GetDataSet();
  std::cout << "Parsing: " << filename << std::endl;
  
  gdcm::MediaStorage ms;
  ms.SetFromFile( RTreader.GetFile() );
  std::cout << "media storage: " << ms << std::endl;

  // (3006,0020) SQ (Sequence with explicit length #=4)      # 370, 1 StructureSetROISequence  
  gdcm::Tag tssroisq(0x3006,0x0020);
  if( !ds.FindDataElement( tssroisq ) )
  {
    std::cout << "Problem locating 0x3006,0x0020 - Is this a valid RT Struct file?" << std::endl;
    return 0;
  }
  gdcm::Tag troicsq(0x3006,0x0039);
  if( !ds.FindDataElement( troicsq ) )
  {
    std::cout << "Problem locating 0x3006,0x0039 - Is this a valid RT Struct file?" << std::endl;
    return 0;
  }

  const gdcm::DataElement &roicsq = ds.GetDataElement( troicsq );

  gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = roicsq.GetValueAsSQ();
  if( !sqi || !sqi->GetNumberOfItems() )
  {
    return 0;
  }
  const gdcm::DataElement &ssroisq = ds.GetDataElement( tssroisq );
  gdcm::SmartPointer<gdcm::SequenceOfItems> ssqi = ssroisq.GetValueAsSQ();
  if( !ssqi || !ssqi->GetNumberOfItems() )
  {
    return 0;
  }

  std::cout << "Number of structures found:" << sqi->GetNumberOfItems() << std::endl;

  //loop through structures
  for(unsigned int pd = 0; pd < sqi->GetNumberOfItems(); ++pd)
  {
    const gdcm::Item & item = sqi->GetItem(pd+1); // Item start at #1
    gdcm::Attribute<0x3006,0x0084> roinumber;
    const gdcm::DataSet& nestedds = item.GetNestedDataSet();
    roinumber.SetFromDataElement( nestedds.GetDataElement( roinumber.GetTag() ) );

    // find structure_set_roi_sequence corresponding to roi_contour_sequence (by comparing id numbers)
    unsigned int spd = 0;
    gdcm::Item & sitem = ssqi->GetItem(spd+1);
    gdcm::DataSet& snestedds = sitem.GetNestedDataSet();
    
    gdcm::Attribute<0x3006,0x0022> sroinumber;

    do
    {
      sitem = ssqi->GetItem(spd+1);
      snestedds = sitem.GetNestedDataSet();

      sroinumber.SetFromDataElement( snestedds.GetDataElement( sroinumber.GetTag() ) );
      
      spd++;

    } while ( sroinumber.GetValue()  != roinumber.GetValue() );

    gdcm::Tag stcsq(0x3006,0x0026);
    if( !snestedds.FindDataElement( stcsq ) )
    {
      std::cout<<"Did not find sttsq data el " << stcsq << "   continuing..." << std::endl;
      continue; //return 0;
    }
    const gdcm::DataElement &sde = snestedds.GetDataElement( stcsq );

    //(3006,002a) IS [255\192\96]                              # 10,3 ROI Display Color
    gdcm::Tag troidc(0x3006,0x002a);
    gdcm::Attribute<0x3006,0x002a> color = {};
    if( nestedds.FindDataElement( troidc) )
    {
      const gdcm::DataElement &decolor = nestedds.GetDataElement( troidc );
      color.SetFromDataElement( decolor );
    }
    //(3006,0040) SQ (Sequence with explicit length #=8)      # 4326, 1 ContourSequence
    gdcm::Tag tcsq(0x3006,0x0040);
    if( !nestedds.FindDataElement( tcsq ) )
    {
      continue; 
    }
    const gdcm::DataElement& csq = nestedds.GetDataElement( tcsq );

    gdcm::SmartPointer<gdcm::SequenceOfItems> sqi2 = csq.GetValueAsSQ();
    if( !sqi2 || !sqi2->GetNumberOfItems() )
    {
      std::cout << "csq: " << csq << std::endl;
      std::cout << "sqi2: " << *sqi2 << std::endl;
      std::cout<<"Did not find sqi2 or no. items == 0   " <<  sqi2->GetNumberOfItems() << "   continuing..." << std::endl;
      continue; 
    }
    unsigned int nitems = sqi2->GetNumberOfItems();
    std::cout << "Structure " << pd << ". Number of regions: " << nitems << std::endl;
    std::string str_currentOrgan(sde.GetByteValue()->GetPointer(), sde.GetByteValue()->GetLength());
    
    //trim to remove spaces in organ name which can cause problems in scripts eg. "CBCT01__BULK  BONE .nii" .  Might need to have this as parameter?
    trim (str_currentOrgan);
    std::cout << pd << ". Structure name: " << str_currentOrgan << std::endl;

    //now loop through each item for this structure (eg one prostate region on a single slice is an item)
    for(unsigned int i = 0; i < nitems; ++i)
    {
      const gdcm::Item & item2 = sqi2->GetItem(i+1); // Item start at #1

      const gdcm::DataSet& nestedds2 = item2.GetNestedDataSet();
      // (3006,0050) DS [43.57636\65.52504\-10.0\46.043102\62.564945\-10.0\49.126537\60.714... # 398,48 ContourData
      gdcm::Tag tcontourdata(0x3006,0x0050);
      const gdcm::DataElement & contourdata = nestedds2.GetDataElement( tcontourdata );

      //const gdcm::ByteValue *bv = contourdata.GetByteValue();
      gdcm::Attribute<0x3006,0x0050> at;
      at.SetFromDataElement( contourdata );
      const double* pts = at.GetValues();
      unsigned int npts = at.GetNumberOfValues() / 3;

      for(unsigned int j = 0; j < npts * 3; j+=3)
      {
        point[0] = pts[j+0];
        point[1] = pts[j+1];
        point[2] = pts[j+2];

        //transform points to image co-ordinates
        if (!(image->TransformPhysicalPointToIndex( point, pixelIndex )))
        {
	  //Are there points outside the image boundary.  This may occur with automatically segmented objects such as benches or external body outlines?  
	  iPointsOutsideBoundary++;
        }
	
	p.SetPosition(pixelIndex[0] ,pixelIndex[1],pixelIndex[2]);

        p.SetRed(1);
        p.SetBlue(1);
        p.SetGreen(1);
        pointList.push_back(p);
      }

      // we have the points for a contour in a single slice.  We need to join these up and insert into the slice as polygon.
      iCurrentSlice = pixelIndex[2];
      
      insertRegion(group, image, iCurrentSlice);
    }
  
    if (iPointsOutsideBoundary > 0) 
        {
	    std::cout <<  " --" << iPointsOutsideBoundary << " contour points detected outside image boundary. Please check the output volume. " ;
	    iPointsOutsideBoundary=0;
	}
  
    std::string strNewVolume;
    strNewVolume= argv[3];
    strNewVolume += "_";
    strNewVolume += str_currentOrgan;
    strNewVolume += ".nrrd";

    writeFile(image, strNewVolume);
    resetImage( image);  //reset the temporary volume ready for next structure (if any)

  } //next structure name

  return 0;
}
