#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkImageFileWriter.h"

#include "itksys/SystemTools.hxx"

#include <boost/algorithm/string.hpp>



int refineString(std::string &in)
{
    boost::trim(in);
    boost::replace_all(in, ",", "");
    boost::replace_all(in, " ", "_");

    return in.length();
}

std::string refineName(std::string &fullName)
{
    std::vector<std::string> strs;
    std::string refineName;

    boost::algorithm::to_lower(fullName);
    boost::split(strs, fullName, boost::is_any_of("^"));

    boost::trim_if(strs[0], ! boost::algorithm::is_alpha());
    boost::trim_if(strs[1], ! boost::algorithm::is_alpha());

    strs[0][0] = std::toupper(strs[0][0]);
    strs[1][0] = std::toupper(strs[1][0]);

    refineName = strs[1] + "_" + strs[0];

    return refineName;
}

std::string initialName(std::string &fullName)
{
    std::vector<std::string> strs;
    std::string initialName;

    boost::algorithm::to_lower(fullName);
    boost::split(strs, fullName, boost::is_any_of("^"));

    boost::trim_if(strs[0], ! boost::algorithm::is_alpha());
    boost::trim_if(strs[1], ! boost::algorithm::is_alpha());

    strs[0][0] = std::toupper(strs[0][0]);
    strs[1][0] = std::toupper(strs[1][0]);

    initialName = strs[1][0] + strs[0];

    return initialName;
}


int main( int argc, char *argv[] )
{
    if ( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " DicomDirectory OutputDicomDirectory"
                  << std::endl;
        return EXIT_FAILURE;
    }

    typedef signed short    PixelType;
    typedef itk::Image< PixelType, 3 >         ImageType;
    typedef itk::Image< PixelType, 2 >         DICOMImageType;

    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails( true );
    nameGenerator->AddSeriesRestriction("0008|0021" );
    nameGenerator->SetRecursive( true );
    nameGenerator->SetDirectory( argv[1] );

    try
    {
        std::cout << std::endl << "The directory: " << std::endl;
        std::cout << std::endl << argv[1] << std::endl << std::endl;
        std::cout << "Contains the following DICOM Series: ";
        std::cout << std::endl << std::endl;


        typedef std::vector< std::string >    SeriesIdContainer;
        const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();
        SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
        SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

        while ( seriesItr != seriesEnd )
        {
            std::cout << seriesItr->c_str() << std::endl;
            ++seriesItr;
        }

        std::string studyDescriptionOld;

        for ( seriesItr = seriesUID.begin(); seriesItr != seriesEnd ; seriesItr++)
        {
            std::string seriesIdentifier;
            seriesIdentifier = seriesItr->c_str();

            std::cout << std::endl << std::endl;
            std::cout << "Now reading series: " << std::endl << std::endl;
            std::cout << seriesIdentifier << std::endl;
            std::cout << std::endl << std::endl;

            typedef itk::ImageSeriesReader< ImageType >        ReaderType;
            ReaderType::Pointer reader = ReaderType::New();

            typedef itk::GDCMImageIO       ImageIOType;
            ImageIOType::Pointer dicomIO = ImageIOType::New();
            reader->SetImageIO( dicomIO );

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


            std::string patientID;
            std::string patientName;
            std::string studyDate;
            std::string studyTime;
            std::string seriesTime;
            std::string aqusitionTime;
            std::string studyDescription;
            std::string seriesDescription;
            std::string timepointID;
            std::string convolutionKernel;

            dicomIO->GetValueFromTag("0010|0010", patientName);
            std::cout << patientName << std::endl;

            dicomIO->GetValueFromTag("0010|0020", patientID);
            std::cout << patientID << std::endl;
            
            dicomIO->GetValueFromTag("0008|0020", studyDate);
            std::cout << studyDate << std::endl;

            dicomIO->GetValueFromTag("0008|0030", studyTime);
            std::cout << studyTime << std::endl;

            dicomIO->GetValueFromTag("0008|0031", seriesTime);
            std::cout << seriesTime << std::endl;

            dicomIO->GetValueFromTag("0008|0032", aqusitionTime);
            std::cout << aqusitionTime << std::endl;

            dicomIO->GetValueFromTag("0008|1030", studyDescription);
            std::cout << studyDescription << std::endl;

            dicomIO->GetValueFromTag("0008|103e", seriesDescription);
            std::cout << seriesDescription << std::endl;

            dicomIO->GetValueFromTag("0012|0050", timepointID);
            std::cout << seriesDescription << std::endl;

            dicomIO->GetValueFromTag("0018|1210", convolutionKernel);
            std::cout << seriesDescription << std::endl;

            refineString(patientID);
            refineString(studyDate);
            refineString(aqusitionTime);
            refineString(studyDescription);
            refineString(seriesDescription);
            refineString(convolutionKernel);

            std::string seriesDirectory = argv[2];
            //
            seriesDirectory.append("/");
            seriesDirectory.append(patientID);// + "_" +refineName(patientName));

            itksys::SystemTools::MakeDirectory(seriesDirectory.c_str());
            //
            seriesDirectory.append("/");

            std::string nrrdFilename = seriesDirectory;
            if(convolutionKernel.length() > 0)
                nrrdFilename.append(patientID + "_" + timepointID + "_" + studyDate + "_" + convolutionKernel);
            else
                nrrdFilename.append(patientID + "_" + timepointID + "_" + studyDate);

            nrrdFilename.append(".nrrd");

            /*
            if (studyDescription.length() > 0)
            {
                seriesDirectory.append(studyDescription);
                studyDescriptionOld = studyDescription;
            }
            else
            {
                seriesDirectory.append(studyDescriptionOld);
            }
            //
            seriesDirectory.append("/");
            std::string nrrdFilename = seriesDirectory;
            if (seriesDescription.length() > 0)
            {
                seriesDirectory.append(seriesDescription);
                nrrdFilename.append(initialName(patientName) + "_" + seriesDescription + ".nrrd");
            }
            else
            {
                seriesDirectory.append(seriesIdentifier);
                nrrdFilename.append(seriesIdentifier + ".nrrd");
            }
            */




            ////////To write Output images /////////////////////
            typedef itk::ImageSeriesWriter< ImageType, DICOMImageType > WriterType;

            nameGenerator->SetOutputDirectory( seriesDirectory );
            FileNamesContainer OutputImageNames;
            OutputImageNames = nameGenerator->GetOutputFileNames();

            WriterType::Pointer writer = WriterType::New();
            writer->SetInput(reader->GetOutput());
            writer->SetImageIO(dicomIO);
            writer->SetFileNames(OutputImageNames);
            writer->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray());

            std::cout  << "Writing the image as " << std::endl << std::endl;
            std::cout  << seriesDirectory << std::endl << std::endl;

            typedef itk::ImageFileWriter< ImageType> SingleOutImageWriterType;
            SingleOutImageWriterType::Pointer OutWriter = SingleOutImageWriterType::New();
            OutWriter->SetUseCompression(true);
            OutWriter->SetInput(reader->GetOutput());
            OutWriter->SetFileName(nrrdFilename);
            OutWriter->Update();

/*
            try
            {
                writer->Update();
            }
            catch (itk::ExceptionObject &ex)
            {
                std::cout << ex << std::endl;
                return EXIT_FAILURE;
            }
*/
        }
    }
    catch (itk::ExceptionObject &ex)
    {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
