#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkImageFileWriter.h"
#include "itkTextOutput.h"

#include "itksys/SystemTools.hxx"

#include <boost/algorithm/string.hpp>


using ShortPixelType = signed short;
using UShortPixelType = unsigned short;
using FloatPixelType = float;
using ShortImageType = itk::Image< ShortPixelType, 3 >;
using ShortDICOMImageType = itk::Image< ShortPixelType, 2 >;
using UShortImageType = itk::Image< UShortPixelType, 3 >;
using UShortDICOMImageType = itk::Image< UShortPixelType, 2 >;
using FloatImageType = itk::Image< FloatPixelType, 3 >;
using FloatDICOMImageType = itk::Image< FloatPixelType, 2 >;

using NamesGeneratorType = itk::GDCMSeriesFileNames;
using FileNamesContainer = std::vector< std::string >;

int refineString(std::string &in)
{
	boost::trim(in);
	boost::replace_all(in, "/", "_");
	boost::replace_all(in, " ", "_");

	return in.length();
}

std::string refineName(std::string &fullName)
{
	std::vector<std::string> strs;
	std::string refineName;

	boost::algorithm::to_lower(fullName);
	boost::split(strs, fullName, boost::is_any_of("^"));

	boost::trim_if(strs[0], !boost::algorithm::is_alpha());
	boost::trim_if(strs[1], !boost::algorithm::is_alpha());

	strs[0][0] = toupper(strs[0][0]);
	strs[1][0] = toupper(strs[1][0]);

	refineName = strs[1] + "_" + strs[0];

	return refineName;
}

std::string initialName(std::string &fullName)
{
	std::vector<std::string> strs;
	std::string initialName;

	boost::algorithm::to_lower(fullName);
	boost::split(strs, fullName, boost::is_any_of("^"));

	boost::trim_if(strs[0], !boost::algorithm::is_alpha());
	boost::trim_if(strs[1], !boost::algorithm::is_alpha());

	strs[0][0] = toupper(strs[0][0]);
	strs[1][0] = toupper(strs[1][0]);

	initialName = strs[1][0] + strs[0];

	return initialName;
}

std::string getModality(std::string fileName)
{
	using ReaderType = itk::ImageFileReader<ShortDICOMImageType>;
	using ImageIOType = itk::GDCMImageIO;

	ImageIOType::Pointer dicomIO = ImageIOType::New();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(dicomIO);

	reader->SetFileName(fileName);
	reader->Update();

	std::string modality;

	dicomIO->GetValueFromTag("0008|0060", modality);
	//std::cout << modality << std::endl;

	return modality;
}

template <typename TImageType>
int readAndWrite(FileNamesContainer fileNames, std::string seriesDirectory)
{
	using ReaderType = itk::ImageSeriesReader<TImageType>;
	using ImageIOType = itk::GDCMImageIO;

	ImageIOType::Pointer dicomIO = ImageIOType::New();
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(dicomIO);


	reader->SetFileNames(fileNames);

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
	std::string modality;
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

	dicomIO->GetValueFromTag("0008|0060", modality);
	std::cout << modality << std::endl;

	dicomIO->GetValueFromTag("0008|1030", studyDescription);
	std::cout << studyDescription << std::endl;

	dicomIO->GetValueFromTag("0008|103e", seriesDescription);
	std::cout << seriesDescription << std::endl;

	dicomIO->GetValueFromTag("0012|0050", timepointID);
	std::cout << timepointID << std::endl;

	dicomIO->GetValueFromTag("0018|1210", convolutionKernel);
	std::cout << convolutionKernel << std::endl;

	refineString(patientID);
	refineString(studyDate);
	refineString(aqusitionTime);
	refineString(studyDescription);
	refineString(seriesDescription);
	refineString(convolutionKernel);


	//
	seriesDirectory.append("/");
	seriesDirectory.append(patientID);// + "_" +refineName(patientName));

	itksys::SystemTools::MakeDirectory(seriesDirectory.c_str());
	//
	seriesDirectory.append("/");

	std::string nrrdFilename = seriesDirectory;
	nrrdFilename.append(patientID + "_" + studyDate + "_" + modality);
	/*if (seriesDescription.length() > 0)
	{
		nrrdFilename.append(seriesDescription);
	}
	else
	{
		if (convolutionKernel.length() > 0)
			nrrdFilename.append(convolutionKernel);
		else
			nrrdFilename.append(studyTime);
	}*/

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
	//using WriterType = itk::ImageSeriesWriter<ImageType>;

	//nameGenerator->SetOutputDirectory(seriesDirectory);
	//FileNamesContainer OutputImageNames;
	//OutputImageNames = nameGenerator->GetOutputFileNames();

	//WriterType::Pointer writer = WriterType::New();
	//writer->SetInput(reader->GetOutput());
	//writer->SetImageIO(dicomIO);
	//writer->SetFileNames(OutputImageNames);
	//writer->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray());

	std::cout << "Writing the image as " << std::endl << std::endl;
	std::cout << seriesDirectory << std::endl << std::endl;

	using SingleOutImageWriterType = itk::ImageFileWriter<TImageType>;
	typename SingleOutImageWriterType::Pointer OutWriter = SingleOutImageWriterType::New();
	OutWriter->SetUseCompression(true);
	OutWriter->SetInput(reader->GetOutput());
	OutWriter->SetFileName(nrrdFilename);
	OutWriter->Update();

	return 0;
}

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory OutputDicomDirectory"
			<< std::endl;
		return EXIT_FAILURE;
	}

	itk::OutputWindow::SetInstance(itk::TextOutput::New());

	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetUseSeriesDetails(true);
	nameGenerator->AddSeriesRestriction("0008|0021");
	nameGenerator->SetRecursive(true);
	nameGenerator->SetDirectory(argv[1]);

	std::string seriesDirectory = argv[2];

	try
	{
		using SeriesIdContainer = std::vector< std::string >;
		const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

		if (seriesUID.size() > 0)
		{
			std::cout << std::endl << "The directory: " << std::endl;
			std::cout << std::endl << argv[1] << std::endl << std::endl;
			std::cout << "Contains the following DICOM Series: ";
			std::cout << std::endl << std::endl;

			SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
			SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

			while (seriesItr != seriesEnd)
			{
				std::cout << seriesItr->c_str() << std::endl;
				++seriesItr;
			}

			std::string studyDescriptionOld;

			for (seriesItr = seriesUID.begin(); seriesItr != seriesEnd; seriesItr++)
			{
				std::string seriesIdentifier;
				seriesIdentifier = seriesItr->c_str();

				FileNamesContainer fileNames;
				fileNames = nameGenerator->GetFileNames(seriesIdentifier);

				std::cout << std::endl << std::endl;
				std::cout << "Now reading series: " << std::endl << std::endl;
				std::cout << seriesIdentifier << std::endl;
				std::cout << std::endl << std::endl;


				if (getModality(fileNames[0]) == "PT") // PET
				{
					readAndWrite<FloatImageType>(fileNames, seriesDirectory);
				}
				else // CT
				{
					readAndWrite<ShortImageType>(fileNames, seriesDirectory);
				}
			}
		}
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return EXIT_FAILURE;
	}
	

	return EXIT_SUCCESS;
}

