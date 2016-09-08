Radiomics tools
===================
[![Build Status](https://travis-ci.org/taznux/radiomics-tools.svg?branch=master)](https://travis-ci.org/taznux/radiomics-tools)

Image processing tools and ruffus based pipeline for radiomics feature analysis

Install software
--------------------------
### Python 3.5 ###
  - Required
  - Script engine and useful modules

  http://conda.pydata.org/miniconda.html - !recommend!  
  https://www.python.org/downloads/

  > Required modules - pandas, ruffus, SimpleITK, scipy, numpy, ipython, matplotlib
  >
  > install_modules.sh or install_modules.cmd is available to install these modules using conda.

### Slicer 4.5 ###
  - Recommended  
  - Image viewer, contour editor, simple image processing tool  

  http://download.slicer.org/


### To build  ###
  - gcc or visual studio
  - cmake
  - ITK


Tools
-----
### 1. DICOMTools ###
1. DICOMTagReader - Display entire DICOM tags
  > DICOMTagReader [DICOM directory]

2. DICOM2NRRDConverter - DICOM to nrrd (Slicer file format)  
  Simple recursive converting for single patient data  
  > DICOM2NRRDConverter [DICOM directory] [nrrd directory]  

  For large data  
  > python DICOM2NRRDConverter.py [DICOM directory] [nrrd directory]

3. DICOM-RT2NRRDConverter - DICOM-RT to nrrd


### 2. ContourTools ###
1. STAPLEComparison - variation comparison on multiple contours
1. ExtractBoundary
1. GTVs2ITV
1. HoleGenerator
1. ROIGenerator
1. ROI2BinImage
1. ROICropImage



### 3. GrowCutSegmentation ###
  NoduleSegmentation - Segment small nodular objects for solid nodule and GGO
  > NoduleSegmentation [InputImageFile] [SeedPoint_x] [SeedPoint_y [SeedPoint_z] [NoduleSize_long] [NoduleSize_short] [OutputImageFile]  


### 4. Feature Extraction ###
  FeatureExtraction - Extract image features from the nodule segmentation
  > FeatureExtraction [InputImage] [LabelImage] [FeatureFile] [Label={1}]


### 5. Python Tools ###
  1. metadata.py - for handling metadata in csv or xls
  2. organize_features.py - for collecting feature data into a single csv file


### 6. MATLAB Tools ###
  1. NRRD4Matlab - for handing nrrd format in MATLAB
  2. PET2SUV - for converting raw PET image to standardized uptake value(SUV)


### 7. ETC ###
  1. RegistrationSITK - simple registration code, required SimpleITK module for python
  2. SlicerPythonExtensions - simple extensions for Slicer
    1. InterpolateROIsEffect.py
    2. LineProfile.py


### 6. LASSO-SVM ###
  TBD - modeling code for radiomics features


Usage
-----
Radiomics feature extraction pipeline example for LUNGx dataset

1. Download DICOM images  
  https://wiki.cancerimagingarchive.net/display/Public/SPIE-AAPM+Lung+CT+Challenge  

  Download all DICOM images to './DATA'  
  You can use the included metadata files for LUNGx (TrainingSet.csv and TestSet.csv)  

2. Environmental parameters  
  Set your parameters in run_CAD.py (recommend default setting).  
  > experiment_set = 'TrainingSet'  
  > \# experiment_set = 'TestSet'  
  > output_path = 'output'  
  > data_path = 'DATA'  
  > dicom_path = data_path + '/DOI'  
  > image_path = data_path + '/' + experiment_set  
  > nodule_info_path = './' + experiment_set + '.csv'  

3. Run radiomics pipeline
  > $ python run.py or ./run.py

4. Analysis feature data
  output files (intermediate images and feature data) will be generated in ./output  
  > TrainingSet              feature_list_TrainingSet.csv
  > TestSet                  feature_list_TestSet.csv
