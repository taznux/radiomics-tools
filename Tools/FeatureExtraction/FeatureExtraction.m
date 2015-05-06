   
    clear; close all; clc; 

    
    SUVfixedImageFile='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/PreSUV';
    SUVMovingOutput='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/SUVRegistered';
    SUVDiffImage='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/SUVDiffImage';
       
    
  
    MaskBySegmentationSUVfixed='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/MaskBySegmentationSUVfixed';
    MaskBySegmentationSUVRegistered='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/MaskBySegmentationSUVRegistered';

    Scale=100.00;
     
    
   InputImage=SUVfixedImageFile;
   MaskImage='/lui_tan/CodeAndData/Data/Patients/Esophagus/1714542OUTPUT/FixedWindowMaskImage';
  
  systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d %s', InputImage, Scale, MaskImage);   
% systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d', InputImage, Scale);   

   
   tic;
   system(systemCallString);
   toc;
%       

%