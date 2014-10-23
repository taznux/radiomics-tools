
clear; close all; clc; 

data_path = '/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/';
bin_path = '/home/wlui/data/work/wchoi/FeatureExtraction/BIN/';

SUVfixedImageFile=[data_path '1714542OUTPUT/PreSUV'];
SUVMovingOutput=[data_path '1714542OUTPUT/SUVRegistered'];
SUVDiffImage=[data_path '1714542OUTPUT/SUVDiffImage'];

MaskBySegmentationSUVfixed=[data_path '1714542OUTPUT/MaskBySegmentationSUVfixed'];
MaskBySegmentationSUVRegistered=[data_path '1714542OUTPUT/MaskBySegmentationSUVRegistered'];

Scale=100.00;

for InputImageCell = {SUVfixedImageFile, SUVMovingOutput, SUVDiffImage}
    InputImage = InputImageCell{1};
    MaskImage=[data_path '1714542OUTPUT/FixedWindowMaskImage'];

    systemCallString = sprintf([bin_path 'FeatureExtraction %s %d %s'], InputImage, Scale, MaskImage);   
    % systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d', InputImage, Scale);   


    tic;
    %system(systemCallString);
    system(['ssh node1 "cd ' pwd ';' systemCallString '"']);
    toc;
end
