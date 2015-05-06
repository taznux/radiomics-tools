 
function MainForFeatureExtraction 
    
  
%    %% For 1714542
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/PreCT/58450157/00205658';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/PostCT/21575436/37294405';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
    folder= '/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT';
    BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);
    
%        
   
      % For 1701386
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/PreCT/04281545/27150748';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/PostCT/30483679/38092153';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder= '/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
   
      
      %% For 1517163
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/PreCT/download20110116162107/93017251/99462603';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/PostCT/download20110116162154/10852211/64913886';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
         %% For 1536392
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/PreCT/download20110114153157/43862064/72501036';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/PostCT/download20110114153244/73099432/39093624';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
      
         %% For 1563576 (NRROI doesent work for this data)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/PreCT/download20110116163510/98605267/08994178';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/PostCT/download20110116163617/10253792/59419185';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

%    
%       
         %% For 1573147
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/PreCT/download20110114153401/92243787/71082175';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/PostCT/download20110114153420/61796401/27234114';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
  folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

%    
            %% For 1574421 (NCROI: -0.931471)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/PreCT/download20110116162806/88046930/92345711';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/PostCT/download20110116162836/15651446/23710765';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);
% 
               %% For 1574991 (NCROI doesnot work for this data)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/PreCT/download20110114153622/33978155/64459386';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/PostCT/download20110114153657/43738900/83330266';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
               %% For 1637365 (NCROI:-0.8304 )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/PreCT/download20110113152338/39397181/71475301';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/PostCT/55965427/28540599';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

      
               %% For 1581225 (NCROI doesnot work for this data -0.528567 )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/PreCT/download20110115123347/98954825/64900683';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/PostCT/download20110115123445/16322672/79038275';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

%       
%                %% For   1652879 (NCROI -0.855128)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/PreCT/download20110116163201/94467804/46345768';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/PostCT/download20110116163232/95733551/07352674';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

% 
               %% For 1693445 (-0.871794)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/PreCT/47497438/77574685';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/PostCT/download20110113143129/31845707/80308312';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
               %% For 1712680 ()
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/PreCT/49185754/03390782';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/PostCT/download20110113140145/88183906/15072673';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

               %% For   1244212 (-0.846467)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/PreCT/74908134/01650767';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/PostCT/download20110113153036/73737401/87365162';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
                  %% For   1155223 ( -0.874222)
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/PreCT/download20110116162937/97739473/86121161';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/PostCT/download20110116163004/57322077/78001466';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

  
   % New five data
   
                     %% For   1250464 ( )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/PreCT/20919678';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/PostCT/download20110509124520/06862234/56206747';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

  
   
                       %% For   1612690 ( )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/PreCT/21032718CT';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/PostCT/12496033';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
      
                       %% For   1628916 ( )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/PreCT/download20110509115729/86267185/97513378';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/PostCT/download20110509120104/75112652/06412995';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
    MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
%                           %% For   1663560 ( )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/PreCT/01427923';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/PostCT/download20110509121455/14689916/65058614';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
  MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    % MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/MaskMaual3DSlicer';

    
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT';
   BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
                          %% For   1602486 ( )
    clear; close all; clc; 
    % Existing Folders
    CTfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/PreCT/download20110509113630/77301598/65870096';
    CTMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/PostCT/27701415';
    
    SUVfixedImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/PreSUV';
    SUVMovingImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/PostSUV';
      
    % ROI for registration   
    CTFixedROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTFixedROIIndexFile/FixedCTROIIndex.txt';  %% The text file including ROI index
    CTMovingROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTMovingROIIndexFile/MovingCTROIIndex.txt';
      
    % ROI for Tumore            
    SUVRoughTumoreROIIndexFile='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/SUVRoughTumoreROIIndexFile/SUVRoughTumoreROIIndex.txt';
    %%%    xStart=70;yStart=70;zStart=130;xExtent=15;yExtent=15;zExtent=16;
      
    % Folders to produce related to Registration
    CTFixedToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTFixedToFixedSUV';
    CTRegisteredToFixedSUV='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTRegisteredToFixedSUV';
    CTDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTDiffImage';

    SUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/SUVRegistered';
    SUVDiffImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/SUVDiffImage';
    CTfixedImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTfixedImageCroped';
    CTmovingImageCroped='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/CTmovingImageCroped';
    
   %Folders to Produce related to Mask Producing
  MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
      %  MaskBySegmentationSUVfixed='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/MaunalMask3DSlicer';

    
    MaskBySegmentationSUVRegistered='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT/MaskBySegmentationSUVRegistered';
    FixedWindowMaskImage='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/Temp';
    
   folder='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT';
   BatchProcessForFeatureExtractionFor1602486(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage);

   
   
   
   i
  
   % --------------------------------------------------------------------
function BatchProcessForFeatureExtraction(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage)
   
   
   if ~(exist( CTFixedToFixedSUV, 'dir') )        mkdir(CTFixedToFixedSUV);     end
   if ~(exist( CTRegisteredToFixedSUV, 'dir') )        mkdir(CTRegisteredToFixedSUV);     end
   if ~(exist( CTDiffImage, 'dir') )        mkdir(CTDiffImage);     end

   if ~(exist( SUVRegistered, 'dir') )        mkdir(SUVRegistered);     end
   if ~(exist( SUVDiffImage, 'dir') )        mkdir(SUVDiffImage);     end
   if ~(exist( CTfixedImageCroped, 'dir') )        mkdir(CTfixedImageCroped);     end
   if ~(exist( CTmovingImageCroped, 'dir') )        mkdir(CTmovingImageCroped);     end
       
   if ~(exist( MaskBySegmentationSUVfixed, 'dir') )        mkdir(MaskBySegmentationSUVfixed);     end
   if ~(exist( MaskBySegmentationSUVRegistered, 'dir') )        mkdir(MaskBySegmentationSUVRegistered);     end
   if ~(exist( FixedWindowMaskImage, 'dir') )        mkdir(FixedWindowMaskImage);     end


    SelectionOfSegmentationMethods=2; 
       
    %  2. Segmentation method with threshold 2.5

    systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/MaskProducing/BIN/MaskProducing %s %s %d %s %s %s %s',SUVfixedImage,SUVRegistered,SelectionOfSegmentationMethods,SUVRoughTumoreROIIndexFile, FixedWindowMaskImage, MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered);   
    tic;
    system(systemCallString);
    toc;
%   

%% For Feature computation
  
   olderfolder=cd(folder);
   Scale=100.00;
    % features based on Masks from Threshold = 2.5 
  
   systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d %s', SUVRegistered, Scale, MaskBySegmentationSUVRegistered);   
%    systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d %s', SUVRegistered, Scale, MaskBySegmentationSUVfixed);   

   system(systemCallString);
   
  cd(olderfolder);

  
  
  
  function BatchProcessForFeatureExtractionFor1602486(folder,CTfixedImage, CTMovingImage, SUVfixedImage,SUVMovingImage,CTFixedROIIndexFile,CTMovingROIIndexFile,SUVRoughTumoreROIIndexFile,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVRegistered,SUVDiffImage,CTfixedImageCroped,CTmovingImageCroped,MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered,FixedWindowMaskImage)
   
   
   if ~(exist( CTFixedToFixedSUV, 'dir') )        mkdir(CTFixedToFixedSUV);     end
   if ~(exist( CTRegisteredToFixedSUV, 'dir') )        mkdir(CTRegisteredToFixedSUV);     end
   if ~(exist( CTDiffImage, 'dir') )        mkdir(CTDiffImage);     end

   if ~(exist( SUVRegistered, 'dir') )        mkdir(SUVRegistered);     end
   if ~(exist( SUVDiffImage, 'dir') )        mkdir(SUVDiffImage);     end
   if ~(exist( CTfixedImageCroped, 'dir') )        mkdir(CTfixedImageCroped);     end
   if ~(exist( CTmovingImageCroped, 'dir') )        mkdir(CTmovingImageCroped);     end
       
   if ~(exist( MaskBySegmentationSUVfixed, 'dir') )        mkdir(MaskBySegmentationSUVfixed);     end
   if ~(exist( MaskBySegmentationSUVRegistered, 'dir') )        mkdir(MaskBySegmentationSUVRegistered);     end
   if ~(exist( FixedWindowMaskImage, 'dir') )        mkdir(FixedWindowMaskImage);     end

%    %%  For Registration
%    MetricOptions=2; %% 1: case 1:NormalizedCorrelation 2: MeanSquares 3: MattesMutualInformation; 4:NormalizedMutualInformation; 5. MutualInformation; 6: GradientDifference
% 		             %% Defult: NormalizedCorrelation
%    CTROIOptions=2; %% Selecting the ROI:  1: Whole image; 2: ROI; 3: ROI+Bone
%     
%     if CTROIOptions==1,
%     systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/PETCTMultiRigRegistraion/BIN/PETCTMultiRigRegistraion %s %s %s %s %s %s %s %s %s %d %d',CTfixedImage,CTMovingImage,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVfixedImage,SUVMovingImage, SUVRegistered, SUVDiffImage, MetricOptions, CTROIOptions);   
%     else    
%     systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/PETCTMultiRigRegistraion/BIN/PETCTMultiRigRegistraion %s %s %s %s %s %s %s %s %s %d %d %s %s %s %s',CTfixedImage,CTMovingImage,CTFixedToFixedSUV,CTRegisteredToFixedSUV,CTDiffImage,SUVfixedImage,SUVMovingImage, SUVRegistered, SUVDiffImage, MetricOptions, CTROIOptions, CTFixedROIIndexFile,  CTfixedImageCroped, CTMovingROIIndexFile, CTmovingImageCroped);   
%     end 
%       
%     tic;
%    system(systemCallString);
%     toc;
    
       
   %%  For Maskproducing
    SelectionOfSegmentationMethods=1;     
     %  1. Manually get a ROI peak region and extract a 3*3 region around
     %  the maximu SUV     
     %  2. Segmentation method with threshold 2.5

    systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/MaskProducingFor1602486/BIN/MaskProducingFor1602486 %s %s %d %s %s %s %s',SUVfixedImage,SUVRegistered,SelectionOfSegmentationMethods,SUVRoughTumoreROIIndexFile, FixedWindowMaskImage, MaskBySegmentationSUVfixed,MaskBySegmentationSUVRegistered);   
    tic;
    system(systemCallString);
    toc; 
% 
%%%  For patient 1602486, the ROI2.5 was gotten manually since the tumore is too close to the heart    
%     
%     

%% For Feature computation

  
   olderfolder=cd(folder);

   Scale=100.00;
  
   systemCallString = sprintf('/lui_tan/CodeAndData/ClearCode/FeatureExtraction/BIN/FeatureExtraction %s %d %s', SUVRegistered, Scale, MaskBySegmentationSUVRegistered);   
   system(systemCallString);
   
 
  cd(olderfolder);

  
