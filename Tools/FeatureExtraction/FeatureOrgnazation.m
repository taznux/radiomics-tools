function FeatureOrgnation
   
clear; close all; clc;

FeatureName=[];
Feature=[];
Response=[];
PatientID=[];

%% For 1714542
Patient=1714542;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1714542/1714542OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];

Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

FeatureName=SubFeatureName;


%% For 1693445
Patient=1693445;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1693445/1693445OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

%% For 1712680
Patient=1712680;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1712680/1712680OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% The three above with Res=0
%% For 1701386

Patient=1701386;  PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1701386/1701386OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];

Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% This patient with Res=2

%% For 1517163

Patient=1517163;  PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1517163/1517163OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];

Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

%    % Res=1
%% For 1536392
Patient=1536392; PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1536392/1536392OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];


Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

%     % Res=2
%% For 1563576
Patient=1563576; PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1563576/1563576OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];     Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% RES=2
%% For 1573147
Patient=1573147; PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1573147/1573147OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];


Response=[Response; PatientRes]; PatientID=[PatientID; Patient];
% with RES=1
%% For 1574421
Patient=1574421; PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574421/1574421OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];


Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% with RES=1
%% For 1574991
Patient=1574991;PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1574991/1574991OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];  Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% RES=2
%% For 1637365
Patient=1637365; PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1637365/1637365OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];  Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% RES=2

%% For 1581225
Patient=1581225; PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1581225/1581225OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];
% RES=0

%% For   1652879
Patient=1652879; PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1652879/1652879OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];
% RES=2
%% For   1244212
Patient=1244212;PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1244212/1244212OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature];  Response=[Response; PatientRes]; PatientID=[PatientID; Patient];
% RES=2
%% For   1155223
Patient=1155223;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1155223/1155223OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

% RES=1


%%% New Patients %%%%

%% For   1250464 ( )
Patient=1250464;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1250464/1250464OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];


%% For   1612690 ( )

Patient=1612690;PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1612690/1612690OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];


%% For   1628916 ( )
Patient=1628916;PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1628916/1628916OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

%    
 %% For   1663560 ( )
   Patient=1663560;PatientRes=0;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1663560/1663560OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];

   
%% For   1602486 ( )
Patient=1602486;PatientRes=1;
[folder]='/home/wlui/data/data-from-lui_tan/Patients/PETCTDATA/1602486/1602486OUTPUT';
[SubFeature,SubFeatureName]=FeatureOrganzation2(folder);
Feature=[Feature SubFeature]; Response=[Response; PatientRes]; PatientID=[PatientID; Patient];



%% ROC analysis

Feature=Feature'; % Feature(find(isnan(Feature)))=0;
RocData=[]; AUCList=zeros(length(Feature),1);SerrorList=zeros(length(Feature),1);

%% Remove redundancy
[B,I,J]=unique(Feature','rows');

%   [B,I,J] = UNIQUE(...) also returns index vectors I and J such
%   that B = A(I) and A = B(J) (or B = A(I,:) and A = B(J,:)).
FeatureName=FeatureName(I);

B1=repmat(B(:,1),1,size(B,2));
B_T=B-B1;
Index=find(sum(abs(B_T),2));
% [Index1, Index2]=isequal(B_T,zeros(size(B,2),1));
FeatureName=FeatureName(Index);
B=B(Index,:);
Feature=B';
RocData = colAUC(Feature, Response);
FullResult=[RocData' Feature'];
[AUCRankFeature,index] =sortrows(FullResult,-1);
AUClist=AUCRankFeature(:,1);
AUCRandofFeatureName=FeatureName(index);


PureFeatureName=AUCRandofFeatureName; 

for k=1:length(AUClist)
       
    
    AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), num2str(AUClist(k)));
    % Test P-value
    x=[AUCRankFeature(k,2:end)' Response];
    [p,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
    
    %
    PValue=strcat('P=', num2str(p));
    AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), PValue);
end



FulllResultPlusID=zeros(size(AUCRankFeature,1)+2,size(AUCRankFeature,2)+1);
FulllResultPlusID(3:end,2:end)=AUCRankFeature;
FulllResultPlusID(1,3:end)=PatientID;
FulllResultPlusID(2,3:end)=Response;
FulllResultPlusID(3:end,1)=[1:length(FeatureName)];
FulllResultPlusID(1:2,1:2)=[NaN NaN; NaN NaN];


%% To get traditionoal features
close all; AUCofTraditionalFeatures={};TraditionalFeatures=[];TradionalFullResult=[];


index_SUV_max_Pre=find(strcmp('PreSUVMaskBySegmentationSUVfixedFeatures_Maximum=', PureFeatureName));
index_SUV_max_Post=find(strcmp('SUVRegisteredMaskBySegmentationSUVfixedFeatures_Maximum=', PureFeatureName));

index_Peakmean_Pre=find(strcmp('PreSUVFixedWindowMaskImageFeatures_Mean=', PureFeatureName));
index_Peakmean_Post=find(strcmp('SUVRegisteredFixedWindowMaskImageFeatures_Mean=', PureFeatureName));

SUV_max_Pre=AUCRankFeature(index_SUV_max_Pre,2:end);
SUV_max_Post=AUCRankFeature(index_SUV_max_Post,2:end);

SUV_Peakmean_pre=AUCRankFeature(index_Peakmean_Pre,2:end);
SUV_Peakmean_Post=AUCRankFeature(index_Peakmean_Post,2:end);

SUV_max_decrease=SUV_max_Pre-SUV_max_Post;
SUV_max_ratio=SUV_max_Post./SUV_max_Pre;

SUV_peak_decrease=SUV_Peakmean_pre-SUV_Peakmean_Post;
SUV_peak_ratio=SUV_Peakmean_Post./SUV_Peakmean_pre;



AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'PreSUVMaskBySegmentationSUVfixedFeatures_Maximum='];
TraditionalFeatures=[TraditionalFeatures;SUV_max_Pre];

AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUVRegisteredMaskBySegmentationSUVfixedFeatures_Maximum='];
TraditionalFeatures=[TraditionalFeatures;SUV_max_Post];

AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'PreSUVFixedWindowMaskImageFeatures_Mean='];
TraditionalFeatures=[TraditionalFeatures;SUV_Peakmean_pre];

AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUVRegisteredFixedWindowMaskImageFeatures_Mean='];
TraditionalFeatures=[TraditionalFeatures;SUV_Peakmean_Post];


AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUV_max_decrease='];
TraditionalFeatures=[TraditionalFeatures;SUV_max_decrease];

AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUV_max_ratio='];
TraditionalFeatures=[TraditionalFeatures;SUV_max_ratio];



AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUV_peak_decrease='];
TraditionalFeatures=[TraditionalFeatures;SUV_peak_decrease];

AUCofTraditionalFeatures=[AUCofTraditionalFeatures;'SUV_peak_ratio='];
TraditionalFeatures=[TraditionalFeatures;SUV_peak_ratio];



TradionalRocData = colAUC(TraditionalFeatures', Response);
[AUCRankofTradionalFeature,AUCTradionalindex] =sort(TradionalRocData');


AUCofTraditionalFeatures=AUCofTraditionalFeatures(AUCTradionalindex);
TraditionalFeatures=TraditionalFeatures(AUCTradionalindex,:);


for k=1:length(TradionalRocData)
       
    
    AUCofTraditionalFeatures(k)=strcat(AUCofTraditionalFeatures(k), num2str(AUCRankofTradionalFeature(k)));
    % Test P-value
    x=[TraditionalFeatures(k,:)' Response];
    [p,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
    
    %
    PValue=strcat('P=', num2str(p));
    AUCofTraditionalFeatures(k)=strcat(AUCofTraditionalFeatures(k), PValue);
end



%% to get three group:
IndexIntensityFeatures=[4;6;9;10;11;12;17];
IndexTexturalFeatures=[1;2;3;5;8;15];
IndexShapeFeatures=[7;13;14;16;18];

IntensityFeaturesName=AUCRandofFeatureName(IndexIntensityFeatures);IntensityFeatures=AUCRankFeature(IndexIntensityFeatures,2:end);
TexturalFeaturesName=AUCRandofFeatureName(IndexTexturalFeatures);TexturalFeatures=AUCRankFeature(IndexTexturalFeatures,2:end);
ShapeFeaturesName=AUCRandofFeatureName(IndexShapeFeatures);ShapeFeatures=AUCRankFeature(IndexShapeFeatures,2:end);

% Corr_IntensityFeatures=corr(IntensityFeatures');
% figure;imagesc(Corr_IntensityFeatures);title('Correlation of Intensity Features');
% 
% Corr_TexturalFeatures=corr(TexturalFeatures');
% figure;imagesc(Corr_TexturalFeatures);title('Correlation of Texture Features');
% 
% Corr_ShapeFeatures=corr(ShapeFeatures');
% figure;imagesc(Corr_ShapeFeatures);title('Correlation of Shape Features');


IntensityFeaturesPlusID=[FulllResultPlusID(1:2,3:end);IntensityFeatures];
[B,IX] = sort(IntensityFeaturesPlusID(2,:),'descend'); 
IntensityFeaturesPlusIDSorted=IntensityFeaturesPlusID(:,IX);


TexturalFeaturesPlusID=[FulllResultPlusID(1:2,3:end);TexturalFeatures];
[B,IX] = sort(TexturalFeaturesPlusID(2,:),'descend'); 
TexturalFeaturesPlusIDSorted=TexturalFeaturesPlusID(:,IX);

ShapeFeaturesPlusID=[FulllResultPlusID(1:2,3:end);ShapeFeatures];
[B,IX] = sort(ShapeFeaturesPlusID(2,:),'descend'); 
ShapeFeaturesPlusIDSorted=ShapeFeaturesPlusID(:,IX);

% Corr_features(find(Corr_features>0.7))=1;
% Corr_features(find(Corr_features<0.7))=0;
% figure;imagesc(Corr_features);




%% To compute the correlation between features with P value smaller than
%% 0.1


IndexWithSmallPValue=23; 
FeaturesWithSmallPValue=AUCRankFeature(1:IndexWithSmallPValue,2:end); 

% i=0;
% while i<length(FeaturesWithSmallPValue)
% for i=1:length(FeaturesWithSmallPValue)
%     for j=i+1:length(FeaturesWithSmallPValue)
%     corr_coef=corr(FeaturesWithSmallPValue(i,:)',FeaturesWithSmallPValue(j,:)');
%     if abs(corr_coef) > 0.7; 
%     FeaturesWithSmallPValue(j,:)=[]; 
%     break;
%     end
%     end
%     
%     if abs(corr_coef) > 0.7; 
%     break;
%     end
% end
% 
% end

Corr_features=corr(FeaturesWithSmallPValue');
figure;imagesc(Corr_features);

Corr_features(find(Corr_features>0.7))=1;
 Corr_features(find(Corr_features<0.7))=0;
 figure;imagesc(Corr_features);

figure;plot(Corr_features(1,:));



%% To compute correlation between features

Intensity_features_index=[2
6
13
19

24
37
38
44
11
];

Texture_features_index=[1
5
16
32
33
43
45
48
50
];

Shape_features_index=[4
12
14
15
18
21
22
25
31
40
46
47
49
55
];


Intensity_features=AUCRankFeature(Intensity_features_index,2:end); 
Corr_Intensity_features=corr(Intensity_features');
figure;imagesc(Corr_Intensity_features);

Texture_features=AUCRankFeature(Texture_features_index,2:end); 
Corr_Texture_features=corr(Texture_features');
figure;imagesc(Corr_Texture_features);

Shape_features=AUCRankFeature(Shape_features_index,2:end);
Corr_Shape_features=corr(Shape_features');
figure;imagesc(Corr_Shape_features);

Features= [Intensity_features ; Texture_features;Shape_features ]' ;

Corr_features=corr(Features);
figure;imshow(Corr_features);
figure;imagesc(Corr_features);
 
figure;plot(Corr_features(1,:));

k



% 
% % figure; colAUC(SUV_max_Pre',Response);
% Auc_SUV_max_decrease=colAUC(SUV_max_Pre',Response);
% 
% TestData=[SUV_max_Pre',Response]; 
% % figure;roc(TestData);
% x=TestData; [p_SUV_max,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% 
% PValue=strcat('P=', num2str(p));  AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), PValue);
%  
% 
% 
% 
% 
% %        
% %     AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), num2str(AUClist(k)));
% %     % Test P-value
% %     x=[AUCRankFeature(k,2:end)' Response];
% %     [p,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% %     
% %     %
% %     PValue=strcat('P=', num2str(p));
% %     AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), PValue);
% 
% 
% 
% 
% 
% 
% SUV_max_decrease=SUV_max_Pre-SUV_max_Post;
% SUV_max_ratio=SUV_max_Post./SUV_max_Pre;
% 
% SUV_peak_decrease=SUV_Peakmean_pre-SUV_Peakmean_Post;
% SUV_peak_ratio=SUV_Peakmean_Post./SUV_Peakmean_pre;
% 
% %
% % figure; colAUC(SUV_max_decrease',Response);
% Auc_SUV_max_decrease=colAUC(SUV_max_decrease',Response);
% 
% TestData=[SUV_max_decrease',Response]; 
% % figure;roc(TestData);
% x=TestData;
% [p_SUV_max_decrease,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% 
% %
% % figure; colAUC(SUV_max_ratio',Response);
% Auc_SUV_max_ratio=colAUC(SUV_max_ratio',Response);
% TestData=[SUV_max_ratio',Response]; 
% % figure;roc(TestData);
% 
%  x=TestData;
%  [p_SUV_max_ratio,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% 
% 
% %
% % figure; colAUC(SUV_peak_decrease',Response);
% Auc_SUV_peak_decrease=colAUC(SUV_peak_decrease',Response);
% 
% TestData=[SUV_peak_decrease',Response]; 
% % figure;roc(TestData);
% 
% x=TestData;
%  [p_SUV_peak_decrease,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% 
% 
% 
% %
% % figure; colAUC(SUV_peak_ratio',Response);
% Auc_SUV_peak_ratio=colAUC(SUV_peak_ratio',Response);
% 
% 
% TestData=[SUV_peak_ratio' Response];
% % figure; roc(TestData);
% 
% x=TestData;
%  [p_SUV_peak_ratio,h]=ranksum(x(x(:,2)==1),x(x(:,2)==0));
% 
% 
%  
%   PValue=strcat('P=', num2str(p));
%     AUCRandofFeatureName(k)=strcat(AUCRandofFeatureName(k), PValue);
%  
%  


% %% ROC using ColAUC (withoutoptimal threshold)
%    j=1;
%    TestData=AUCRankFeature(j,2:end);
%    figure;
%    colAUC(TestData',Response);
% %    FeatName=AUCRandofFeatureName(j);
% %    title([ FeatName 'ROC']);
%  

%% ROC using roc (with optimal threshold) 
   j=16;
   Test=AUCRankFeature(j,2:end);
   TestData=[Test' Response]; figure;
   roc(TestData);
   
   
   j=5;
   Test=TraditionalFeatures(j,:);
   TestData=[Test' Response]; figure;
   roc(TestData);
   
   
%    FeatName=AUCRandofFeatureName(j);
%    title([ FeatName 'ROC']);


%% Drawing Trence analysis

Test=AUCRankFeature(j,2:end);

% Test= SUV_max_decrease;
% 
% Test=SUV_max_ratio;
% 
% Test=SUV_peak_decrease;
% 
% Test=SUV_peak_ratio;


x=[Test' Response];

lu=length(x(x(:,2)==1)); %number of responder subjects
lh=length(x(x(:,2)==0)); %number of non-responder subjects
z=sortrows(x,1);
%find unique values in z
labels=unique(z(:,1));
ll=length(labels); %count unique value
a=zeros(ll,2); %array preallocation
PatientsOverThresh=zeros(ll,2);
NumofResponse=zeros(ll,2);
ubar=mean(x(x(:,2)==1),1); %responder mean value
hbar=mean(x(x(:,2)==0),1); %non-responder mean value
for K=1:ll
    if hbar<ubar
        TP=length(x(x(:,2)==1 & x(:,1)>labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)>labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)<=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)<=labels(K)));
    else
        TP=length(x(x(:,2)==1 & x(:,1)<labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)<labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)>=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)>=labels(K)));
    end
    a(K,:)=[TP/(TP+FN) TN/(TN+FP)]; %Sensitivity and Specificity
    PatientsOverThresh(K,:)=TP+FP;
    NumofResponse(K,:)=TP;
end

proportion = NumofResponse(:,1) ./ PatientsOverThresh(:,1);
[logitCoef,dev] = glmfit(labels,[NumofResponse(:,1) PatientsOverThresh(:,1)],'binomial','logit');
logitFit = glmval(logitCoef,labels,'logit');
figure;plot(labels,proportion,'bs', labels,logitFit,'r-','LineWidth',3);
ax1 = gca; set(ax1,'fontsize',16,'LineWidth',3)
xlabel(AUCRandofFeatureName(j),'FontSize',16,'FontWeight', 'bold'); ylabel('Positive Predictive Ratio','FontSize',16,'FontWeight', 'bold');
xlabel('Ratio of SUV-peak','FontSize',16,'FontWeight', 'bold'); ylabel('Positive Predictive Ratio','FontSize',16,'FontWeight', 'bold');

% 
% %% boxplot for each feature
% 
% 
% for j=1:20
%     j=2;
%     Test=AUCRankFeature(j,2:end);AUCRandofFeatureName(k)
%     x=[Test' Response];
%     
%     figure;
%     boxplot(x(x(:,2)==1));hold on;
%     boxplot(x(x(:,2)==0));
%     
%     lu=length(x(x(:,2)==1)); %number of responder subjects
%     lh=length(x(x(:,2)==0)); %number of non-responder subjects
%     z=sortrows(x,1);
%     %find unique values in z
%     labels=unique(z(:,1));
%     ll=length(labels); %count unique value
%     a=zeros(ll,2); %array preallocation
%     PatientsOverThresh=zeros(ll,2);
%     NumofResponse=zeros(ll,2);
%     ubar=mean(x(x(:,2)==1),1); %responder mean value
%     hbar=mean(x(x(:,2)==0),1); %non-responder mean value
%     for K=1:ll
%         if hbar<ubar
%             TP=length(x(x(:,2)==1 & x(:,1)>labels(K)));
%             FP=length(x(x(:,2)==0 & x(:,1)>labels(K)));
%             FN=length(x(x(:,2)==1 & x(:,1)<=labels(K)));
%             TN=length(x(x(:,2)==0 & x(:,1)<=labels(K)));
%         else
%             TP=length(x(x(:,2)==1 & x(:,1)<labels(K)));
%             FP=length(x(x(:,2)==0 & x(:,1)<labels(K)));
%             FN=length(x(x(:,2)==1 & x(:,1)>=labels(K)));
%             TN=length(x(x(:,2)==0 & x(:,1)>=labels(K)));
%         end
%         a(K,:)=[TP/(TP+FN) TN/(TN+FP)]; %Sensitivity and Specificity
%         
%         PatientsOverThresh(K,:)=TP+FP;
%         NumofResponse(K,:)=TP;
%     end
%     
%     
%     proportion = NumofResponse(:,1) ./ PatientsOverThresh(:,1);
%     
%     [logitCoef,dev] = glmfit(labels,[NumofResponse(:,1) PatientsOverThresh(:,1)],'binomial','logit');
%     logitFit = glmval(logitCoef,labels,'logit');
%     figure;plot(labels,proportion,'bs', labels,logitFit,'r-','LineWidth',3);
%     ax1 = gca; set(ax1,'fontsize',16,'LineWidth',3)
%     xlabel(AUCRandofFeatureName(j),'FontSize',16,'FontWeight', 'bold'); ylabel('Positive Predictive Ratio','FontSize',16,'FontWeight', 'bold');
%     
%     xlabel('Bounding Box Volume of ROI-SUV2.5','FontSize',16,'FontWeight', 'bold'); ylabl('Positive Predictive Ratio','FontSize',16,'FontWeight', 'bold');
%     
% end
% 







    function [SubFeature,SubFeatureName]=FeatureOrganzation2(folder)
        
        oldfolder=cd(folder);
        A=[];
        B=[];
        C=[];
        SubFeature=[];
        SubFeatureName=[];
        
        %
        [subfilename]=  'PreSUVFixedWindowMaskImageFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        
        FeatureSize=length(B);
        
        [subfilename]=  'PreSUVMaskBySegmentationSUVfixedFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        
        %
        [subfilename]=  'SUVDiffImageFixedWindowMaskImageFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        %
        [subfilename]=  'SUVDiffImageMaskBySegmentationSUVfixedFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        
        
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        
        %
        [subfilename]=  'SUVRegisteredFixedWindowMaskImageFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        
        %   [subfilename]=  'SUVRegisteredMaskBySegmentationSUVRegisteredFeatures.txt';
        %   [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        %   if length(A)==FeatureSize,
        %         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        %   SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        %   else
        %       [C, B] = textread([folder '/' 'SUVRegisteredFixedWindowMaskImageFeatures.txt'], '%s %f', -1);
        %   [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        %   SubFeatureName=[SubFeatureName; strcat(name,C)];   SubFeature=[SubFeature; B];
        %   end
        
        %
        [subfilename]=  'SUVRegisteredMaskBySegmentationSUVfixedFeatures.txt';
        [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
        [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
        SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
        
        
%         
%         [subfilename]=  'CTDiffImageFixedWindowMaskImageFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
%         [subfilename]=  'CTDiffImageMaskBySegmentationSUVfixedFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
%         %
%         [subfilename]=  'CTFixedToFixedSUVFixedWindowMaskImageFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
%         [subfilename]=  'CTFixedToFixedSUVMaskBySegmentationSUVfixedFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
%         %
%         [subfilename]=  'CTRegisteredToFixedSUVFixedWindowMaskImageFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
%         
%         %    [subfilename]=  'CTRegisteredToFixedSUVMaskBySegmentationSUVRegisteredFeatures.txt';
%         %   [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         %   if length(A)==FeatureSize,
%         %         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         %   SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         %   else
%         %       [C, B] = textread([folder '/' 'CTRegisteredToFixedSUVFixedWindowMaskImageFeatures.txt'], '%s %f', -1);
%         %   [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         %   SubFeatureName=[SubFeatureName; strcat(name,C)];   SubFeature=[SubFeature; B];
%         %   end
%         
%         %
%         [subfilename]=  'CTRegisteredToFixedSUVMaskBySegmentationSUVfixedFeatures.txt';
%         [fullfilename]=[folder '/' subfilename];  [A, B] = textread(fullfilename, '%s %f', -1);
%         [pathstr, name, ext] = fileparts(subfilename);    name=strcat(name,'_');
%         SubFeatureName=[SubFeatureName; strcat(name,A)];   SubFeature=[SubFeature; B];
%         
        cd(oldfolder);
        
   
   
