% interrun variation (Evaluating the Robustness of segmentation methods to the ROI variation)

%All DSI values should be centralized into a matrix named DSI, in which the columns represent different ROIs and rows represent different images of a dataset.

%The parameters N and Q in each M-file shoule be modified according to the actual number of images and ROI.

clc;
clear all; 
close all;

N=20;   %%%%%%%%%%%%the number of images (for example, if the dataset included twenty patients, then N=20 )
Q=8;     %%%%%%%%%%%%% the number of ROI

Number=N*(Q*(Q-1)/2);
SumDifference=0;


load DSI  %%%%%load DSI data (In the matrix of DSI, columns represent different ROI and row represent different images of a dataset.)


for i=1:N
    for q=2:Q
        for p=1:q-1
            SumDifference=SumDifference+abs((DSI(i,p)-DSI(i,q))/(max(DSI(i,p),DSI(i,q))+0.000001));
        end
    end
end

InterrunVariation=SumDifference*100/Number;








            
