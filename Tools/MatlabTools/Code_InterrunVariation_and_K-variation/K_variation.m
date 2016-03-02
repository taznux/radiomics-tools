
%the generalized Cohen¡¯s kappa (K) (Evaluating the Robustness of segmentation methods to the ROI variation)

%All DSI values should be centralized into a matrix named DSI, in which the columns represent different ROIs and rows represent different images of a dataset.

%The parameters N and Q in each M-file shoule be modified according to the actual number of images and ROI.

clc; 
clear all; 
close all;

N=20;    %%%%%%%%%%%%the number of images (for example, if the dataset included twenty patients, then N=20 )
Q=8;   %%%%%%%%%%%%%the number of ROI

Number1=N*(Q*(Q-1)/2);
Number2=N*N*(Q*(Q-1)/2);

SumDifference1=0;   %%%%%
SumDifference2=0; 
% K=0;
% delta=0;
% mu=0;


load DSI  %%%%%load DSI data (In the matrix of DSI, columns represent different ROI and row represent different images of a dataset.)


for i=1:N
    for q1=2:Q
        for p1=1:q1-1
             SumDifference1=SumDifference1+abs(DSI(i,p1)-DSI(i,q1));
        end 
    end
    
    for j=1:N
        for q2=2:Q
        for p2=1:q2-1
             SumDifference2=SumDifference2+abs(DSI(i,p2)-DSI(j,q2));
        end 
        end
    end
end

delta=SumDifference1/Number1;
mu=SumDifference2/Number2;

K=1-delta/mu;




    
            
    