function SUV=calc_suv_Philips(dicomhd,slice)
% SUV=calc_suv(dicomhd,slice)
%
% Calcualtes the SUV*100 for the passed dicomHeader (dicomhd) and PET image
% slice(slice)
%
%Written IEN
%Modified for SUV*100 WL 3/27/09

SUV_SCALE = 100; %WL% multiply the SUV by 100 so that SUV can be represented as 2-bytes integer type


% Get calibration factor which is the Rescale slope Attribute Name in DICOM
m = dicomhd.RescaleSlope;
f = str2num(sprintf('%s', dicomhd.Private_7053_1000));
% f = double (dicomhd.Private_7053_1000(1));
b = dicomhd.RescaleIntercept; 


SUV= (slice * m + b) * f; %WL%

% SUV = max(0, min(round(SUV*SUV_SCALE), intmax('int16'))); %WL%
SUV = round(SUV*SUV_SCALE);
SUV = int16(SUV); %WL%



return