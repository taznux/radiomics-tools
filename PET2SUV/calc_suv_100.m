function SUV=calc_suv(dicomhd,slice)
% SUV=calc_suv(dicomhd,slice)
%
% Calcualtes the SUV*100 for the passed dicomHeader (dicomhd) and PET image
% slice(slice)
%
%Written IEN
%Modified for SUV*100 WL 3/27/09

SUV_SCALE = 100; %WL% multiply the SUV by 100 so that SUV can be represented as 2-bytes integer type

% Get Patient weight in grams
ptweight=dicomhd.PatientWeight*1000;  % in grams

% Get Scan time
scantime=dcm_hhmmss(dicomhd.AcquisitionTime);
% Get calibration factor which is the Rescale slope Attribute Name in DICOM
calibration_factor=dicomhd.RescaleSlope;
intercept=dicomhd.RescaleIntercept; 


% Start Time for the Radiopharmaceutical Injection
injection_time=dcm_hhmmss(dicomhd.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime);
% Half Life for Radionuclide
half_life=dicomhd.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;
% Total dose injected for Radionuclide
injected_dose=dicomhd.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;

% Calculate the decay
%   decayFactor = e^(t1-t2/halflife);
decay=exp(-log(2)*(scantime-injection_time)/half_life);
%Calculate the dose decayed during procedure
injected_dose_decay=injected_dose*decay; % in Bq

% Calculate SUV.
% SUV  = (2DSlice x Calibration factor x Patient Weight) / Dose after decay
SUV=slice*calibration_factor*ptweight/injected_dose_decay; %WL%

% SUV = max(0, min(round(SUV*SUV_SCALE), intmax('int16'))); %WL%
SUV = round(SUV*SUV_SCALE);
SUV = int16(SUV); %WL%

return