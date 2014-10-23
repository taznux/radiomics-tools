% Convert all dicom PET iamges into SUV*100


%clc; clear; close all;

% %%%% For 1701386 
% 
% % [dicom_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1701386/PrePET/04281545/41571870';
% % [dest_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1701386/1701386OUTPUT/PreSUV';
% 
% [dicom_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1701386/PostPET/29246687/07387738';
% [dest_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1701386/1701386OUTPUT/PostSUV';
% 

[dicom_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1602486/PrePET/download20110512112839/77104478/36354565';
[dest_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/PETCTDATA/1602486/1602486OUTPUT/PreSUV';

[dicom_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/NewPETCTDATA/PET001/PrePET';
[dest_dir] = '/data/research/lu_lab/work/wchoi/DATA/Patients/NewPETCTDATA/PET001/PET001OUTPUT/PreSUV';
         

% getthedir = dir([dicom_dir '/*DCM']);  % If all images end with .img
getthedir = dir([dicom_dir '/*dcm']); 
getthedir = dir([dicom_dir '/IM*']); 

endofthem = length(getthedir);
for i=1:endofthem;
% for i=1:2;
    src_temp = getthedir(i).name;
  %  info = dicominfo(src_temp);
    info = dicominfo([dicom_dir '/' src_temp]);
    
    %fprintf(1,'%d %s %d %d %f %f\n', i, src_temp, info.InstanceNumber, info.Private_7053_1000(1), info.RescaleSlope, info.RescaleIntercept);    
    temp = dicomread(info);    
    SUV_100 = calc_suv_100(info,double(temp));
%    SUV_100 = calc_suv_Philips_100(info,double(temp));
    info.RescaleSlope = 1.0;
    info.Modality = 'PT_SUV100';
    dest_temp = strcat(dest_dir, '/', src_temp);
    dicomwrite(SUV_100, dest_temp, info, 'CreateMode', 'copy');

end;
i

