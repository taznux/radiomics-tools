#!/usr/bin/python

import os
import sys

data_list = os.listdir('../DATA/NLST')
print data_list

for data in data_list[2:]:
    print data
    os.system('./DICOMReadAndWriteNLST.py ../DATA/NLST/'+data+' ../DATA/Patients/NLST_nrrd/'+data)


