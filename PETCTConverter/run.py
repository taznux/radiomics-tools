#!/home/wlui/bin/python

import os
import sys
import csv
from collections import defaultdict

columns = defaultdict(list)
data_path = '../DATA/Patients/PETCTDATA/'
output_path = 'output'
patient_info_path = data_path+'/patient.csv'

def data_dir_recursion(data_path):
    data_list = os.listdir(data_path)
    data_list.sort()
    print data_list
    for data in data_list:
        new_path = data_path+'/'+data
        if data.find('Post') > 0:
            os.system('./BIN/DICOMReadAndWrite '+new_path+' '+output_path)
        elif os.path.isdir(new_path):
            data_dir_recursion(new_path)
        elif data.find('dcm') > 0:
            return


def create_path_to_dir(path):
    if not os.path.exists(path):
        l=[]
        p = ""
        l = path.split("/")
        i = 0
        while i < len(l):
            p = p + l[i] + "/"
            i = i + 1
            if not os.path.exists(p):
                os.mkdir(p, 0755)



if __name__ == "__main__":
    data_dir_recursion(data_path)
