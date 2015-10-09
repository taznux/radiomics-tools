#!/usr/bin/python

import os
import sys


def recur_dir(input_dir,output_dir):
  input_list = os.listdir(input_dir)
  nodir = True;
  for input in input_list:
    if os.path.isdir(input) == True:
      recur_dir(input)
      nodir = False
  if nodir == True:
    cmd = 'DICOM2NRRDConverter '+input_dir+' '+output_dir
    print cmd
    os.system(cmd)


#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

if len(sys.argv) < 3:
  print sys.argv[0], ' [input dir] [output dir]'
  sys.exit(-1)

  recur_dir(sys.argv[1],sys.argv[2])
