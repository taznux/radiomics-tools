#!/usr/bin/python

import os
import sys

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

if len(sys.argv) < 3:
    print sys.argv[0], ' [input dir] [output dir]'
    sys.exit(-1)

input_list = os.listdir(sys.argv[1])

for input in input_list:
    cmd = './BIN/DICOM2NRRDConverter '+sys.argv[1]+'/'+input+' '+sys.argv[2]
    print cmd
    os.system(cmd)


