#!/usr/bin/python

import sys

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) < 3:
    print sys.argv[0], ' [input dir] [output dir]'
    sys.exit(-1)
