#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is for image boolean calculation
"""

import SimpleITK as sitk
import sys


def STAPLE_to_labelmap(input_image, th=0.5):
    output_image = sitk.BinaryThreshold(input_image, th, 255.0)
    return output_image


if __name__ == "__main__":
    th = 0.5
    if len(sys.argv) < 2:
        print(sys.argv[0] + " input output [threshold]\n")
    else:
        inputFilename = sys.argv[1]
        outputFilename = sys.argv[2]
        if len(sys.argv) == 4:
            th = float(sys.argv[3])

        inputImage = sitk.ReadImage(inputFilename)
        outputImage = STAPLE_to_labelmap(inputImage, th)

        sitk.WriteImage(outputImage, outputFilename, True)
