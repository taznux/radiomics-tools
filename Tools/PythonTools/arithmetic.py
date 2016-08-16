#!/usr/bin/env python

import SimpleITK as sitk
import sys

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(sys.argv[0] + " [add|sub|div|mul] input1 input2 output")
    else:
        operator = sys.argv[1]
        inputFilename1 = sys.argv[2]
        inputFilename2 = sys.argv[3]
        outputFilename = sys.argv[4]

        if operator == "add":
            opFilter = sitk.AddImageFilter()
        elif operator == "sub":
            opFilter = sitk.SubtractImageFilter()
        elif operator == "div":
            opFilter = sitk.DivideImageFilter()
        elif operator == "mul":
            opFilter = sitk.MultiplyImageFilter()
        else:
            print("Not supporting operator " + sys.argv[1])
            exit(-1)

        inputImage1 = sitk.ReadImage(inputFilename1)
        inputImage2 = sitk.ReadImage(inputFilename2)

        outputImage = opFilter.Execute(inputImage1, inputImage2)

        sitk.WriteImage(outputImage, outputFilename)
