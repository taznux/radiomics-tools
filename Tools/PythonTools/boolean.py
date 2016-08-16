#!/usr/bin/env python

import SimpleITK as sitk
import sys

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(sys.argv[0] + " [and|or|xor] input1 input2 output \n" +
              sys.argv[0] + " not input output")
    else:
        operator = sys.argv[1]
        inputFilename1 = sys.argv[2]
        inputFilename2 = sys.argv[3]
        outputFilename = sys.argv[4]

        if operator == "and":
            opFilter = sitk.AndImageFilter()
        elif operator == "or":
            opFilter = sitk.OrImageFilter()
        elif operator == "xor":
            opFilter = sitk.XorImageFilter()
        elif operator == "not": # unary operator
            opFilter = sitk.NotImageFilter()
        else:
            print("Not supporting operator " + sys.argv[1])
            exit(-1)

        inputImage1 = sitk.ReadImage(inputFilename1, sitk.sitkUInt16)
        if operator != "not":
            inputImage2 = sitk.ReadImage(inputFilename2, sitk.sitkUInt16)
            outputImage = opFilter.Execute(inputImage1, inputImage2)
        else:
            outputFilename = inputFilename2
            outputImage = opFilter.Execute(inputImage1)

        sitk.WriteImage(outputImage, outputFilename)
