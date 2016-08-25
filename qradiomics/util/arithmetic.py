#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Four arithmetic operations for images
"""

import SimpleITK as sitk
import sys

from qradiomics.io import *


def image_add(input_image1, input_image2):
    """Add operation between two input images: input_image1 + input_image2

    :param input_image1: first input image
    :param input_image2: second input image
    :return: result image
    """
    return sitk.Add(input_image1, input_image2)


def image_sub(input_image1, input_image2):
    """Subtracting second input image from first input image: input_image1 - input_image2

    :param input_image1: first input image
    :param input_image2: second input image
    :return: result image
    """
    return sitk.Subtract(input_image1, input_image2)


def image_div(input_image1, input_image2):
    """Dividing first input image by second input image: input_image1 / input_image2

    :param input_image1: first input image
    :param input_image2: second input image
    :return: result image
    """
    return sitk.Divide(input_image1, input_image2)


def image_mul(input_image1, input_image2):
    """Multiplying first input image and second input image: input_image1 * input_image2

    :param input_image1: first input image
    :param input_image2: second input image
    :return: result image
    """
    return sitk.Multiply(input_image1, input_image2)


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(sys.argv[0] + " [add|sub|div|mul] input1 input2 output")
    else:
        operator = sys.argv[1]
        inputFilename1 = sys.argv[2]
        inputFilename2 = sys.argv[3]
        outputFilename = sys.argv[4]

        inputImage1 = image_read(inputFilename1)
        inputImage2 = image_read(inputFilename2)

        if operator == "add":
            outputImage = image_add(inputImage1, inputImage2)
        elif operator == "sub":
            outputImage = image_sub(inputImage1, inputImage2)
        elif operator == "div":
            outputImage = image_div(inputImage1, inputImage2)
        elif operator == "mul":
            outputImage = image_mul(inputImage1, inputImage2)
        else:
            print("Not supporting operator " + sys.argv[1])
            exit(-1)

        image_write(outputImage, outputFilename)
