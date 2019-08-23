#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script file for image resample with reference image.
"""

import SimpleITK as sitk
import sys
import numpy as np


def image_resample(input_image, reference, is_label):
    """

    Image resampling with reference information

    :param input_image: input image
    :param reference: reference image or reference spacing
    :param is_label: input image type - label map or typical image
    :return: resampled image
    """
    resample_filter = sitk.ResampleImageFilter()

    input_spacing = input_image.GetSpacing()
    input_origin = input_image.GetOrigin()
    input_size = input_image.GetSize()
    input_direction = input_image.GetDirection()

    if type(reference) == sitk.Image:
        ref_image = reference
        output_spacing = ref_image.GetSpacing()
        output_direction = ref_image.GetDirection()
        origin = ref_image.GetOrigin()
        #origin = np.asarray(ref_image.GetOrigin()) - np.asarray(output_spacing)/2 + np.asarray(input_spacing)/2
        dist = (np.asarray(input_origin) - origin) / np.asarray(input_spacing)
        close_interp = np.floor(
            dist * np.asarray(input_spacing) / np.asarray(output_spacing))
        output_origin = ref_image.GetOrigin() + close_interp * \
            np.asarray(output_spacing)
    elif type(reference) == list:
        output_spacing = reference
        output_direction = input_direction
        output_origin = input_origin
        #output_origin = np.asarray(input_origin) - np.asarray(input_spacing)/2 + np.asarray(output_spacing)/2

    output_size = np.ceil(np.asarray(
        input_size) * np.asarray(input_spacing) / np.asarray(output_spacing)).astype(int)

    print("output spacing: " + output_spacing.__str__())
    print("output origin: " + output_origin.__str__())
    print("output size: " + output_size.__str__())
    print("output direction: " + output_direction.__str__())

    resample_filter.SetOutputSpacing(output_spacing)
    resample_filter.SetOutputOrigin(output_origin)
    resample_filter.SetSize(output_size.tolist())
    resample_filter.SetOutputDirection(output_direction)

    if is_label:
        resample_filter.SetInterpolator(sitk.sitkNearestNeighbor)
        output_image = resample_filter.Execute(input_image)
    else:
        resample_filter.SetInterpolator(sitk.sitkLinear)
        output_image = resample_filter.Execute(input_image)

    return output_image


if __name__ == "__main__":
    is_label = False
    is_space = True

    try:
        float(sys.argv[3])
    except:
        is_space = False

    if len(sys.argv) < 4 or (is_space and len(sys.argv) == 5):
        print("Usage: ")
        print(sys.argv[0] +
              " InputImage OutputImage TargetSpacingX Y Z [is_label]")
        print(sys.argv[0] +
              " InputImage OutputImage ReferenceImage [is_label]")

        exit(-1)

    if len(sys.argv) > 6 or (not(is_space) and len(sys.argv) == 5):
        is_label = True

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    # load input image file
    print(input_filename)
    input_image = sitk.ReadImage(input_filename)
    input_pixel_ID = input_image.GetPixelID()

    input_spacing = input_image.GetSpacing()
    input_origin = input_image.GetOrigin()
    input_size = input_image.GetSize()
    input_direction = input_image.GetDirection()

    print("input spacing: " + input_spacing.__str__())
    print("input origin: " + input_origin.__str__())
    print("input size: " + input_size.__str__())
    print("input direction: " + input_direction.__str__())

    ref_spacing = list(input_image.GetSpacing())
    if is_space:
        ref_spacing[0] = float(sys.argv[3])
        ref_spacing[1] = float(sys.argv[4])
        ref_spacing[2] = float(sys.argv[5])
        output_image = image_resample(input_image, ref_spacing, is_label)
    else:
        ref_filename = sys.argv[3]
        print(ref_filename)
        ref_image = sitk.ReadImage(ref_filename)
        output_image = image_resample(input_image, ref_image, is_label)

    sitk.WriteImage(sitk.Cast(output_image, input_pixel_ID),
                    output_filename, True)
    #sitk.WriteImage(output_image, output_filename, True)
