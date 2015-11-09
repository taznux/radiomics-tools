__author__ = 'wchoi'


import SimpleITK as sitk
import numpy as np
#import os

#from ipywidgets import interact, fixed

import registration_callbacks as rc
#import registration_utilities as ru


inputImageList = [
    ["UMD0003","UMD0003_20050121_PT","UMD0003_20050120_CT"],
    ["UMD0012","UMD0012_20050321_PT","UMD0012_20050321_CT"],
    ["UMD0053","UMD0053_20050203_PT","UMD0053_20050202_CT"]
]


#%matplotlib qt

# This is the registration configuration which we use in all cases. The only parameter that we vary
# is the initial_transform.
def multires_registration(fixed_image, moving_image, initial_transform):
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=100)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.04)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=200, estimateLearningRate=registration_method.EachIteration)
    registration_method.SetOptimizerScalesFromPhysicalShift()
    registration_method.SetInitialTransform(initial_transform)
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = [2,1,0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    registration_method.AddCommand(sitk.sitkStartEvent, rc.metric_start_plot)
    registration_method.AddCommand(sitk.sitkEndEvent, rc.metric_end_plot)
    registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, rc.metric_update_multires_iterations)
    registration_method.AddCommand(sitk.sitkIterationEvent, lambda: rc.metric_plot_values(registration_method))

    final_transform = registration_method.Execute(fixed_image, moving_image)
    print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
    print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))
    return final_transform



def save_transform_and_image(transform, fixed_image, moving_image, outputfile_prefix):
    """
    Write the given transformation to file, resample the moving_image onto the fixed_images grid and save the
    result to file.

    Args:
        transform (SimpleITK Transform): transform that maps points from the fixed image coordinate system to the moving.
        fixed_image (SimpleITK Image): resample onto the spatial grid defined by this image.
        moving_image (SimpleITK Image): resample this image.
        outputfile_prefix (string): transform is written to outputfile_prefix.tfm and resampled image is written to
                                    outputfile_prefix.mha.
    """
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(fixed_image)

    # SimpleITK supports several interpolation options, we go with the simplest that gives reasonable results.
    resample.SetInterpolator(sitk.sitkLinear)
    resample.SetTransform(transform)
    sitk.WriteImage(resample.Execute(moving_image), outputfile_prefix+'.mha')
    sitk.WriteTransform(transform, outputfile_prefix+'.tfm')


for inputImages in inputImageList:
    crop_str = "-subvolume-scale_1"
    print("Load Images")
    fixed_image = sitk.ReadImage("D:/WFUBMC_nrrd/"+inputImages[0]+"/"+inputImages[1]+crop_str+".nrrd", sitk.sitkFloat32)
    print(fixed_image)
    moving_image = sitk.ReadImage("D:/WFUBMC_nrrd/"+inputImages[0]+"/"+inputImages[2]+crop_str+".nrrd", sitk.sitkFloat32)
    print(moving_image)

    initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    print(initial_transform)

    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.02)
    registration_method.SetInterpolator(sitk.sitkLinear)
    # The order of parameters for the Euler3DTransform is [angle_x, angle_y, angle_z, t_x, t_y, t_z]. The parameter
    # sampling grid is centered on the initial_transform parameter values, that are all zero for the rotations. Given
    # the number of steps and their length and optimizer scales we have:
    # angle_x = -pi, 0, pi
    # angle_y = 0
    # angle_z = -pi, -pi/2, 0, pi/2, pi
    registration_method.SetOptimizerAsExhaustive(numberOfSteps=[1,0,2,0,0,0], stepLength = np.pi)
    registration_method.SetOptimizerScales([1,1,0.5,1,1,1])

    registration_method.AddCommand(sitk.sitkStartEvent, rc.metric_start_plot)
    registration_method.AddCommand(sitk.sitkEndEvent, rc.metric_end_plot)
    registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, rc.metric_update_multires_iterations)
    registration_method.AddCommand(sitk.sitkIterationEvent, lambda: rc.metric_plot_values(registration_method))

    #Perform the registration in-place so that the initial_transform is modified.
    registration_method.SetInitialTransform(initial_transform, inPlace=True)
    registration_method.Execute(fixed_image, moving_image)

    print(initial_transform)

    final_transform = multires_registration(fixed_image, moving_image, initial_transform)

    print(final_transform)
    final_transform.WriteTransform("D:/WFUBMC_nrrd/"+inputImages[0]+"/LinearTransform_5.h5")