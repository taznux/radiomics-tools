#!/usr/bin/env python

import os
import os.path as osp
import sys
import subprocess
import multiprocessing
from multiprocessing import freeze_support
from ruffus import *

sys.path.append(".")
sys.path.append("..")

from qradiomics import *
from qradiomics.io import *
from qradiomics.util import organize_features as organize_features

################################################################################
# Set your environmental parameters
################################################################################
experiment_set = 'TrainingSet'
#experiment_set = 'TestSet'
output_path = 'output'
data_path = 'DATA'
dicom_path = data_path + '/DOI'
image_path = data_path + '/' + experiment_set
nodule_info_path = './' + experiment_set + '.csv'
################################################################################

metadata.load_metadata(nodule_info_path)

parameter_list = ['No.', 'PID', 'name']
shape_feature_list = ['2DEquivalentEllipsoidDiameter1', '2DEquivalentEllipsoidDiameter2', 'Volume', 'Centroid1',
                      'Centroid2', 'Centroid3', 'AxesLength1', 'AxesLength2', 'AxesLength3', 'MajorAxisLength',
                      'MinorAxisLength', 'Eccentricity', 'Elongation', 'Orientation', 'BoundingBoxVolume',
                      'BoundingBoxSize1', 'BoundingBoxSize2', 'BoundingBoxSize3', 'OrientedBoundingBoxVolume',
                      'OrientedBoundingBoxSize1', 'OrientedBoundingBoxSize2', 'OrientedBoundingBoxSize3',
                      'Eigenvalues1', 'Eigenvalues2', 'Eigenvalues3', 'Eigenvectors1', 'Eigenvectors2', 'Eigenvectors3',
                      'Eigenvectors4', 'Eigenvectors5', 'Eigenvectors6', 'Eigenvectors7', 'Eigenvectors8',
                      'Eigenvectors9']
intensity_feature_list = ['Minimum', 'Maximum', 'Median', 'Mean', 'Variance', 'Sum', 'StandardDeviation', 'Skewness',
                          'Kurtosis', ]
texture_feature_list = ['MeanOfEnergy', 'MeanOfEntropy', 'MeanOfCorrelation', 'MeanOfInverseDifferenceMoment',
                        'MeanOfInertia', 'MeanOfClusterShade', 'MeanOfClusterProminence', 'MeanOfHaralickCorrelation',
                        'StandardDeviationOfEnergy', 'StandardDeviationOfEntropy', 'StandardDeviationOfCorrelation',
                        'StandardDeviationOfInverseDifferenceMoment', 'StandardDeviationOfInertia',
                        'StandardDeviationOfClusterShade', 'StandardDeviationOfClusterProminence',
                        'StandardDeviationOfHaralickCorrelation', 'MeanOfShortRunEmphasis', 'MeanOfLongRunEmphasis',
                        'MeanOfGreyLevelNonuniformity', 'MeanOfRunLengthNonuniformity', 'MeanOfLowGreyLevelRunEmphasis',
                        'MeanOfHighGreyLevelRunEmphasis', 'MeanOfShortRunLowGreyLevelEmphasis',
                        'MeanOfShortRunHighGreyLevelEmphasis', 'MeanOfLongRunLowGreyLevelEmphasis',
                        'MeanOfLongRunHighGreyLevelEmphasis', 'StandardDeviationOfShortRunEmphasis',
                        'StandardDeviationOfLongRunEmphasis', 'StandardDeviationOfGreyLevelNonuniformity',
                        'StandardDeviationOfRunLengthNonuniformity', 'StandardDeviationOfLowGreyLevelRunEmphasis',
                        'StandardDeviationOfHighGreyLevelRunEmphasis',
                        'StandardDeviationOfShortRunLowGreyLevelEmphasis',
                        'StandardDeviationOfShortRunHighGreyLevelEmphasis',
                        'StandardDeviationOfLongRunLowGreyLevelEmphasis',
                        'StandardDeviationOfLongRunHighGreyLevelEmphasis']
feature_list = shape_feature_list + intensity_feature_list + texture_feature_list

# paths for tools
dicom2nrrd_converter = osp.abspath('bin/DICOM-RT2NRRDConverter')
nodule_segmentation = osp.abspath('bin/NoduleSegmentation')
feature_extraction = osp.abspath('bin/FeatureExtraction')


def task_load_dicom_list():
    # retrieve dicom dirs and preparation of the input image sets
    patient_path_list = [osp.join(dicom_path, fn) for fn in next(os.walk(dicom_path))[1]]

    for pt_dicom_path in patient_path_list:
        pid = osp.basename(pt_dicom_path)

        pt = metadata.get_patient(pid)
        if len(pt) == 0:
            continue

        print(os.listdir(pt_dicom_path)[0])
        study = os.listdir(pt_dicom_path)[-1]
        series = os.listdir(osp.join(pt_dicom_path, study))[-1]

        input_file = osp.join(pt_dicom_path, study, series)
        output_file = osp.join(image_path, pid + ".nrrd")

        print(input_file, output_file)
        yield [input_file, output_file]


def task_dicom_to_nrrd_convert(input_file, output_file):
    image_write(image_read(input_file), output_file)


def task_nodule_segmentation(input_file, output_files, pid, output_prefix):
    # # IMPORTANT: cleanup rubbish from previous run first
    # for oo in output_files:
    #     os.unlink(oo)

    pt = metadata.get_patient(pid)

    for idx in pt:
        pt_row = pt[idx]
        x = pt_row['X']
        y = pt_row['Y']
        slice_number = pt_row['Z']
        large_diameter = pt_row['LD']
        short_diameter = pt_row['PD']
        print(idx, slice_number, large_diameter, short_diameter)

        output_file = output_prefix + str(idx) + "-label.nrrd"
        open(output_file, "w")
        p = subprocess.Popen([nodule_segmentation, input_file, str(x), str(y), str(slice_number), str(large_diameter),
                              str(short_diameter), str(output_file)])
        p.wait()


def task_feature_extraction(input_file, output_file, mask_file):
    open(output_file, "w")
    p = subprocess.Popen([feature_extraction, input_file, mask_file, output_file])
    p.wait()


def task_feature_organization(input_files, output_file):
    organize_features.organize(input_files, output_file, parameter_list, feature_list)


def make_pipeline_lungx():
    pipeline_name = experiment_set
    pipeline = Pipeline(pipeline_name)

    # TODO : change to originate operator
    pipeline.files(task_dicom_to_nrrd_convert,
                   task_load_dicom_list) \
        .follows(mkdir(image_path))

    pipeline.subdivide(name="task_nodule_segmentation",
                       task_func=task_nodule_segmentation,
                       input=output_from("task_dicom_to_nrrd_convert"),
                       filter=formatter(),
                       output=output_path + "/" + experiment_set + "/{basename[0]}_*[0-9]-label.nrrd",
                       # after '_' will be index of nodule which is some digit
                       extras=["{basename[0]}",
                               output_path + "/" + experiment_set + "/{basename[0]}_"]) \
        .follows(mkdir(output_path), mkdir(output_path + "/" + experiment_set))

    pipeline.transform(name="task_feature_extraction",
                       task_func=task_feature_extraction,
                       input=output_from("task_nodule_segmentation"),
                       filter=formatter(),
                       output=output_path + "/" + experiment_set + "/{basename[0]}.txt",
                       extras=["{path[0]}/{basename[0]}.nrrd"])

    pipeline.merge(name="task_feature_organization",
                   task_func=task_feature_organization,
                   input=output_from("task_feature_extraction"),
                   output=osp.join(output_path, "feature_list_" + experiment_set + ".csv"))

    return pipeline


if __name__ == "__main__":
    freeze_support()

    pipeline_lungx = make_pipeline_lungx()

    # pipeline_lungx._printout_graph("flowchart.png")
    pipeline_lungx.printout(sys.stdout, verbose=6)
    pipeline_lungx.run(multiprocess=multiprocessing.cpu_count()-1)

# TODO : bugfix ruffus task.py      5774:           job_result = ii.next(timeout=999999990->9999)
