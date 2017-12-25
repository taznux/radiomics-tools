#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script file for organizing features from each individual into a single csv file.
"""

import os
import sys
import glob

parameter_list = ['No.', 'PID', 'Image', 'Tags']
shape_feature = ['PhysicalSize', 'BoundingBoxVolume', 'BoundingBoxSize1', 'BoundingBoxSize2', 'BoundingBoxSize3',
                 'OrientedBoundingBoxVolume', 'OrientedBoundingBoxSize1', 'OrientedBoundingBoxSize2',
                 'OrientedBoundingBoxSize3', 'EquivalentEllipsoidDiameter1', 'EquivalentEllipsoidDiameter2',
                 'EquivalentEllipsoidDiameter3', 'EquivalentSphericalPerimeter', 'EquivalentSphericalRadius',
                 'FeretDiameter', 'NumberOfLines', 'NumberOfPixels', 'Perimeter', 'PrincipalAxes1', 'PrincipalAxes2',
                 'PrincipalAxes3', 'PrincipalAxes4', 'PrincipalAxes5', 'PrincipalAxes6', 'PrincipalAxes7',
                 'PrincipalAxes8', 'PrincipalAxes9', 'PrincipalMoments1', 'PrincipalMoments2', 'PrincipalMoments3',
                 'Eccentricity', 'Elongation', 'Flatness', 'Orientation', 'Roundness']
shape_intensity_feature = ['WeightedElongation', 'WeightedFlatness', 'WeightedPrincipalAxes1', 'WeightedPrincipalAxes2',
                           'WeightedPrincipalAxes3', 'WeightedPrincipalAxes4', 'WeightedPrincipalAxes5',
                           'WeightedPrincipalAxes6', 'WeightedPrincipalAxes7', 'WeightedPrincipalAxes8',
                           'WeightedPrincipalAxes9', 'WeightedPrincipalMoments1', 'WeightedPrincipalMoments2',
                           'WeightedPrincipalMoments3']
intensity_feature = ['Kurtosis', 'Maximum', 'Mean', 'Median', 'Minimum', 'Skewness', 'StandardDeviation', 'Sum',
                     'Variance']
shape_feature_2d = ['2DPhysicalSize', '2DBoundingBoxVolume', '2DBoundingBoxSize1', '2DBoundingBoxSize2',
                    '2DOrientedBoundingBoxVolume', '2DOrientedBoundingBoxSize1', '2DOrientedBoundingBoxSize2',
                    '2DEquivalentEllipsoidDiameter1', '2DEquivalentEllipsoidDiameter2',
                    '2DEquivalentSphericalPerimeter', '2DEquivalentSphericalRadius', '2DFeretDiameter',
                    '2DNumberOfLines', '2DNumberOfPixels', '2DPerimeter', '2DPrincipalAxes1', '2DPrincipalAxes2',
                    '2DPrincipalAxes3', '2DPrincipalAxes4', '2DPrincipalMoments1', '2DPrincipalMoments2',
                    '2DEccentricity', '2DElongation', '2DFlatness', '2DOrientation', '2DRoundness']
shape_intensity_feature_2d = ['2DWeightedElongation', '2DWeightedFlatness', '2DWeightedPrincipalAxes1',
                              '2DWeightedPrincipalAxes2', '2DWeightedPrincipalAxes3', '2DWeightedPrincipalAxes4',
                              '2DWeightedPrincipalMoments1', '2DWeightedPrincipalMoments2']
intensity_feature_2d = ['2DKurtosis', '2DMaximum', '2DMean', '2DMedian', '2DMinimum', '2DSkewness',
                        '2DStandardDeviation', '2DSum', '2DVariance']
texture_feature = ['MeanOfEnergy', 'MeanOfEntropy', 'MeanOfCorrelation', 'MeanOfInverseDifferenceMoment',
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
                   'StandardDeviationOfHighGreyLevelRunEmphasis', 'StandardDeviationOfShortRunLowGreyLevelEmphasis',
                   'StandardDeviationOfShortRunHighGreyLevelEmphasis', 'StandardDeviationOfLongRunLowGreyLevelEmphasis',
                   'StandardDeviationOfLongRunHighGreyLevelEmphasis']

feature_list = shape_feature + shape_intensity_feature + intensity_feature + shape_feature_2d + shape_intensity_feature_2d + intensity_feature_2d + texture_feature


def find(inputlist, value):
    """Find target feature in feature list

    :param inputlist: feature list
    :param value: feature name
    :return: target feature index
    """
    try:
        idx = inputlist.index(value)
    except ValueError:
        idx = -1

    return idx

def feature_parsing(filename, in_feature_list):
    """Parsing feature text to value, vector, or matrix

    :param filename: a feature text filename from feature extraction code
    :param in_feature_list: target feature list
    :return: parsed feature values
    """
    values = [''] * len(in_feature_list)
    try:
        f = open(filename, 'r')
    except IOError:
        return values
    row_idx = 0
    matrix = [[], [], []]
    for line in f:
        strs = line.split('=')
        if len(strs) > 1:
            feature_name = strs[0].strip()
            value_str = strs[1].strip()
        else:
            value_str = strs[0].strip()

        if len(value_str) == 0 or feature_name == 'Histogram':
            continue
        if feature_name[0] == '2':
            matrix_size = 2
        else:
            matrix_size = 3

        if value_str[0] == '[':  # vector
            vector = []
            idx = 1
            for v in value_str[1:-1].split(','):
                new_feature_name = feature_name + str(idx)
                v_str = v.strip()
                value = v_str
                fidx = find(in_feature_list, new_feature_name)
                if fidx > -1:
                    values[fidx] = value
                vector.append(value)
                idx += 1
                print(new_feature_name + "=" + value)

        elif value_str.find(' ')>0:  # matrix
            if row_idx == 0:
                idx = 1
                matrix[0] = []
                matrix[1] = []
                matrix[2] = []
            for v in value_str.split(' '):
                new_feature_name = feature_name + str(idx)
                v_str = v.strip()
                value = v_str
                fidx = find(in_feature_list, new_feature_name)
                if fidx > -1:
                    values[fidx] = value
                matrix[row_idx].append(value)
                print(new_feature_name + "=" + value)
                idx += 1
            row_idx += 1
            if row_idx == matrix_size:
                row_idx = 0
        else:
            value = value_str

            fidx = find(in_feature_list, feature_name)
            if fidx > -1:
                values[fidx] = value
            print(feature_name + "=" + value)

    return values


def convert_name_to_code(prefix, in_feature_list):
    """
    Feature name conversion from original name to code
    T01, T02 ... T99

    :param prefix: code prefix, eg. T for texture, S for shape, or I for intensity
    :param in_feature_list: target feature list
    :return: coded feature names
    """
    idx = []
    for i in range(len(in_feature_list)):
        idx.append(prefix + '%02d' % (i + 1))

    return idx


def organize(input_files, output, parameter_list, feature_list):
    """Feature organization from feature text files to a single csv file for analysis

    :param input_files: feature text files in list
    :param output: csv file path
    :param parameter_list: patient parameter list, eg. PID and No.
    :param feature_list: target feature list
    """
    fw = open(output, 'w')

    column_names = ','.join(parameter_list + feature_list)
    fw.write(column_names + '\n')

    print(input_files)

    no = 1
    for f in input_files:
        if f.find('.txt') < 0:
            continue

        print(f)
        path = os.path.dirname(f)
        strs = os.path.basename(f).replace(".txt", "").split("_")
        pid = strs[0]
        image = strs[1]


        fw.write(str(no) + ',')
        fw.write(pid + ',')
        fw.write(image + ',')
        if len(parameter_list) > 3:
            tags = '_'.join(strs[2:len(strs)])
            fw.write(tags + ',')
        no += 1

        file_name = f
        values = ','.join(feature_parsing(file_name, feature_list))
        fw.write(values+'\n')

    fw.close()


if __name__ == "__main__":
    if len(sys.argv) > 2:
        input_files = sys.argv[1]
        output = sys.argv[2]
        if len(sys.argv) > 4:
            parameter_list = sys.argv[3]
            feature_list = sys.argv[4]
    else:
        os.mkdir("output")
        input_files = glob.glob("output/*.txt")
        output = "output/feature_list.csv"

    organize(input_files, output, parameter_list, feature_list)
