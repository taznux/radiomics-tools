# -*- coding: utf-8 -*-

"""
This module is metadata using pandas data frame
"""

import pandas

metadata = pandas.DataFrame()


def load_metadata(metadata_path):
    """load metadata file

    :param metadata_path: a path for metadata file
    """
    global metadata
    # load metadata
    metadata = pandas.read_csv(metadata_path)
    try:
        metadata.set_index(['No.'], inplace=True)
    except KeyError:
        metadata.set_index(['No'], inplace=True)
    #print(metadata)


def get_patient(pid):
    """To access an individual patient's information with pid

    :param pid: unique patient ID
    :return: metadata of the patient in dictionary structure
    """
    res = ''
    #print(metadata)
    try:
        df = metadata.loc[metadata['PID'] == pid]
        res = df.to_dict('index')

        #print(res)

    except KeyError as e:
        print('KeyError: '+str(e))

    return res
