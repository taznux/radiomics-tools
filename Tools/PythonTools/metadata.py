# -*- coding: utf-8 -*-

"""
This module is metadata using pandas data frame
"""

import pandas

class metadata():
    data = None
    path = None


    def __init__(self,metadata_path=None):
        if metadata_path!=None:
            self.path = metadata_path
            self.load()


    def load(self):
        """
        load metadata file
        """
        # load metadata
        self.data = pandas.read_csv(self.path)
        try:
            self.data.set_index(['No.'], inplace=True)
        except KeyError:
            self.data.set_index(['No'], inplace=True)
        #print(self.data)

    def save(self):
        """
        save metadata file
        """
        # save metadata
        self.data.to_csv(self.path)
        #print(self.data)


    def getRows(self,column,kwd):
        """
        To access an individual patient's information with specific coulumn matching

        :param column: column name
        :param kwd: keyword
        :return: metadata of the patient in dictionary structure
        """
        res = ''
        #print(self.data)
        try:
            df = self.data.loc[self.data[column] == kwd]
            res = df.to_dict('index')

            #print(res)

        except KeyError as e:
            print('KeyError: '+str(e))

        return res


    def getPatient(self,pid):
        """
        To access an individual patient's information with pid

        :param pid: unique patient ID
        :return: metadata of the patient in dictionary structure
        """
        res = self.getRows('PID', pid)

        return res
    

    def getSeries(self,uid):
        """
        To access an individual series's information with uid

        :param pid: unique patient ID
        :return: metadata of the patient in dictionary structure
        """
        res = self.getRows('Series UID', uid)

        return res
