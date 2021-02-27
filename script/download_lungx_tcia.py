from tciaexplorer import TciaExplorer
import os
import json
import requests
import sys
from zipfile import ZipFile

dicom_path = 'DATA/DOI'

tcia = TciaExplorer()

patients = ["CT-Training-LC001", "CT-Training-LC002", "CT-Training-LC003"]
#patients = ["CT-Training-LC001", "CT-Training-LC002", "CT-Training-LC003",
#            "CT-Training-LC008", "CT-Training-LC009", "CT-Training-BE001",
#            "CT-Training-BE002", "CT-Training-BE006", "CT-Training-BE007",
#            "CT-Training-BE010",]

for patient in patients:
    ###################Fetch study###############
    print(patient)
    try:
        patientStudy  = json.loads(tcia.get_patient_study(patientID=patient).text)
    except requests.exceptions.RequestException as e:
        print(e)
        sys.exit(1)

    for study in patientStudy:
        ##########Fetch series#############
        print("Fetching series for the studyInstanceUID: "+ study['StudyInstanceUID'])
        try:
            patientSeries = json.loads(tcia.get_series(studyInstanceUID=study['StudyInstanceUID']).text)
        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)

        for series in patientSeries:
            print("Fetching images for the seriesInstanceUID: "+series['SeriesInstanceUID'])
            ##################Fetch images#############
            try:
                images = (tcia.get_image(seriesInstanceUID=series['SeriesInstanceUID']))
            except requests.exceptions.RequestException as e:
                print(e)
                sys.exit(1)

            #print(images)
            print("Writing image zipfile for the seriesInstanceUID: "+series['SeriesInstanceUID'])
            series_path = os.path.join(dicom_path,patient,study['StudyInstanceUID'],series['SeriesInstanceUID'])

            try:
                os.makedirs(series_path)
            except:
                pass

            fileName = os.path.join(series_path,patient+".zip")
            f = open(fileName, "wb")
            f.write(images.content)
            f.close()

            print("Unzipping images for the seriesInstanceUID: "+series['SeriesInstanceUID'])
            with ZipFile(fileName, 'r') as myzip:
                myzip.extractall(series_path)
            os.unlink(fileName)
