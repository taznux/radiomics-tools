import SimpleITK as sitk

def load_dicom(input_directory, **args):
    # Find the Dicom series
    reader = sitk.ImageSeriesReader()

    target_series = ""

    series_found = reader.GetGDCMSeriesIDs(input_directory)
    # Process each Dicom series
    if len(series_found):
        for serie in series_found:
            print("Series:", serie )

            # Get the Dicom filename corresponding to the current series
            dicom_names = reader.GetGDCMSeriesFileNames(input_directory, serie)
            #print( "\nFiles in series: ", dicom_names )

            if len(dicom_names):
                if (target_series == "" or target_series == serie):
                    reader.SetFileNames(dicom_names)
                    image = reader.Execute()
                    print("Image size: ", image.GetSize())
                    return image
    else:
        raise(AttributeError)
