import pandas

metadata = pandas.DataFrame()

def load_metadata(metadata_path):
    global metadata
    # load metadata
    metadata = pandas.read_csv(metadata_path)
    try:
        metadata.set_index(['No.'], inplace=True)
    except KeyError:
        metadata.set_index(['No'], inplace=True)
    #print(metadata)


def getPatient(pid):
    res = ''
    #print(metadata)
    try:
        df = metadata.loc[metadata['PID'] == pid]
        res = df.to_dict('index')

        #print(res)

    except KeyError as e:
        print('KeyError: '+str(e))

    return res
