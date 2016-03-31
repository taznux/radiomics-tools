import csv
from collections import defaultdict

metadata = defaultdict(list)


def load_metadata(metadata_path):
    # load metadata
    f = open(metadata_path, 'r')
    nodule_info = csv.DictReader(f.read().splitlines(), dialect='excel')
    for n in nodule_info:
        for (k, v) in n.items():
            metadata[k].append(v)


def get(col, pid):
    res = ''
    try:
        idx = metadata['PID'].index(pid)
        res = metadata[col][idx]
    except ValueError as e:
        print(e)

    return res
