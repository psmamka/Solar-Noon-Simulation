# from pathlib import Path
import os
import csv
import numpy as np
import re

def read_noaa_data(fname, skip_lines=8, verbose=False, delimiter="\t"):
    data_ar = np.zeros((0,12))
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)

        for idx, row in enumerate(csv_reader):
            if idx < skip_lines:
                if verbose: print(row)
            else:
                # print(len(row))
                data_ar = np.append(data_ar, np.array(row[1:]).reshape(1,12), axis=0)
    
    pattern = '\d\d\:\d\d\:\d\d'    # hh:mm:ss
    noon_ar = np.array([elm for elm in data_ar.transpose().ravel() if bool(re.search(pattern, elm))])
    return noon_ar


if __name__ == "__main__":
    # fname = "../data/noaa-phoenix-solar-noon.txt"

    cur_path = os.path.dirname(__file__)
    fname = os.path.join(cur_path, '..', 'data', 'noaa-phoenix-solar-noon.txt')
    print(fname)

    data = read_noaa_data(fname)
    print(data)
    print(data.shape)
