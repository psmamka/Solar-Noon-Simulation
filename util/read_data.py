# from pathlib import Path
import os
import csv
import numpy as np
import re

def read_noaa_data(fname, skip_lines=8, verbose=False, delimiter="\t"):
    data_ar = np.zeros((0,12))  # 12 columns, one for each month
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)

        for idx, row in enumerate(csv_reader):
            if idx < skip_lines:
                if verbose: print(row)
            else:
                # print(len(row))
                data_ar = np.append(data_ar, np.array(row[1:]).reshape(1,12), axis=0)
    
    pattern = '\d\d\:\d\d\:\d\d'    # hh:mm:ss format for noaa solar noon
    noon_ar = np.array([elm for elm in data_ar.transpose().ravel() if bool(re.search(pattern, elm))])
    return noon_ar

def read_navy_data(fname, skip_lines=9, verbose=False):
    data_ar = np.zeros((0,12))  # 12 columns, one for each month
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        for idx, row in enumerate(csv_reader):
            if idx < skip_lines or idx >= skip_lines + 31:
                if verbose: print(row)
            else:
                #                               '0725 1759'     2 spaces              '01'
                row_data = process_line(row[0], element_len=9, seperator_len=2, skip_chars=2, verbose=verbose)
                if verbose: print(row_data)
                data_ar = np.append(data_ar, np.array(row_data[:]).reshape(1,12), axis=0)
    
    pattern = '\d{4} \d{4}'    # hhmm hhmm format for navy sunrise/sunset
    srss_ar = np.array([elm for elm in data_ar.transpose().ravel() if bool(re.search(pattern, elm))])

    noon_ar = np.array([noon_from_sunrise_sunset(elm[0:2], elm[2:4], elm[5:7], elm[7:9]) for elm in srss_ar])
    return noon_ar


def process_line(str, element_len, seperator_len, skip_chars=0, verbose=False):
    line_ar = []
    idx = skip_chars
    step_len = seperator_len + element_len
    while idx + step_len <= len(str):
        line_ar.append(str[idx + seperator_len:idx + step_len])
        idx += step_len
    if verbose: print(line_ar)
    return line_ar

# convert sunrise sunset to noon data
def noon_from_sunrise_sunset(sr_hr, sr_mn, ss_hr, ss_mn):
    # average of the two time points
    noon_time = (int(sr_hr) * 3600 + int(sr_mn) * 60 + int(ss_hr) * 3600 + int(ss_mn) * 60) / 2.0

    noon_hr, mn_rem = np.divmod(noon_time, 3600)
    noon_mn, sc_rem = np.divmod(mn_rem, 60)
    noon_sc = round(sc_rem)
    return f"{int(noon_hr):02}:{int(noon_mn):02}:{int(noon_sc):02}"     # pad with zeros to length 2


if __name__ == "__main__":

    cur_path = os.path.dirname(__file__)
    fname1 = os.path.join(cur_path, '..', 'data', 'noaa-phoenix-solar-noon.txt')
    fname2 = os.path.join(cur_path, '..', 'data', 'navy-phoenix-sunrise-sunset.txt')
    print(fname1, fname2)

    data1 = read_noaa_data(fname1)
    data2 = read_navy_data(fname2, verbose=False)
    print(data1, data1.shape)
    print(data2, data2.shape)

