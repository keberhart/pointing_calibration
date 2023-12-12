#!/usr/bin/env python3

import csv
import datetime as D
import time

def get_csv_data(file_path):
    ''' Open a CSV file, return as a list. Full path expected.
    '''
    data_list = []
    with open(file_path, newline='') as data_file:
        data_set = csv.reader(data_file, delimiter=',')
        for row in data_set:
            data_list.append(row)

    return data_list

def peak_sort(center_point, span, data_set):
    ''' Find the point with the highest power.

        center_point is the expected peak time from ephemeris, python datetime obj
        span is the duration of the dwell in minutes
    '''
    max_bin = None
    first_time = center_point - (D.timedelta(minutes=span) / 2)
    last_time = center_point + (D.timedelta(minutes=span) / 2)
    date = first_time.date()
    for item in data_set:
        line_time = D.datetime(*(time.strptime(item[0], "%H:%M:%S.%f")[0:6]))
        line_time.replace(tzinfo=D.UTC)
        line_dtg = D.datetime.combine(first_time.date(), line_time.time(), tzinfo=D.UTC)
        if line_dtg < first_time:
            continue
        if line_dtg > last_time:
            break
        if max_bin == None:
            new_row = item
            new_row.append(line_dtg)
            max_bin = new_row
        if float(max_bin[1]) < float(item[1]):
            new_row = item
            new_row.append(line_dtg)
            max_bin = new_row
    return max_bin

def test_times(time_set, span, data_set):
    ''' Test the data for peaks surrounding the center times.
    '''
    peak_values = []
    for item in time_set:
        center_time = D.datetime(*(time.strptime(item[0], "%Y:%m:%d:%H:%M:%S")[0:6]),tzinfo=D.UTC)
        max_bin = peak_sort(center_time, span, data_set)
        if max_bin:
            max_bin.extend(item[1:3])
            max_bin.extend([center_time])
            peak_values.append(max_bin)
    return peak_values


if __name__ == "__main__":
    DATA_FILE_PATH = 'data/Band_Power1.csv'
    TIME_FILE_PATH = 'data/26M_CAS_A_samples.txt'
    data_set = get_csv_data(DATA_FILE_PATH)
    time_set = get_csv_data(TIME_FILE_PATH)
    results = test_times(time_set, 15, data_set)
    # now we have our test data lined up. We know the time when CAS_A was in
    # the beam and the angles we were pointed at. So we've got to figure out;
    # What is the difference in angles between actual CAS_A and the antenna
    # readouts.
    for item in results:
        print(item)
#    test_time = D.datetime(year=2023, month=11, day=2, hour=4, minute=15, second=48)
#    max_bin = peak_sort(test_time, 15, data_set)
#    print(max_bin)
