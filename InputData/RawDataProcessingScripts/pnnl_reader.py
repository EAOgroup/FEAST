import csv
import numpy as np
from InputData.input_data_classes import PNNLData
import pickle


def pnnl_data_reader(file_in=None, notes=None, file_out=None):
    """
    reads a csv file containing pnnl methane data.
    file_in:      Path to a raw data file
    notes:        Comments on the data
    file_out:     Path to store the data at
    """

    file = open(file_in)
    data = csv.reader(file, delimiter=',')
    n_lines = sum(1 for row in data) - 1
    file.seek(0)
    index = 0
    km, nu = np.zeros(n_lines), np.zeros(n_lines)
    for row in data:
        if row[0][0].isdigit():
            km[index] = float(row[1])
            nu[index] = float(row[0])
            index += 1
    if notes is None:
        notes = "Data from the PNNL spectral library"
    pnnl = PNNLData(notes, 'pnnl_methane.csv')
    pnnl.define_data(km, nu)
    pickle.dump(pnnl, open(file_out, 'wb'))
