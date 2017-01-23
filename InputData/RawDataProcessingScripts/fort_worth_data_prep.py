""" fortWorthDataPrep loads data from the Fort Worth emissions data spreadsheet and processes it before saving the
    results to a DataFile object. The Fort Worth emissions data are available from the following resource:
    City of Fort Worth Natural Gas Air Quality Study. Tech. rep. Fort Worth, TX:
    Eastern Research Group and Sage Environmental Consulting LP. for the City
    of Fort Worth, 2011. <http://fortworthtexas.gov/gaswells/air-quality-study/final> (Emissions Calculations Workbook)
    For convenience, the sheet titled '2.2 Emission Data' was exported from the original document, commas were stripped
    and the file is saved as a csv in the DataFiles directory of PyFEAST.
"""

# -------------- reading the csv file --------------
import csv
from InputData.input_data_classes import LeakData
import pickle

# -------------- Hard coded values --------------
# In some cases, a leak would be detected with an FID or IR camera, but no flux could be measured with the HI-FLOW
# sampler. In these cases, the study team assigned a flux of 0.001 cfm to the leak. These data are omitted from FEAST.
cfm_unmeasured_value = 0.001  # cfm
# Number of wells surveyed with an IR camera and FID in the Fort Worth study
n_wells_IR, n_wells_FID = 1138, 114
# Unit conversion from cfm to g/s (assuming standard conditions and pure methane)
cfm_to_gps = 0.0283/60*1e5/8.314/293*16

flux_IR, flux_FID = [], []  # g/s
with open('InputData/RawData/FortWorth.csv') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in data:
        if row[0].isdigit():
            flux = float(row[28])
            if flux != cfm_unmeasured_value:
                flux *= cfm_to_gps
                if row[29] == '--':
                    flux_FID.append(flux)
                else:
                    flux_IR.append(flux)

notes = \
    """Data extracted from City of Fort Worth Natural Gas Air Quality Study. Tech. rep. Fort Worth, TX:
    Eastern Research Group and Sage Environmental Consulting LP. for the City
    of Fort Worth, 2011. Flux data are recorded in grams/second."""

fort_worth_leaks = LeakData(notes=notes, raw_file_name='FortWorth.csv', data_prep_file='fortWorthDataPrep')

leak_data = {'IR': flux_IR, 'FID': flux_FID}
well_counts = {'IR': n_wells_IR, 'FID': n_wells_FID}
fort_worth_leaks.define_data(leak_data=leak_data, well_counts=well_counts)

pickle.dump(fort_worth_leaks, open('InputData/DataObjectInstances/fort_worth_leaks.p', 'wb'))
