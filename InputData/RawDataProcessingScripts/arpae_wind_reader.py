from InputData.RawDataProcessingScripts.wind_data_reader import wind_data_reader

file_in = 'InputData/RawData/ARPAEWind.csv'
file_out = 'InputData/DataObjectInstances/arpae_wind.p'
Notes = 'Data are the US Department of Energy Methane Observation Networks with Innovative Technology to Obtain' +\
        'Reductions – MONITOR wind speed data, 2014.'
speed_col = 5
time_col = 2

wind_data_reader(file_in, speed_col=5, direction_col=None, time_col=2, notes=Notes, file_out=file_out)
