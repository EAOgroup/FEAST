from InputData.RawDataProcessingScripts.repair_cost_data_reader import repair_cost_data_reader

file_in = 'InputData/RawData/FernandezRepairCost.csv'
file_out = 'InputData/DataObjectInstances/fernandez_leak_repair_costs_2006.p'
Notes = 'Data are from Fernandez, R. "Cost Effective Directed Inspection and Maintenance..." National Gas Machinery' + \
        'Laboratory. 2006'

repair_cost_data_reader(file_in, Notes, file_out)
