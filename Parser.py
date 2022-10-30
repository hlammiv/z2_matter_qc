import csv
from itertools import count
import json
import ast
import numpy as np

measurements = []
new_measurements = []
nt = 20

with open('D:\z2_matter_qc\simulation_main\simulation_production_run_on_ibm_nairobi_data=2022-09-12_7qubits_ccfpd42pr43cku82ufkg.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        measurements.append(row['counts bare'])
#res = json.loads(measurement)
#counts = {'0x0': 23, '0x1': 13, '0x10': 16, '0x11': 25, '0x12': 17, '0x13': 15, '0x14': 18, '0x15': 24, '0x16': 14, '0x17': 10, '0x18': 14, '0x19': 15, '0x1a': 13, '0x1b': 15, '0x1c': 16, '0x1d': 19, '0x1e': 9, '0x1f': 13, '0x2': 11, '0x20': 12, '0x21': 15, '0x22': 10, '0x23': 14, '0x24': 19, '0x25': 16, '0x26': 7, '0x27': 14, '0x28': 19, '0x29': 26, '0x2a': 11, '0x2b': 12, '0x2c': 12, '0x2d': 14, '0x2e': 12, '0x2f': 8, '0x3': 14, '0x30': 24, '0x31': 21, '0x32': 16, '0x33': 14, '0x34': 11, '0x35': 22, '0x36': 18, '0x37': 11, '0x38': 21, '0x39': 22, '0x3a': 17, '0x3b': 17, '0x3c': 20, '0x3d': 18, '0x3e': 9, '0x3f': 14, '0x4': 16, '0x40': 22, '0x41': 20, '0x42': 12, '0x43': 15, '0x44': 19, '0x45': 31, '0x46': 15, '0x47': 11, '0x48': 28, '0x49': 28, '0x4a': 14, '0x4b': 10, '0x4c': 16, '0x4d': 11, '0x4e': 11, '0x4f': 12, '0x5': 17, '0x50': 10, '0x51': 11, '0x52': 16, '0x53': 13, '0x54': 27, '0x55': 13, '0x56': 11, '0x57': 14, '0x58': 22, '0x59': 19, '0x5a': 16, '0x5b': 15, '0x5c': 22, '0x5d': 24, '0x5e': 6, '0x5f': 15, '0x6': 11, '0x60': 17, '0x61': 16, '0x62': 18, '0x63': 16, '0x64': 20, '0x65': 17, '0x66': 10, '0x67': 8, '0x68': 25, '0x69': 16, '0x6a': 7, '0x6b': 15, '0x6c': 15, '0x6d': 24, '0x6e': 8, '0x6f': 6, '0x7': 16, '0x70': 18, '0x71': 14, '0x72': 9, '0x73': 14, '0x74': 19, '0x75': 15, '0x76': 16, '0x77': 10, '0x78': 14, '0x79': 19, '0x7a': 12, '0x7b': 8, '0x7c': 16, '0x7d': 16, '0x7e': 8, '0x7f': 10, '0x8': 17, '0x9': 21, '0xa': 12, '0xb': 12, '0xc': 20, '0xd': 27, '0xe': 20, '0xf': 16}
for measurement in measurements:
    counts_dict = ast.literal_eval(measurement)
    counts_dict2 = {}
    for key in counts_dict:
        counts_dict2[bin((int(key, 16)))[2:].zfill(7)] = counts_dict[key]
    new_measurements.append(counts_dict2)

counts_states = {}
counts_states_list = np.array([])
states = ['0000000', '0011100', '1101011']


for i in range(30):
    value = new_measurements[i]
    temp = {}
    for key in value:
        if key in states:
            temp[key] = new_measurements[i][key]


                
    
    counts_states_list = np.append(counts_states_list,temp)
    
    #counts_states_list = temp + counts_states_list
    
    




print(counts_states_list[28])