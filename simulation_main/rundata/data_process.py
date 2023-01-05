#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 08:44:49 2022

@author: egustafs
"""

import pandas 
import numpy as np
import gvar as gv
import matplotlib.pyplot as plt


import matplotlib as mpl
from  matplotlib import rc
import matplotlib.pyplot as plt


colorsrgb = [(0, 0, 0), (230, 159, 0), (86, 180, 233),
          (0, 158, 115), (240, 228, 66), (0, 114, 178),
          (213, 94, 0), (204, 121, 167)]

fmts = ['o', 'd', 's', '^', '<', '8', '6', 'd', 'x']
colors = [tuple([el / 255 for el in c]) for c in colorsrgb]

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = "Times New Roman"
mpl.rcParams.keys()
# mpl.rcParams["text.fontset"] = "Times New Roman"
mpl.rcParams["mathtext.fontset"] = "stix"



class DataAnalysis(object):
    
    def __init__(self):
        d3_1 = self.load_data_set(3, 1)
        d3_2 = self.load_data_set(3, 2)
        d7_1 = self.load_data_set(7, 1)
        d7_2 = self.load_data_set(7, 2)
        self.data = {(7, 1): d7_1, (3, 1): d3_1, (7, 2): d7_2, (3, 2): d3_2}
        # counts_list31 = construct_readout_counts_dictionary(d3_1, 3)
        # counts_list32 = construct_readout_counts_dictionary(d3_1, 3)
        # counts_list71 = construct_readout_counts_dictionary(d7_2, 7)
        # counts_list72 = construct_readout_counts_dictionary(d7_2, 7)
        # counts1 = {(3, 1): counts_list31[0], (3, 2): counts_list32[0],
        #            (7, 1): counts_list71[0], (7, 2): counts_list72[0]}
        # counts2 = {(3, 1): counts_list31[1], (3, 2): counts_list32[1],
        #            (7, 1): counts_list71[1], (7, 2): counts_list72[1]}
        # cal_mats = construct_readout_confusion_matrix([counts1, counts2])


    def load_data_set(self, nqubits: int, paramnum: int):
        '''
        
    
        Parameters
        ----------
        nqubits : int
            whether to load the 3 qubit or 7 qubit data sets
        paramnum : int
            whether to load the param 1 or param 7 data sets
    
        Returns a list of the data files
        [no twirling or dynamic decoupling,
         twirling no dynamic decoupling,
         dynamic decoupling + twirling,
         dynamic decoupling + no twirling]
        -------
    
        '''
        
        if nqubits != 3 and nqubits != 7:
            print(f'nqubits was {nqubits}. Nqubits must be 3 or 7')
        elif nqubits == 7 and paramnum == 2:
            data0 = pandas.read_csv('NormanData/no_twirl_simulation_production_run_on_ibm_nairobi_data=2023-01-03_7qubits_ceq41lju4l1fkdjq3dlg.csv')
            data1 = pandas.read_csv('NormanData/FinalData/simulation_production_run_on_ibm_nairobi_data=2022-09-12_7qubits_ccfpd42pr43cku82ufkg.csv')
            data2 = pandas.read_csv('NormanData/FinalData/simulation_production_run_on_ibm_nairobi_data=2022-09-14_7qubits_6320fab55ccea78b019435f8.csv')
            data3 = [] #pandas.read_csv('NormanData/no_twirl_dd_simulation_production_run_on_ibm_nairobi_data=2023-01-03_3qubits_ceq42870c79flplhsis0.csv')
        elif nqubits == 3 and paramnum == 1:
            data0 = pandas.read_csv('NormanData/no_twirl_simulation_production_run_on_ibm_nairobi_data=2023-01-03_3qubits_ceq42870c79flplhsis0.csv')
            data1 = pandas.read_csv('NormanData/FinalData/simulation_production_run_on_ibm_nairobi_data=2022-09-10_3qubits_cbqi25b4fms1nkdes8fg.csv')
            data2 = pandas.read_csv('NormanData/FinalData/[DD]simulation_production_run_on_ibm_nairobi_data=2022-09-10_3qubits_63067c74e9f73b717085437c.csv')
            data3 = pandas.read_csv('NormanData/no_twirl_dd_simulation_production_run_on_ibm_nairobi_data=2023-01-04_3qubits_63b4a7bbe9c45ab44a611871.csv')
        elif nqubits == 7 and paramnum == 1:
            data0 = pandas.read_csv('ElizabethData/notwirl_simulation_production_run_on_ibmq_jakarta_data=2023-01-03_7qubits_ceq5r5u3ckj0h3h00te0.csv')
            data1 = pandas.read_csv('ElizabethData/7qubit_slapdashrun_main_ibmq_jakarta.csv')
            data2 = pandas.read_csv('ElizabethData/7qubit_slapdashrun_with_dynamic_decouple_ibmq_jakarta.csv')
            data3 = pandas.read_csv('ElizabethData/no_twirl_dd_simulation_production_run_on_ibmq_jakarta_data=2023-01-04_7qubits_63b58995ccb3671fe6d75d2a.csv')
        elif nqubits == 3 and paramnum == 2:
            data0 = pandas.read_csv('SarahData/notwirl_simulation_production_run_on_ibmq_jakarta_data=2023-01-03_3qubits_ceq5tnbu4l1fkdjq6esg.csv')
            data1 = pandas.read_csv('SarahData/simulation_production_run_on_ibmq_jakarta_data=2022-08-27_3qubits_cc4k6pr9k0hjph7mmrd0.csv')
            data2 = pandas.read_csv('SarahData/simulation_production_run_on_ibmq_jakarta_data=2022-09-22_3qubits_632ccd160cfeefcb807de983.csv')
            data3 = pandas.read_csv('SarahData/no_twirl_dd_simulation_production_run_on_ibmq_jakarta_data=2023-01-04_3qubits_63b58edd38d550e80d676096.csv')
        return [data0, data1, data2, data3]
        
    
    
    def construct_readout_confusion_matrix(self, counts_list : list):
        '''
        

        Parameters
        ----------
        counts_list : list
            list containing the counts to build the confusion matrix

        Returns
        -------
        cal_mats : dictionary of numpy arrays
            dictionary containing the confusion matrices

        '''
        cal_mats = {}
        counts1, counts2 = counts_list
        for key0 in counts1.keys():
            matrices = np.zeros((key0[0], 2, 2))
            for key in counts1[key0].keys():
                norm = sum(counts1[key0].values())
                # iterate through the qubits
                for j in range(key0[0]):
                    if key[j] == '1':
                        matrices[j, 0, 1] += counts1[key0][key] / norm
                    else:
                        matrices[j, 0, 0] += counts1[key0][key] / norm
            for key in counts2[key0].keys():
                norm = sum(counts2[key0].values())
                # iterate through the qubits
                for j in range(key0[0]):
                    if key[j] == '1':
                        matrices[j, 1, 1] += counts2[key0][key] / norm
                    else:
                        matrices[j, 1, 0] += counts2[key0][key] / norm
            invmats = np.array([np.linalg.inv(matrices[j]) 
                                for j in range(key0[0])])
        
            calibration_matrix = np.identity(1)
            for matrix in invmats:
                calibration_matrix = np.kron(calibration_matrix, matrix)
            cal_mats[key0] = calibration_matrix
        return cal_mats


    def get_no_twirl(self, key):
        data = self.data[key][0]
        nsteps = len(data['nt']) // 2
        obs_raw = np.zeros(nsteps)
        obs_readout = np.zeros(nsteps)
        for i in range(20):
            bare_counts = eval(data['counts bare'][2 * i])
            readout_counts = eval(data['counts mitigated'][2 * i])
            bc = {}
            rc = {}
            for key in bare_counts.keys():
                key2 = np.binary_repr(int(key, 16), width=3)
                bc[key2] = bare_counts[key]
            for key in readout_counts.keys():
                key2 = np.binary_repr(int(key, 16), width=3)
                rc[key2] = readout_counts[key]
            normalization = sum(bare_counts.values())
            for key in bc.keys():
                if key.count('1') % 2 == 0:
                    obs_raw[i] += bc[key] / normalization
                else:
                    obs_raw[i] -= bc[key] / normalization
            for key in rc.keys():
                if key.count('1') % 2 == 0:
                    obs_readout[i] += rc[key]
                else:
                    obs_readout[i] -= rc[key]
        return obs_raw, obs_readout


    def get_readout_shift(self):
        
        
        readout_shift = {}
        for key in self.data:
            obs_raw, obs_readout = self.get_no_twirl(key)
            readout_shift[key] = np.abs((obs_readout - obs_raw) / obs_raw)
        return readout_shift
            
            
        
        
        
        
        
if __name__ == "__main__":
    analyzer = DataAnalysis()
    
    readout_errors = True
    if readout_errors:
        rshifts = analyzer.get_readout_shift()
        counter = 0
        plt.figure(figsize=(3.375, 3.375))
        fill_style = ['none', 'full', 'full', 'none']
        for key in rshifts.keys():
            plt.plot(np.linspace(1, 20, 20),
                     rshifts[key], marker=fmts[counter], color=colors[counter],
                     label=r'$N_{q}$ = '+f'{key[0]}' 
                     + r' $m_0 = $' + f'{key[1]}',
                     fillstyle=fill_style[counter])
            counter += 1
        plt.ylabel(r'$R_{shift}$')
        plt.xlabel(r'$N_t$')
        plt.legend(ncol=1, fontsize='small', framealpha=0)
        plt.tight_layout()
        plt.savefig('readoutshift.pdf')






