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
        print('loading data files...')
        d3_1 = self.loadDataSet(3, 1)
        d3_2 = self.loadDataSet(3, 2)
        d7_1 = self.loadDataSet(7, 1)
        d7_2 = self.loadDataSet(7, 2)
        self.data = {(7, 1): d7_1, (3, 1): d3_1, (7, 2): d7_2, (3, 2): d3_2}
        
        dlist = [{(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1200]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1200])},
                 {(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1201]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1201])}]
        
        # print('processing unmitigated data...')
        # self.gen_no_twirl()
        # print('processing twirled data...')
        # self.gen_twirl_data()
        
        
        
    '''
    ========================================================================
    Data loading functions
    ========================================================================
    '''
    
    def load_simulation_data(self):
        '''
        load in the data to begin analysis

        Returns
        -------
        None.

        '''
        d3_1 = self.loadDataSet(3, 1)
        d3_2 = self.loadDataSet(3, 2)
        d7_1 = self.loadDataSet(7, 1)
        d7_2 = self.loadDataSet(7, 2)
        self.data = {(7, 1): d7_1, (3, 1): d3_1, (7, 2): d7_2, (3, 2): d3_2}
        
            
    def loadDataSet(self, nqubits: int, paramnum: int):
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
        drc = '../../'#'../../'
        if nqubits != 3 and nqubits != 7:
            print(f'nqubits was {nqubits}. Nqubits must be 3 or 7')
        elif nqubits == 7 and paramnum == 2:
            data0 = pandas.read_csv(drc+'ClementData/no_twirl_simulation_production_run_on_ibm_nairobi_data=2023-01-03_7qubits_ceq41lju4l1fkdjq3dlg.csv')
            data1 = pandas.read_csv(drc+'ClementData/simulation_production_run_on_ibm_nairobi_data=2022-09-12_7qubits_ccfpd42pr43cku82ufkg.csv')
            data2 = pandas.read_csv(drc+'ClementData/simulation_production_run_on_ibm_nairobi_data=2022-09-14_7qubits_6320fab55ccea78b019435f8.csv')
            data3 = [] #pandas.read_csv('NormanData/no_twirl_dd_simulation_production_run_on_ibm_nairobi_data=2023-01-03_3qubits_ceq42870c79flplhsis0.csv')
            data4 = pandas.read_csv(drc+'ClementData/simulation_production_run_on_ibmq_jakarta_data=2023-01-04_7qubits_cer2slve2dgblp3nfd80.csv')
            ret_list = [data0, data1, data2, data3, data4]
        elif nqubits == 3 and paramnum == 1:
            data0 = pandas.read_csv(drc+'NormanData/FinalData/no_twirl_simulation_production_run_on_ibm_nairobi_data=2023-01-03_3qubits_ceq42870c79flplhsis0.csv')
            data1 = pandas.read_csv(drc+'NormanData/FinalData/simulation_production_run_on_ibm_nairobi_data=2022-09-10_3qubits_cbqi25b4fms1nkdes8fg.csv')
            data2 = pandas.read_csv(drc+'NormanData/FinalData/[DD]simulation_production_run_on_ibm_nairobi_data=2022-09-10_3qubits_63067c74e9f73b717085437c.csv')
            data3 = pandas.read_csv(drc+'NormanData/FinalData/no_twirl_dd_simulation_production_run_on_ibm_nairobi_data=2023-01-04_3qubits_63b4a7bbe9c45ab44a611871.csv')
            data4 = pandas.read_csv(drc+'NormanData/simulation_production_run_on_ibmq_jakarta_data=2023-01-04_3qubits_cer3g3oc7aa76laosflg.csv')
            ret_list = [data0, data1, data2, data3, data4]
        elif nqubits == 7 and paramnum == 1:
            data0 = pandas.read_csv(drc+'ElizabethData/notwirl_simulation_production_run_on_ibmq_jakarta_data=2023-01-03_7qubits_ceq5r5u3ckj0h3h00te0.csv')
            data1 = pandas.read_csv(drc+'ElizabethData/7qubit_slapdashrun_main_ibmq_jakarta.csv')
            data2 = pandas.read_csv(drc+'ElizabethData/7qubit_slapdashrun_with_dynamic_decouple_ibmq_jakarta.csv')
            data3 = pandas.read_csv(drc+'ElizabethData/no_twirl_dd_simulation_production_run_on_ibmq_jakarta_data=2023-01-04_7qubits_63b58995ccb3671fe6d75d2a.csv')
            data4 = pandas.read_csv(drc+'ElizabethData/simulation_production_run_on_ibm_nairobi_data=2023-01-05_7qubits_cer4e5gc7aa76lap2j7g.csv')
            ret_list = [data0, data1, data2, data3, data4]
        elif nqubits == 3 and paramnum == 2:
            data0 = pandas.read_csv(drc+'SarahData/notwirl_simulation_production_run_on_ibmq_jakarta_data=2023-01-03_3qubits_ceq5tnbu4l1fkdjq6esg.csv')
            data1 = pandas.read_csv(drc+'SarahData/simulation_production_run_on_ibmq_jakarta_data=2022-08-27_3qubits_cc4k6pr9k0hjph7mmrd0.csv')
            data2 = pandas.read_csv(drc+'SarahData/simulation_production_run_on_ibmq_jakarta_data=2022-09-22_3qubits_632ccd160cfeefcb807de983.csv')
            data3 = pandas.read_csv(drc+'SarahData/no_twirl_dd_simulation_production_run_on_ibmq_jakarta_data=2023-01-04_3qubits_63b58edd38d550e80d676096.csv')
            data4 = pandas.read_csv(drc+'SarahData/simulation_production_run_on_ibm_nairobi_data=2023-01-04_3qubits_cer3v5v5otekiij83360.csv')
            ret_list = [data0, data1, data2, data3, data4]
        return ret_list
        

    def poissonError(self, value : float, norm : int):
        '''
        

        Parameters
        ----------
        value : float
            observable value
        norm : int
            number of samples

        Returns
        -------
        the gvar object with the correct value

        '''
        std = np.sqrt(value - value ** 2) / np.sqrt(norm)
        return gv.gvar(value, std)


    '''
    ========================================================================
    Bootstrapping dictionary 
    ========================================================================
    '''
    def bootstrapCountsDict(self, counts_dict: dict, nqubits: int,
                            nboot: int = 1000, samplesize: int=-1,
                            retvec=False, trim=False, normalize=False):
        boot_dicts = []
        # convert the counts dictionary to an array of probabilities
        vector = np.zeros(2 ** nqubits)
        for key in counts_dict.keys():
            if trim:
                keya = np.binary_repr(int(key, 16), width=nqubits)[:nqubits]
            else:
                keya = np.binary_repr(int(key, 16),
                                      width=2 *nqubits)[:nqubits]
            index = int(keya, 2)
            vector[index] += counts_dict[key]
        norm = sum(vector)
        if samplesize < 0:
            samplesize = int(norm)
        vector /= norm
        # if flagged to return the vector representation return vector rep
        if retvec:
            for i in range(nboot):
                boot_vec = np.random.choice(2 ** nqubits, size=samplesize,
                                            p=vector)
                boot_vec = np.array([len(np.where(boot_vec == j)[0])
                            for j in range(2 ** nqubits)], dtype='float64')
                if normalize:
                    boot_vec /= np.sum(boot_vec)
                boot_dicts.append(boot_vec)
            ret = np.array(boot_dicts)
            return ret
        # otherwise return the dictionary representation
        else:
            for i in range(nboot):
                boot_vec = np.random.choice(2 ** nqubits, size=samplesize,
                                            p=vector)
                boot_dict = {}
                for j in range(2 ** nqubits):
                    key = np.binary_repr(j, width=nqubits)
                    boot_dict[key] = len(np.where(boot_vec == j)[0])
                boot_dicts.append(boot_dict)
            return boot_dicts
            
    '''
    ========================================================================
    Readout Confusion Matrix Construction Functions
    ========================================================================
    '''
    
    def buildReadoutMatrixFromCountsHex(self, counts1: dict, counts2: dict,
                                        key0: tuple):
        '''
        
        builds the readout calibration matrices from a dictionary that has
        hexidecimal keys
        
        Parameters
        ----------
        counts1 : dict
            DESCRIPTION.
        counts2 : dict
            DESCRIPTION.
        key0 : tuple
            DESCRIPTION.

        Returns
        -------
        matrices : numpy array
            matrices containing the confusion matrices for each qubit
        norm : int
            sum of the number of counts

        '''
        matrices = np.zeros((key0[0], 2, 2))
        for key in counts1.keys():
            keya = np.binary_repr(int(key, 16),
                                  width=key0[0] * 2)[:key0[0]]
            norm = sum(counts1.values())
            # iterate through the qubits
            for j in range(key0[0]):
                if keya[j] == '1':
                    matrices[j, 0, 1] += counts1[key] / norm
                else:
                    matrices[j, 0, 0] += counts1[key] / norm
        for key in counts2.keys():
            norm = sum(counts2.values())
            # iterate through the qubits
            keya = np.binary_repr(int(key, 16),
                                  width=key0[0] * 2)[:key0[0]]
            for j in range(key0[0]):
                if keya[j] == '1':
                    matrices[j, 1, 1] += counts2[key] / norm
                else:
                    matrices[j, 1, 0] += counts2[key] / norm
        return matrices, norm
    
    
    def buildReadoutMatrixFromCountsBits(self, counts1: dict, counts2: dict,
                                           key0: tuple):
        '''
        
        builds the readout calibration matrices from a dictionary that has
        a vector
        
        Parameters
        ----------
        counts1 : dict
            DESCRIPTION.
        counts2 : dict
            DESCRIPTION.
        key0 : tuple
            DESCRIPTION.

        Returns
        -------
        matrices : numpy array
            matrices containing the confusion matrices for each qubit
        norm : int
            sum of the number of counts

        '''
        matrices = np.zeros((key0[0], 2, 2))
        for key in counts1.keys():
            keya = key
    
            norm = sum(counts1.values())
            # iterate through the qubits
            for j in range(key0[0]):
                if keya[j] == '1':
                    matrices[j, 0, 1] += counts1[key] / norm
                else:
                    matrices[j, 0, 0] += counts1[key] / norm
        for key in counts2.keys():
            norm = sum(counts2.values())
            # iterate through the qubits
            keya = key
            for j in range(key0[0]):
                if keya[j] == '1':
                    matrices[j, 1, 1] += counts2[key] / norm
                else:
                    matrices[j, 1, 0] += counts2[key] / norm
        return matrices, norm


    def constructReadoutConfusionMatrix(self, mat_type: str,
                                           counts_list: list, nboot = 1000):
        if mat_type == 'reg':
            return self.constructReadoutConfusionMatrixRegular(counts_list)
        if mat_type == 'gvar':
            return self.constructReadoutConfusionMatrixGvar(counts_list)
        if mat_type == 'bootstrap':
            return self.constructReadoutConfusionMatrixBootstrap(counts_list,
                                                                 nboot)
    
    
    def constructReadoutConfusionMatrixBootstrap(self, counts_list : list,
                                                 nboot : int):
        counts1, counts2 = counts_list
        cal_mats = {}
        for key in counts1.keys():
            cal_mats[key] = []
            c1_dict = self.bootstrapCountsDict(counts1[key], key[0], nboot)
            c2_dict = self.bootstrapCountsDict(counts2[key], key[0], nboot)
            for c1, c2 in zip(c1_dict, c2_dict):
                matrices, _ = self.buildReadoutMatrixFromCountsBits(c1, c2,
                                                                    key)
                invmats = np.array([np.linalg.inv(matrices[j]) 
                                    for j in range(key[0])])
            
                calibration_matrix = np.identity(1)
                for matrix in invmats:
                    calibration_matrix = np.kron(calibration_matrix, matrix)
                cal_mats[key].append(calibration_matrix)
        self.cal_mats = cal_mats


    def constructReadoutConfusionMatrixGvar(self, counts_list : list):
        
        
        # generate a dictionary to hold the calibration matrices
        cal_mats = {}
        # the counts lists for measuring all 0 and measuring all 1
        counts1, counts2 = counts_list
        # iterate through the keys in counts1
        for key0 in counts1.keys():
            matrices = np.zeros((key0[0], 2, 2), dtype=gv._gvarcore.GVar)
            # build first direction values
            vec1 = np.zeros((2 ** key0[0]))
            norm1 = sum(counts1[key0].values())
            for key in counts1[key0].keys():
                keya = np.binary_repr(int(key, 16),
                                      width=key0[0] * 2)[:key0[0]]
                index = int(keya, 2)
                vec1[index] += counts1[key0][key] / norm1
            sdev1 = np.diag(vec1)
            cov1 = sdev1 - np.outer(vec1, vec1)
            cov1 /= norm1
            vec_1 = gv.gvar(vec1, cov1)
            # build second direction vector
            vec1 = np.zeros((2 ** key0[0]))
            norm1 = sum(counts2[key0].values())
            for key in counts2[key0].keys():
                keya = np.binary_repr(int(key, 16),
                                      width=key0[0] * 2)[:key0[0]]
                index = int(keya, 2)
                print(key, index)
                vec1[index] += counts2[key0][key] / norm1
            sdev1 = np.diag(vec1)
            cov1 = sdev1 - np.outer(vec1, vec1)
            cov1 /= norm1
            vec_2 = gv.gvar(vec1, cov1)
            
            for j in range(len(vec_2)):
                key = np.binary_repr(j, width=key0[0])
                for k in range(key0[0]):
                    if key[k] == '1':
                        matrices[k, 0, 1] += vec_1[j]
                        matrices[k, 1, 1] += vec_2[j]
                    else:
                        matrices[k, 0, 0] += vec_1[j]
                        matrices[k, 1, 0] += vec_2[j]
            gvm_inv = np.array([gv.linalg.inv(matrices[j])
                                for j in range(matrices.shape[0])])
            calibration_matrices = np.identity(1, dtype=gv._gvarcore.GVar)
            for matrix in gvm_inv:
                calibration_matrices = np.kron(calibration_matrices, matrix)
            cal_mats[key0] = calibration_matrices
        self.cal_mats = cal_mats
            # print(gv.evalcov(matrices))
        #     matrices, norm = self.buildReadoutMatrixFromCounts(counts1[key0],
        #                                                         counts2[key0],
        #                                                         key0)
            
        #     std = np.sqrt(matrices - matrices ** 2) / np.sqrt(norm)
        #     gvm = [gv.gvar(matrices[j], std[j]) for j in range(key0[0])]
        #     gvm_inv = np.array([gv.linalg.inv(m) for m in gvm])
        #     calibration_matrix = np.identity(1, dtype=gv._gvarcore.GVar)
        #     for matrix in gvm_inv:
        #         calibration_matrix = np.kron(calibration_matrix, matrix)
        #     cal_mats[key0] = calibration_matrix
        # self.cal_mats = cal_mats
        
        
    # def constructReadoutConfusionMatrixGvar(self, counts_list : list):
        
        
    #     # generate a dictionary to hold the calibration matrices
    #     cal_mats = {}
    #     # the counts lists for measuring all 0 and measuring all 1
    #     counts1, counts2 = counts_list
    #     # iterate through the keys in counts1
    #     for key0 in counts1.keys():
    #         matrices, norm = self.buildReadoutMatrixFromCounts(counts1[key0],
    #                                                             counts2[key0],
    #                                                             key0)
            
    #         std = np.sqrt(matrices - matrices ** 2) / np.sqrt(norm)
    #         gvm = [gv.gvar(matrices[j], std[j]) for j in range(key0[0])]
    #         gvm_inv = np.array([gv.linalg.inv(m) for m in gvm])
    #         calibration_matrix = np.identity(1, dtype=gv._gvarcore.GVar)
    #         for matrix in gvm_inv:
    #             calibration_matrix = np.kron(calibration_matrix, matrix)
    #         cal_mats[key0] = calibration_matrix
    #     self.cal_mats = cal_mats

    
    def constructReadoutConfusionMatrixRegular(self, counts_list : list):
        '''
        Parameters
        ----------
        counts_list : list
            list containing the counts to build the confusion matrix

        Returns
        -------
        None
        '''
        cal_mats = {}
        counts1, counts2 = counts_list
        for key0 in counts1.keys():
            matrices, norm = self.buildReadoutMatrixFromCounts(counts1[key0],
                                                                counts2[key0],
                                                                key0)
            invmats = np.array([np.linalg.inv(matrices[j]) 
                                for j in range(key0[0])])
        
            calibration_matrix = np.identity(1)
            for matrix in invmats:
                calibration_matrix = np.kron(calibration_matrix, matrix)
            cal_mats[key0] = calibration_matrix
        self.cal_mats
     

    '''
    ========================================================================
    observable calculation functions
    ========================================================================
    '''
    
    def getTwirlDDGvar(self, keya: tuple):
        '''
        
        return the twirled circuits as gvars with and without readout
        error mitigation
        
        Parameters
        ----------
        keya : tuple
            key indicating which simulation run to select from
        
        as_gvar : bool
            whether to return the array as a gvar or as a single vector
        
        Returns
        -------
        
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
        
        '''
        data1 = self.data[keya][2]
        nsteps = len(data1['nt']) // 60
        
        vec_twirled_raw_1 = np.zeros((30, nsteps, 2 ** keya[0]),
                                     dtype=gv._gvarcore.GVar)
        vec_twirled_read_1 = np.zeros((30, nsteps, 2 ** keya[0]),
                                      dtype=gv._gvarcore.GVar)
        for i in range(nsteps):
            for j in range(30):
                print(i, j)
                bc = eval(data1['counts bare'][60 * i + j])
                vec = np.zeros((2 ** keya[0]))
                norm = sum(bc.values())
                for key in bc.keys():
                    index = int(key, 16)
                    vec[index] += bc[key] / norm
                sdev = np.diag(vec)
                cov = sdev - np.outer(vec, vec)
                cov /= norm
                vec = gv.gvar(vec, cov)
                vec_twirled_raw_1[j, i] += vec
                vec_twirled_read_1[j, i] += self.cal_mats[keya] @ vec
        return vec_twirled_read_1, vec_twirled_raw_1
    
    
    
    
    def getTwirlDDGvarFull(self, keya: tuple):
        '''
        
        return the twirled circuits as gvars with and without readout
        error mitigation
        
        Parameters
        ----------
        keya : tuple
            key indicating which simulation run to select from
        
        as_gvar : bool
            whether to return the array as a gvar or as a single vector
        
        Returns
        -------
        
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
        
        '''
        data1 = self.data[keya][2]
        nsteps = len(data1['nt']) // 60
        
        vec_twirled_raw_1 = np.zeros((2, 30, nsteps, 2 ** keya[0]),
                                     dtype=gv._gvarcore.GVar)
        vec_twirled_read_1 = np.zeros((2, 30, nsteps, 2 ** keya[0]),
                                      dtype=gv._gvarcore.GVar)
        for i in range(nsteps):
            for j in range(60):
                print(i, j)
                bc = eval(data1['counts bare'][60 * i + j])
                vec = np.zeros((2 ** keya[0]))
                norm = sum(bc.values())
                for key in bc.keys():
                    index = int(key, 16)
                    vec[index] += bc[key] / norm
                sdev = np.diag(vec)
                cov = sdev - np.outer(vec, vec)
                cov /= norm
                vec = gv.gvar(vec, cov)
                vec_twirled_raw_1[j // 30, j % 30, i] += vec
                cvec = self.cal_mats[keya] @ vec
                vec_twirled_read_1[j // 30, j % 30, i] += cvec
        return vec_twirled_read_1, vec_twirled_raw_1
    
    
    
    def getTwirlDDGvarOld(self, keya: tuple):
        '''
        
        return the twirled circuits as gvars with and without readout
        error mitigation
        
        Parameters
        ----------
        keya : tuple
            key indicating which simulation run to select from
        
        as_gvar : bool
            whether to return the array as a gvar or as a single vector
        
        Returns
        -------
        
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
            
        vec_twirled_raw_1 : numpy array
            vector of the twirling circuits counts without going through
            the calibration matrix
        
        '''
        data1 = self.data[keya][2]
        nsteps = len(data1['nt']) // 60
        
        vec_twirled_raw_1 = np.zeros((30, nsteps, 2 ** keya[0]),
                                     dtype=gv._gvarcore.GVar)
        vec_twirled_read_1 = np.zeros((30, nsteps, 2 ** keya[0]),
                                      dtype=gv._gvarcore.GVar)
        for i in range(nsteps):
            for j in range(30):
                bc = eval(data1['counts bare'][60 * i + j])
                vec = np.zeros((2 ** keya[0]), dtype=gv._gvarcore.GVar)
                norm = sum(bc.values())
                for key in bc.keys():
                    index = int(key, 16)
                    vec[index] += self.poissonError(bc[key] / norm, norm)
                vec_twirled_raw_1[j, i] += vec
                vec_twirled_read_1[j, i] += self.cal_mats[keya] @ vec
        return vec_twirled_read_1, vec_twirled_raw_1
        


    def getTwirlDDBootstrapForReadout(self, key: tuple, nboot: int=1000,
                                      nsize: int=-1):
        
        data1 = self.data[key][2]
        #calculate the number of trotter steps and number of twirling circuits
        nsteps = 20
        ntwirl = 30
        print(nsteps, ntwirl)
        oper = np.array([1, -1, -1, 1, -1, 1, 1, -1])
        if key[0] == 7:
            oper = np.kron(np.ones(4),
                           np.kron(oper, np.ones(4)))
        observables_readout = np.zeros((nboot, nsteps, ntwirl))
        for num in range(nsteps * ntwirl):
            i = num // ntwirl
            j = num % ntwirl
            bc = eval(data1['counts bare'][60 * i + j])
            boot_vecs = self.bootstrapCountsDict(bc, key[0], nboot=nboot, 
                                                 retvec=True, trim=True,
                                                 samplesize=nsize,
                                                 normalize=True)
            cal_mats = np.array(self.cal_mats[key])
            for k in range(nboot):
                vector = cal_mats[k] @ boot_vecs[k]
                obs = oper @ vector
                observables_readout[k, i, j] += obs
        
        return observables_readout
        
                    
  
    '''
    ========================================================================
    Analysis functions
    ========================================================================
    '''
    def run_gvar_analysis(self):
        dlist = [{(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1200]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1200])},
                 {(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1201]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1201])}]
        self.constructReadoutConfusionMatrixGvar(dlist)
        gvar_objects = {}
        for key in dlist[0].keys():
            vec, _ = self.getTwirlDDGvarFull(key)
            oper = np.array([1, -1, -1, 1, -1, 1, 1, -1])
            if key[0] == 7:
                oper = np.kron(np.ones(4),
                                np.kron(oper, np.ones(4)))
            gvar_objects[key] = vec @ oper
        gv.dump(gvar_objects, 'fullsetreadoutcorrelations.pkl')
        
    def run_readout_correlation_analysis_gvar(self):
        dlist = [{(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1200]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1200])},
                 {(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1201]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1201])}]
        self.constructReadoutConfusionMatrixGvar(dlist)
        gvar_objects = {}
        for key in dlist[0].keys():
            vec, _ = self.getTwirlDDGvar(key)
            oper = np.array([1, -1, -1, 1, -1, 1, 1, -1])
            if key[0] == 7:
                oper = np.kron(np.ones(4),
                                np.kron(oper, np.ones(4)))
            gvar_objects[key] = vec @ oper
        gv.dump(gvar_objects, 'readoutcorrelations.pkl')
            
    def run_readout_correlation_analysis_bootstrap(self, nboot: int=1000,
                                                   DD: bool=True):
        '''
        Parameters
        ----------
        nboot : int, optional
            number of bootstrap samples. The default is 1000.
        DD : bool, optional
            flag to use dynamic decoupling data. The default is True.

        Returns
        -------
        None.

        '''
        
        # construct a list containing the readout correction data
        dlist = [{(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1200]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1200]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1200])},
                 {(3, 1): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 1): eval(self.data[(7, 2)][2]['counts bare'][1201]),
                  (3, 2): eval(self.data[(3, 1)][2]['counts bare'][1201]),
                  (7, 2): eval(self.data[(7, 2)][2]['counts bare'][1201])}]
        # construct a list of the calibration matrices using nboot samples
        self.constructReadoutConfusionMatrix('bootstrap', dlist, nboot=nboot)
        # construct a list of the observables via bootstrapping
        # using the number of twirls and the number of time steps
        # obs_bootstrap = {}
        for key in dlist[0].keys():
            print('analyzing ' + str(key))
            self.obsred = self.getTwirlDDBootstrapForReadout(key, nboot=nboot)
            nq, param = key
            filename = f'bootstrapreadoutnq={nq}pnum={param}.npy'
            np.save(filename, self.obsred)
        
        
    def getCorrelationMatrix(self):
        pass


    '''
    ========================================================================
    backup functions
    ========================================================================
    '''


    def get_no_twirl(self, keya):
        
        data = self.data[keya][0]
        nsteps = len(data['nt']) // 2
        obs_raw = np.zeros(nsteps)
        obs_readout = np.zeros(nsteps)
        for i in range(20):
            bare_counts = eval(data['counts bare'][2 * i])
            readout_counts = eval(data['counts mitigated'][2 * i])
            bc = {}
            rc = {}
            obs_readout[i] += self.get_observable(readout_counts,
                                                  keya[0])
            obs_raw[i] += self.get_observable(bare_counts, keya[0])
        return obs_raw, obs_readout


    def gen_no_twirl(self):
        self.no_twirl_data = {}
        for key in self.data.keys():
            o1, o2 = self.get_no_twirl(key)
            self.no_twirl_data[key] = [o1, o2]


    def get_observable(self, counts, nqubits):
        count_dict = {}
        for key in counts.keys():
            key2 = np.binary_repr(int(key, 16), nqubits)
            count_dict[key2] = counts[key]
        ret_val = 0
        if nqubits == 3:
            for key in count_dict.keys():
                if key.count('1') % 2 == 0:
                    ret_val += count_dict[key]
                else:
                    ret_val -= count_dict[key]
        elif nqubits == 7:
            for key in count_dict.keys():
                if key[2:5].count('1') % 2 == 0:
                    ret_val += count_dict[key]
                else:
                    ret_val -= count_dict[key]
        ret_val /= sum(count_dict.values())
        return ret_val
    

    def get_twirled(self, keya):
        data1 = self.data[keya][1]
        data2 = self.data[keya][4]
        nsteps = 20
        print(nsteps)
        obs_twirled_raw_1 = np.zeros((30, nsteps))
        obs_twirled_readout_1 = np.zeros((30, nsteps))
        obs_twirled_raw_2 = np.zeros((30, nsteps))
        obs_twirled_readout_2 = np.zeros((30, nsteps))
        for i in range(nsteps):
            for j in range(30):
                bc = eval(data1['counts bare'][60 * i + j])
                rc = eval(data1['counts mitigated'][60 * i + j])
                obs_twirled_raw_1[j, i] += self.get_observable(bc, keya[0])
                obs_twirled_readout_1[j, i] += self.get_observable(rc,
                                                                   keya[0])
                bc = eval(data2['counts bare'][60 * i + j])
                rc = eval(data2['counts mitigated'][60 * i + j])
                obs_twirled_raw_2[j, i] += self.get_observable(bc, keya[0])
                obs_twirled_readout_2[j, i] += self.get_observable(rc,
                                                                   keya[0])
        return [obs_twirled_raw_1, obs_twirled_raw_2, obs_twirled_readout_1, 
                obs_twirled_readout_2]
            
    def gen_twirl_data(self):
        twirl_data = {}
        for key in self.data.keys():
            oraw1, oraw2, oread1, oread2 = self.get_twirled(key)
            twirl_data[key] = [oraw1, oraw2, oread1, oread2]
        self.twirl_data = twirl_data


    def get_readout_shift(self):
        readout_shift = {}
        for key in self.data.keys():
            print(key)
            oraw1, oraw2, oread1, oread2 = self.get_twirled(key)
            shift1 = np.abs(oread1 - oraw1)
            shift2 = np.abs(oread2 - oraw2)
            readout_shift[key] = [[shift1, oraw1], [shift2, oraw2]]
        return readout_shift
    

    # def get_readout_shift(self):
        
        
    #     readout_shift = {}
    #     for key in self.data:
    #         obs_raw, obs_readout = self.get_no_twirl(key)
    #         readout_shift[key] = np.abs((obs_readout - obs_raw) / obs_raw)
    #     return readout_shift
            
            
def make_readout_shift_versus_machine_plot(rshifts):
    fig, ax = plt.subplots(ncols=2, nrows=2)
    
    counter = 0
    xpts = np.abs(rshifts[(3, 1)][0][1].flatten())
    ax[0, 0].plot(np.abs(rshifts[(3, 1)][0][1].flatten()),
             rshifts[(3, 1)][0][0].flatten() / xpts,
             marker=fmts[counter], color=colors[counter],
             fillstyle='none',
             linestyle='none', markersize=3, label='ibm_nairobi')
    counter = 1
    # ax[0, 0].plot(np.abs(rshifts[(3, 1)][1][1].flatten()),
    #          rshifts[(3, 1)][1][0].flatten(),
    #          marker=fmts[counter], color=colors[counter],
    #          fillstyle='full',
    #          linestyle='none', markersize=3, label='ibm_jakarta')
    ax[0, 0].annotate(r'$N_s = 2$ $m_0 = 1$',(0.75, 0.15), 
                      xycoords='axes fraction',
                      horizontalalignment='center')
    
    counter = 0
    xpts = np.abs(rshifts[(7, 1)][1][1].flatten())
    ax[1, 0].plot(np.abs(rshifts[(7, 1)][1][1].flatten()),
             rshifts[(7, 1)][1][0].flatten() / xpts,
             marker=fmts[counter], color=colors[counter],
             fillstyle='none',
             linestyle='none', markersize=3, label='ibm_nairobi')
    counter = 1
    # ax[1, 0].plot(np.abs(rshifts[(7, 1)][0][1].flatten()),
    #          rshifts[(7, 1)][0][0].flatten(),
    #          marker=fmts[counter], color=colors[counter],
    #          fillstyle='full',
    #          linestyle='none', markersize=3, label='ibm_jakarta')
    ax[1, 0].annotate(r'$N_s = 4$ $m_0 = 1$',(0.75, 0.25), 
                      xycoords='axes fraction',
                      horizontalalignment='center')
    
    
    counter = 0
    xpts = np.abs(rshifts[(3, 2)][1][1].flatten())
    ax[0, 1].plot(np.abs(rshifts[(3, 2)][1][1].flatten()),
              rshifts[(3, 2)][1][0].flatten() / xpts,
              marker=fmts[counter], color=colors[counter],
              fillstyle='none',
              linestyle='none', markersize=3, label='ibm_nairobi')
    counter = 1
    # ax[0, 1].plot(np.abs(rshifts[(3, 2)][1][1].flatten()) +
    #               np.random.rand(600) * 0.05,
    #           rshifts[(3, 2)][1][0].flatten() + 
    #                         0.07 * np.abs(rshifts[(3, 2)][1][1].flatten()),
    #           marker=fmts[counter], color=colors[counter],
    #           fillstyle='full',
    #           linestyle='none', markersize=3, label='ibm_jakarta')
    ax[0, 1].yaxis.tick_right()
    ax[0, 1].annotate(r'$N_s = 2$ $m_0 = 2$',(0.75, 0.25), 
                      xycoords='axes fraction',
                      horizontalalignment='center')
    
    counter = 0
    xpts = np.abs(rshifts[(7, 2)][0][1].flatten())
    ax[1, 1].plot(np.abs(rshifts[(7, 2)][0][1].flatten()),
             rshifts[(7, 2)][0][0].flatten() / xpts,
             marker=fmts[counter], color=colors[counter],
             fillstyle='none',
             linestyle='none', markersize=3, label='ibm_nairobi')
    # counter = 1
    # ax[1, 1].plot(np.abs(rshifts[(7, 2)][1][1].flatten()),
    #          rshifts[(7, 2)][1][0].flatten(),
    #          marker=fmts[counter], color=colors[counter],
    #          fillstyle='full',
    #          linestyle='none', markersize=3, label='ibm_jakarta')
    ax[1, 1].yaxis.tick_right()
    # ax[1, 1].legend(framealpha=0)
    ax[1, 1].annotate(r'$N_s = 4$ $m_0 = 2$',(0.75, 0.25), 
                      xycoords='axes fraction',
                      horizontalalignment='center')
    for i in range(4):
        ax[i // 2, i % 2].set_ylim(0, 0.4)
    fig.text(0.05, 0.5, r'$\mathcal{R}$',
             horizontalalignment='center', verticalalignment='center',
             rotation=90)
    
    fig.text(0.5, 0.05, r'$\langle XZX \rangle_{unmitigated}$',
             horizontalalignment='center', verticalalignment='center')
    fig.savefig('readoutshiftsversusmachines.pdf')        
        
        
def make_readout_shift_versus_circuit_depth(rshifts):
    xpts = np.linspace(1, 20, 20)
    key = (3, 1)
    ypts1 = gv.gvar(np.mean(rshifts[key][0][0], axis=0),
             np.std(rshifts[key][0][0], axis=0))
    
    ypts1 /= gv.gvar(np.abs(np.mean(rshifts[key][0][1], axis=0)),
             np.std(rshifts[key][0][1], axis=0))
    
    key = (7, 1)
    ypts2 = gv.gvar(np.mean(rshifts[key][1][0], axis=0),
              np.std(rshifts[key][1][0], axis=0))
    
    ypts2 /= gv.gvar(np.abs(np.mean(rshifts[key][1][1], axis=0)),
              np.std(rshifts[key][1][1], axis=0))
    
    key = (3, 2)
    ypts3 = gv.gvar(np.mean(rshifts[key][1][0], axis=0),
             np.std(rshifts[key][1][0], axis=0))
    ypts3 /= gv.gvar(np.abs(np.mean(rshifts[key][1][1], axis=0)),
             np.std(rshifts[key][1][1], axis=0))

    key = (7, 2)
    ypts4 = gv.gvar(np.mean(rshifts[key][0][0], axis=0),
              np.std(rshifts[key][0][0], axis=0))
    ypts4 /= gv.gvar(np.abs(np.mean(rshifts[key][0][1], axis=0)),
              np.std(rshifts[key][0][1], axis=0))
    counter = 0
    plt.figure()
    plt.errorbar(xpts, gv.mean(ypts1), yerr=gv.sdev(ypts1), fmt=fmts[counter],
                 color=colors[counter], label=r'$N_s = 2$ $m_0 = 1$',
                 capsize=5)
    counter += 1
    # plt.errorbar(xpts, gv.mean(ypts2), yerr=gv.sdev(ypts2), fmt=fmts[counter],
    #               color=colors[counter], label=r'$N_s = 4$ $m_0 = 1$',
    #               capsize=7)
    counter += 1
    plt.errorbar(xpts, gv.mean(ypts3), yerr=gv.sdev(ypts3), fmt=fmts[counter],
                  color=colors[counter], label=r'$N_s = 2$ $m_0 = 2$',
                  capsize=9)
    counter += 1
    # plt.errorbar(xpts, gv.mean(ypts4), yerr=gv.sdev(ypts4), fmt=fmts[counter],
    #               color=colors[counter], label=r'$N_s = 4$ $m_0 = 2$',
    #               capsize=11)
    plt.legend(framealpha=0, ncol=2)
    # plt.yscale('log')
    plt.ylim(0.01, 1)
    plt.xlabel(r'$N_t$')
    plt.ylabel(r'Relative Shift')
    plt.savefig('relativereadoutshiftversuscircuitdepth.pdf')
    
        
if __name__ == "__main__":
    analyzer = DataAnalysis()
    # analyzer.run_gvar_analysis()
    analyzer.run_readout_correlation_analysis_gvar()
    # a, b = analyzer.getTwirlDDGvar((3, 1))
    # analyzer.run_readout_correlation_analysis_bootstrap(nboot=2000)
    readout_errors = False
    if readout_errors:
        rshifts = analyzer.get_readout_shift()
        make_readout_shift_versus_machine_plot(rshifts)
        # make_readout_shift_versus_circuit_depth(rshifts)
        






