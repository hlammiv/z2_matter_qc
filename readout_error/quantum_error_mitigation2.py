#################################################################################################################
#                           Authors: Sara Starecheski    &   Clement Charles                                    #
#                                      -----Read Out Error Mitigation-----                                      #
#                                                                                                               #
#       Description: This python code creates the quibit states |00>, |01>, |10>, |11>, and counts the          # 
#       states for a specified number shots from a quantum computer. Using this data a calibration matrix is    #
#       constructed. This matrix is manipulated such that the diagonal element are 1's.                         #
#                                                                                                               #
#                                               Dated: 3 Dec 2021                                               #
#################################################################################################################


#import modules here
import time
from qiskit import IBMQ, Aer, assemble, transpile, QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.visualization import array_to_latex
from qiskit.tools.monitor import job_monitor 
import numpy as np
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
import matplotlib.pyplot as plt
import pytz
from pytz import timezone
from datetime import datetime
from csv import writer
from scipy.stats import norm
from collections import defaultdict
import json


###---------Global Variables----------#
times = []
matrix_array_for_noise = []
normalized_matrix_array_for_noise = []
calibration_det_for_noise = []
datet = []

#Simulation noise, (not required when backends in use)
def get_noise(p):                                                       #Not used with quantum backends, only with simulation
    error_meas = pauli_error([('X',p), ('I', 1 - p)])
    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_meas, "measure")      # measurement error is applied to measurements
       
    return noise_model

noise_model = [0.01] #Not used


#IBMQ.save_account("d494816e01a737eb4ea87d4c9f6e2dcfdc7ed10d72ad68be30ecc7b5f66db445359005a848fc923f7ee96ba8909bffe76351ae2780cf53fd6a5f8d0ecc4fb491")
IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-education', group='fermilab-1', project='qjs-for-hep')
backend = provider.get_backend('ibm_perth')
#print("backend: ", backend)

shots = 5000

def time_of_interval(value,unit="s"):  #Function to create time interval depending on unit, set as seconds by default 
    if unit == "s":
        return value
    if unit == "m" or "min":
        return value * 60
    if unit == "h" or "hr" or "hour":
        return value * 60 * 60
 
    
t_i = time_of_interval(1)   #time in specififed unit of each interval where the function will execute

num_of_intervals = 75        #number of interval to be loop in t_i

    

def calibration(noise):                                  #This function creates calibration matrix for a set of state vectors
    aer_sim = Aer.get_backend('aer_simulator')
    calibration_matrix = []
    for state in ['00','01','10','11']:                  #Four qubit states
                results = []
                qc = QuantumCircuit(2,2)
                if state[0]=='1':
                    qc.x(1)
                if state[1]=='1':
                    qc.x(0)
                qc.measure([0, 1],[0,1])
                
                #copied from summer school code
                transpiled_qc = transpile(qc, backend=backend)
                assembled_qc = assemble(transpiled_qc)
                job = backend.run(transpiled_qc, shots=shots)
                job_monitor(job, interval=2)
                results = job.result().get_counts()
                #results = aer_sim.run(assembled_qc, noise_model=get_noise(noise), shots=shots).result().get_counts()
                for state2 in ['00','01','10','11']:
                    check = state2 in results.keys()
                    if check == False:
                        results.update({state2 : 0})
                    all_results = []
                list00 = []
                all_results.append(results)
                for n in all_results:                               #Construct calibration matrix
                    for m in ['00','01','10','11']:
                        list00.append(n.get(m)/shots) 
                calibration_matrix.append(list00)    
    return(calibration_matrix)

for n in range(num_of_intervals):               
    times.append(n*t_i)  

#Run's code automatically on a stated time-interval
for noise_vals in noise_model:
    matrix_array = []
    for n in range(num_of_intervals):                              #Loop for the number of intervals dersired
        close_time = time.time() + t_i                             #Adds delay time
        run = True
        datet.append(time.localtime())                             #Appends times to a list
        while run == True:
            if time.time() > close_time:
                matrix_array.append(calibration(noise_vals))       #Runs function which return calibration matrix
                run = False
    matrix_array_for_noise.append(matrix_array)



for matrix_array2 in matrix_array_for_noise:                       #Loop runs once with no simulated noise
    #Loop over each matrix
    normalized_matrix_array2 = []
    for matrix in matrix_array2: 
        normalized_matrix = []
        for i in range(len(matrix)):
            normalized_matrix00 = []
            for j in range(len(matrix[i])):
                normalized_matrix00.append(matrix[i][j]/np.sqrt(matrix[i][i]*matrix[j][j]))
            normalized_matrix.append(normalized_matrix00)
        normalized_matrix_array2.append(normalized_matrix)
    normalized_matrix_array_for_noise.append(normalized_matrix_array2)


#Save determinants
"""
for normalized_matrix_array2 in normalized_matrix_array_for_noise:
    calibration_det = []
    for mat in normalized_matrix_array2:
        calibration_det.append(1-np.linalg.det(np.array(mat)))
    calibration_det_for_noise.append(calibration_det)
"""
#off_diagonals = []

"""Defining arrays for each matrix element"""
all_off_diagonals = defaultdict(list)                                   #Instantiate defaultdict object

for matrix2 in normalized_matrix_array_for_noise:                       #Loop runs once with no simulated noise
    #Loop over each matrix
    normalized_matrix_array3 = []
    for matrix in matrix2: 
        normalized_matrix = []
        for i in range(len(matrix)):
            normalized_matrix00 = []
            for j in range(len(matrix[i])):
                if i != j:
                    all_off_diagonals[str(i)+str(j)].append(matrix[i][j])
                    if matrix[i][j] == 0:
                        all_off_diagonals[str(i)+str(j)].append(0)

for key in all_off_diagonals:
    ###-----Normal Distribution-----###
    std = np.std(all_off_diagonals[key],ddof=1)
    mean = np.mean(all_off_diagonals[key])
    domain = np.linspace(np.min(all_off_diagonals[key]),np.max(all_off_diagonals[key]))
    #print(off_diagonals) 
    ###-----Histogram Plotting-----###
    plt.plot(domain, norm.pdf(domain,mean,std))
    n, bins, patches = plt.hist(x=all_off_diagonals[key], bins=25, color='red',
                                alpha=0.7, rwidth=0.85,density=True)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of $M_{i,j}$ elements where ij = '+key)
    maxfreq = n.max()   
    plt.show()

### Save Matrix Data in txt File ###

with open('Documents/data.txt', 'w') as file:
    file.write(json.dumps(all_off_diagonals)) #Record dictionary
    file.writelines(["\n","times: ", str(datet)]) #Record times









### Code to Save Data to CSV ###


"""
with open("qc_data_histogram.csv", "a") as f_object:
    writer_object = writer(f_object)
    #writer_object.writerow(calibration_det_for_noise[0])
    #writer_object.writerow(datet)
    #writer_object.writerow(normalized_matrix_array_for_noise[0])
    writer_object.writerow(off_diagonals)
    #writer_object.writerow(all_off_diagonals)"""
"""
print(calibration_det_for_noise)
print(times)
"""
"""
with open("qc_data.csv", "a") as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(calibration_det_for_noise[0])
    writer_object.writerow(datet)
    writer_object.writerow(normalized_matrix_array_for_noise[0])
    writer_object.writerow(matrix_array_for_noise[0])
    writer_object.writerow(times)
    f_object.close()

for dets in calibration_det_for_noise:
    i = i + 1
    plt.plot(datet,dets)
plt.title("Variation of Norm-Calibration Matrix Determinants Over Time")
plt.xlabel("Time (hours:minutes:seconds)")
plt.ylabel("Matrix determinant")
#plt.legend()
plt.show()

"""
