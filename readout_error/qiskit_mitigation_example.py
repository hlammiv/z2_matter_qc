import numpy as np
import matplotlib.pyplot as plt

# Importing standard Qiskit libraries
from qiskit import QuantumCircuit, QuantumRegister, transpile, assemble, Aer, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from ibm_quantum_widgets import *
from qiskit.providers.aer import QasmSimulator
from qiskit.ignis.mitigation.measurement import complete_meas_cal, CompleteMeasFitter
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
from qiskit.tools.monitor import job_monitor

# Loading your IBM Quantum account(s)
provider = IBMQ.load_account()

# Prepare a simulator to run calibration
aer_sim = Aer.get_backend('aer_simulator')

# Introduce a noise model
def get_noise(p):

    error_meas = pauli_error([('X',p), ('I', 1 - p)])

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_meas, "measure") # measurement error is applied to measurements
        
    return noise_model

noise_model = get_noise(0.3)


# Define a quantum register with 3 qubits and prepare the circuits for the readout calibration
# These two lines will set up all required circuits for the calibration
qr = QuantumRegister(3)
meas_calibs, state_labels = complete_meas_cal(qr=qr, circlabel='mcal')

# Run simulation with noise
t_qc = transpile(meas_calibs, aer_sim)
cal_results = aer_sim.run(t_qc, noise_model=noise_model, shots=10000).result()

# This object contains the calibration matrix
meas_fitter = CompleteMeasFitter(cal_results, state_labels, circlabel='mcal')

# Print the calibration matrix
print(meas_fitter.cal_matrix)


# Ok, let's prepare a circuit that prepares the |000> + |111> state and measures
# This is the bit where our Z2 simulation code comes in lateron
qc = QuantumCircuit(3,3)
qc.h(0)
qc.cx(0,1)
qc.cx(0,2)  
qc.measure([0, 1,2], [0, 1,2])

t_qc = transpile(qc, aer_sim)
results = aer_sim.run(t_qc, noise_model=noise_model, shots=10000).result()
noisy_counts = results.get_counts()

# Look at our results with noise
print(noisy_counts)

# In the following we use the calibration matrix and apply it
# Get the filter object
meas_filter = meas_fitter.filter

# Results with mitigation (apply basically multiplies on the matrix from above)
mitigated_results = meas_filter.apply(results)
mitigated_counts = mitigated_results.get_counts()

print(mitigated_counts)
