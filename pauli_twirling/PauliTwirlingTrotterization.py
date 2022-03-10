from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit, execute
from qiskit import Aer, IBMQ
from qiskit import tools
from qiskit.test import mock 
import qiskit
from qiskit.quantum_info.operators import Operator
from qiskit.compiler import transpile
from qiskit.visualization import plot_gate_map
from qiskit.providers.aer import noise
from qiskit.visualization import array_to_latex
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor
from qiskit.quantum_info import Operator 
import qiskit.quantum_info as qi

import numpy
import qiskit

TWIRL_CNOT_DICT = {(0, 0): (0, 0), (0, 1): (0, 1),
    (0, 2): (3, 2), (0, 3): (3, 3), (1, 0): (1, 1),
    (1, 1): (1, 0), (1, 2): (2, 3), (1, 3): (2, 2),
    (2, 0): (2, 1), (2, 1): (2, 0), (2, 2): (1, 3),
    (2, 3): (1, 2), (3, 0): (3, 0), (3, 1): (3, 1),
    (3, 2): (0, 2), (3, 3): (0, 3)}

def generate_twirl(quantum_circuit, nqubits):
    ''' Applies Random Pauli Gates to the circuit
    
    Parameters:
    quantum_circuit (qiskit.QuantumCircuit): quantum circuit to twirl_trotter_
    nqubits (int): the number of qubits 
    '''
def twirl_trotter_two_site():

def twirl_trotter_ 