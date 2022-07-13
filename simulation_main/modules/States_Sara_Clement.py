#import
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

def one_meson_state_prep(nqubits):                #function that creates meson at the center of N qubits
    qc = QuantumCircuit(nqubits, nqubits)         #Create Quantum Circuit with N qubits
    pos = nqubits//2
    qc.x(pos-1)                                   
    qc.x(pos)
    qc.x(pos+1)
    for i in range(1, nqubits, 2):
        qc.h(i)
    #print(qc)
    return qc

def all_zero_state_prep(nqubits):                 #all zero state for N qubits
    qc = QuantumCircuit(nqubits, nqubits)
    for i in range(1, nqubits, 2):
        qc.h(i)
    #print(qc)
    return qc

def gauge_transformation(qc, pos=None):                    #function applies gauge transformation XZX
    pos=(qc.num_qubits//2)+1
    #print(pos)
    originalstate = Statevector.from_instruction(qc)  #Record state vector of state before gauge transformation           
    qc.x(pos-1)
    qc.z(pos)
    qc.x(pos+1)
    transformedstate = Statevector.from_instruction(qc) #Record transformed state
    if originalstate == transformedstate:               #check if gauge transformation altered state vector
        gaugeinvariant = True
    else: 
        gaugeinvariant = False
    print("Gauge Invariance : ", gaugeinvariant)
    #print(qc)
    return qc