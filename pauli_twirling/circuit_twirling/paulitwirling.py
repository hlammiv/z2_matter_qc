"""from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit, execute
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
import qiskit"""

import numpy
from . import twirlingconstants

def twirl_cnot(quantum_circuit, target_qubit, control_qubit):
    """
    Applies a random pair of Paulis before and complementary Pauli Gates
    after a CNOT that is equivalent to a CNOT upto a global phase
    
    
    ----|A_c|--- . ---|B_C|----
                 |  
    ----|A_t|--- + ---|B_t|----
    
    
    Parameters
    
    ----------
    
    quantum_circuit: qiskit.QuantumCircuit
        The quantum circuit that will have the twirled CNOT
    target_qubit : int
        The target qubit for the CNOT gate
    control_qubit : int
        The control qubit for the CNOT gate
      
    
    Returns
    
    -------
    
    quantum_circuit : qiskit.QuantumCircuit
        the quantum circuit with the twirled CNOT applied
    
    """
    
    # use numpy to generate an array of two random numbers
    
    paulis = None
    
    # this is a list which will apply the quantum gates
    gates = [lambda x: 0, quantum_circuit.x, quantum_circuit.y,
             quantum_circuit.z]

    """
    for each element in paulis
        use the list gates to apply the pauli gate onto 
        the respective qubit
    """
    
    # apply the cnot gate on control_qubit and target_qubit
    pass
    
    """
    use twirling_constants.TWIRL_CNOT_DICT to extract the complementary gate. You will need to turn the list into a 
    tuple via (A_c, A_t) the returned pair will be called
    comp
    """
    pauli_key = None
    
    # the complementary gate index (like paulis) returned from dict
    comp = None
    
    """
    for each element in comp
        use the list gates to apply the pauli gate onto 
        the respective qubit
    """
    pass
    
    # return the quantum circuit
    return quantum_circuit


def twirl_hard_cycle(quantum_circuit, num_qubits, tc_pairs):
    """
    Apply twirls across a CNOT hard cycle (simultaneously applied CNOTS) 
    including spectator qubits
    
    --- |A_c| -- . -- |B_c| ---
                 |             
    --- |A_t| -- + -- |B_t| --- 
    
    ----|A_1| ------- |A_1+| --
    
    --- |C_c| -- . -- |D_c| ---
                 |
    --- |C_t| -- x -- |D_t| ---
    
    
    Parameters
    
    ----------
    
    quantum_circuit : qiskit.QuantumCircuit
        the quantum circuit to apply the cnot hard cycle Pauli Twirl tobytes
    num_qubits : int
        the number of qubits in the quantum register
    tc_pairs : list of tuples
        a list containing tuples whose first element is the control qubit
        and second element is the target qubit
    
    Returns
    
    -------
    
    quantum_circuit : qiskit.QuantumCircuit
        The quantum circuit 
    
    """

    # this is a list which will apply the quantum gates
    gates = [lambda x: 0, quantum_circuit.x, quantum_circuit.y,
             quantum_circuit.z]
    
    # generate a list of the spectator qubits
    # this should be all qubits not in tc_pairs
    
    spectators = None
    
    """
    for each qubit in spectators:
        generate a random integer between 0 and 3 (inclusive)
        apply the pauli gate on the spectator 
        apply a barrier
        apply the same pauli gate on the spectator
    """
    
    pass
    
    """
    for each pair in tc_pairs:
        apply each pauli gate for the tc_pair
    """
    
    pass 
    
    # return the quantum circuit
    
    return quantum_circuit


