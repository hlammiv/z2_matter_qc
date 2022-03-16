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

import random
import numpy
from . import twirlingconstants

def twirl_cnot(quantum_circuit, control_qubit, target_qubit):
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
    
    paulis =[random.randint(0,3),random.randint(0,3)]

    # this is a list which will apply the quantum gates
    gates = [lambda x: 0, quantum_circuit.x, quantum_circuit.y,
             quantum_circuit.z]

    """
    for each element in paulis
        use the list gates to apply the pauli gate onto 
        the respective qubit
    """

    # apply the cnot gate on control_qubit and target_qubit

    gates[paulis[0]](control_qubit)
    gates[paulis[1]](target_qubit)

    pass

    quantum_circuit.barrier()
    quantum_circuit.cx(control_qubit,target_qubit)
    quantum_circuit.barrier()
    
    """
    use twirlingconstants.TWIRL_CNOT_DICT to extract the complementary gate. You will need to turn the list into a 
    tuple via (A_c, A_t) the returned pair will be called
    comp
    """

    pauli_key = (paulis[0],paulis[1])

    # the complementary gate index (like paulis) returned from dict
    comp =  twirlingconstants.TWIRL_CNOT_DICT[pauli_key]
    
    """
    for each element in comp
        use the list gates to apply the pauli gate onto 
        the respective qubit
    """

    gates[comp[0]](control_qubit)
    gates[comp[1]](target_qubit)

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
    
    spectators = []
    cxqubits = []

    for pair in tc_pairs:
        for qubit in pair:
            cxqubits.append(qubit)

    for qubit in range(num_qubits):
        if qubit not in cxqubits:
            spectators.append(qubit)
    
    """
    for each qubit in spectators:
        generate a random integer between 0 and 3 (inclusive)
        apply the pauli gate on the spectator 
        apply a barrier
        apply the same pauli gate on the spectator
    """
    

    for qubit in spectators:
        specrand = random.randint(0,3)
        gates[specrand](qubit)
        quantum_circuit.barrier()
        gates[specrand](qubit)

    pass
    
    """
    for each pair in tc_pairs:
        apply each pauli gate for the tc_pair
    """

    for pair in tc_pairs:
        quantum_circuit = twirl_cnot(quantum_circuit,pair[0],pair[1])
    
    pass 
    
    # return the quantum circuit
    
    return quantum_circuit


