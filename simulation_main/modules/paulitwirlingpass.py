# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:57:40 2022

@author: Erik Gustafson
"""

from qiskit import QuantumRegister, transpile
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.dagcircuit import DAGCircuit
from qiskit.circuit.library import CXGate, XGate, YGate, ZGate, IGate
import numpy as np



class RandomizedCompile(TransformationPass):
    """Transforms all the CNOTs in circuit into different representations
    of CNOTs by dressing the circuit with randomly selected Pauli gates"""

    def __init__(self):

        super().__init__()
        self.twirldict = {(0, 0): (0, 0), (0, 1): (0, 1), (0, 2): (3, 2),
                          (0, 3): (3, 3), (1, 0): (1, 1), (1, 1): (1, 0),
                          (1, 2): (2, 3), (1, 3): (2, 2), (2, 0): (2, 1),
                          (2, 1): (2, 0), (2, 2): (1, 3), (2, 3): (1, 2),
                          (3, 0): (3, 0), (3, 1): (3, 1), (3, 2): (0, 2),
                          (3, 3): (0, 3)}
        # pauli gates
        self.gates = [IGate(), XGate(), YGate(), ZGate()]
        # list of equivalent DAG circuits
        self.twirls = []
        for i in range(16):
            a = i // 4
            b = i % 4
            ap, bp = self.twirldict[(a, b)]
            q = QuantumRegister(2)
            ap, bp = self.twirldict[(a, b)]
            mini_dag = DAGCircuit()
            mini_dag.add_qreg(q)
            if a > 0:
                mini_dag.apply_operation_back(self.gates[a], qargs=[q[0]])
            if b > 0:
                mini_dag.apply_operation_back(self.gates[b], qargs=[q[1]])
            mini_dag.apply_operation_back(CXGate(), qargs=[q[0], q[1]])
            if ap > 0:
                mini_dag.apply_operation_back(self.gates[ap], qargs=[q[0]])
            if bp > 0:
                mini_dag.apply_operation_back(self.gates[bp], qargs=[q[1]])
            self.twirls.append((mini_dag, q))

    def run(self, new_dag):
        """Runs the RandomizedCompile pass on 'dag'.

        args:
            dag (DAGCircuit): DAG to randomly compile

        Returns:
            DAGCircuit: a randomly compiled DAG
        """
        cx_list = new_dag.op_nodes(op=CXGate)
        for j in range(len(cx_list)):
            node = cx_list[j]
            number = np.random.randint(0, 16)
            mini_dag, reg = self.twirls[number]
            new_dag.substitute_node_with_dag(node=node, input_dag=mini_dag,
                                             wires=[reg[0], reg[1]])
        return new_dag

def randomly_compile(circuit, ncopy=20, backend=None, initial_layout=None):
    """ generate ncopy random compiled circuits

    """
    pm = PassManager()
    pm.append(RandomizedCompile())
    output_circuits = []
    for i in range(ncopy):
        out_circ = pm.run(circuit.copy())
        output_circuits.append(out_circ)
    if backend is not None:
        if initial_layout is None:
            output_circuits = transpile(output_circuits, backend=backend, basis_gates=['cx', 'sx', 'rz'],
                                                optimization_level=1)
        else:
            output_circuits = transpile(output_circuits, backend=backend, basis_gates=['cx', 'sx', 'rz'],
                                                optimization_level=1, initial_layout=initial_layout)

    return output_circuits


if __name__ == "__main__":
    qc = QuantumCircuit(4)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.cx(2, 3)
    transpiled = transpile(qc, basis_gates=['sx', 'rz', 'cx'])
    print(transpiled)
    twirled = randomly_compile(transpiled.copy())
    print(twirled[0])


