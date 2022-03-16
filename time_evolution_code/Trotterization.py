# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:06:12 2022

@author: Erik Gustafson, Henry Lamm, Ruth Van der Water
Mike Wagman, Norman, Clement, Sarah, Elizabeth
"""

import Z2gates
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator as Operator


def trotter_evolution_2site(quantum_circuit, epsilon: float, mass: float,
                            ntrotter: int, twirl=False):
    """
    apply the Suzuki-Trotter evolution operator ntrotter times on the quantum
    circuit:

        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2) U_{matter-gauge}(dt)
        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2)

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit
        quantum circuit for time evolution
    epsilon : float
        the Trotter step size
    mass : float
        the fermion mass
    ntrotter : int
        the number of times the trotter operator is applied
    twirl : boolean (optional)
        whether to apply pauli twirling to the cnot gates

    Returns
    -------
    quantum_circuit : TYPE
        quantum circuit object after the gates have been applied

    """
    # apply the rotations corresponding to the mass operators
    quantum_circuit = Z2gates.apply_mass_terms(quantum_circuit, 2, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2gates.apply_gauge_terms(quantum_circuit, 2, epsilon)
    # apply the fermion hopping term across the qubits
    quantum_circuit = Z2gates.apply_fermion_hopping_2sites(quantum_circuit,
                                                           epsilon, eta=1.0,
                                                           twirl=twirl)
    # apply the rotations corresponding to the mass operators
    quantum_circuit = Z2gates.apply_mass_terms(quantum_circuit, 2, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2gates.apply_gauge_terms(quantum_circuit, 2, epsilon)

    return quantum_circuit

def trotter_evolution_4sites(quantum_circuit, epsilon: float, mass: float,
                             ntrotter: int, twirl=False):
    """
    apply the Suzuki-Trotter evolution operator ntrotter times on the quantum
    circuit:

        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2) U_{matter-gauge}(dt)
        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2)

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit
        quantum circuit for time evolution
    epsilon : float
        the Trotter step size
    mass : float
        the fermion mass
    ntrotter : int
        the number of times the trotter operator is applied
    twirl : boolean (optional)
        whether to apply pauli twirling to the cnot gates

    Returns
    -------
    quantum_circuit : TYPE
        quantum circuit object after the gates have been applied

    """
    # apply the rotations corresponding to the mass operators
    quantum_circuit = Z2gates.apply_mass_terms(quantum_circuit, 4, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2gates.apply_gauge_terms(quantum_circuit, 4, epsilon)
    # apply the fermion hopping term across the qubits
    quantum_circuit = Z2gates.apply_fermion_hopping_4sites(quantum_circuit,
                                                            epsilon, eta=1.0,
                                                            twirl=twirl)
    # apply the rotations corresponding to the mass operators
    quantum_circuit = Z2gates.apply_mass_terms(quantum_circuit, 4, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2gates.apply_gauge_terms(quantum_circuit, 4, epsilon)

    return quantum_circuit



if __name__ == "__main__":
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    fhop = (np.kron(x, np.kron(z, x)) + np.kron(y, np.kron(z, y))) / 4
    dt = 0.1
    mass = 1.0
    fhop = expm(-1.0j * dt * fhop)
    rmass = expm(0.25j * dt * mass * z)
    rgauge = expm(-0.25j * dt * x)
    rgauges = [np.kron(np.identity(2 ** i),
                       np.kron(rgauge, np.identity(2 ** (6 - i))))
               for i in range(1, 7, 2)]
    rmasses = [np.kron(np.identity(2 ** i),
                       np.kron(rmass, np.identity(2 ** (6 - i))))
               for i in range(0, 7, 2)]

    r1q = np.identity(2 ** 7)
    for op in rgauges:
        r1q = op @ r1q
    for j in range(4):
        if j % 2 == 1:
            r1q = rmasses[j].conjugate().transpose() @ r1q
        else:
            r1q = rmasses[j] @ r1q
    fhops = [np.kron(fhop, np.identity(2 ** 4)),
             np.kron(np.identity(16), fhop),
             np.kron(np.identity(4), np.kron(fhop, np.identity(4)))]

    oper = np.identity(2 ** 7).dot(r1q)
    for op in fhops:
        oper = op @ oper
    oper = r1q @ oper
    oper = Operator(oper)
    opert = Operator(trotter_evolution_4sites(QuantumCircuit(7), dt, mass, 1))
    # print(np.array(oper)[:10, :10], np.array(opert)[:10, :10])
    print(oper.equiv(opert))