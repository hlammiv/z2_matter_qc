# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:06:12 2022

@author: Erik Gustafson, Henry Lamm, Ruth Van der Water
Mike Wagman, Norman, Clement, Sarah, Elizabeth
"""

from . import Z2Gates


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
    quantum_circuit = Z2Gates.apply_mass_terms(quantum_circuit, 2, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2Gates.apply_gauge_terms(quantum_circuit, 2, epsilon)
    # apply the fermion hopping term across the qubits
    quantum_circuit = Z2Gates.apply_fermion_hopping_2sites(quantum_circuit, epsilon, eta=1.0,
                                         twirl=twirl)
    # apply the rotations corresponding to the mass operators
    quantum_circuit = Z2Gates.apply_mass_terms(quantum_circuit, 2, mass, epsilon)
    # apply the rotations corresponding to the gauge dynamics operators
    quantum_circuit = Z2Gates.apply_gauge_terms(quantum_circuit, 2, epsilon)

    return quantum_circuit

