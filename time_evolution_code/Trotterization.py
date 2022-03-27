# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:06:12 2022
Last edited on Thu Mar 26 2022

@authors: Erik Gustafson, Elizabeth Hardt, Norman Hogan,
          Henry Lamm, Ruth Van de Water, Mike Wagman
          
Build circuits for Z2 gauge theory simulation
using gates defined in ./Z2gates.py
"""

import Z2gates
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator as Operator


def trotter_evolution(num_sites: int, epsilon: float, mass: float,
                      ntrotter: int, twirl=False, qsim=True,
                      richardson_level: int=1, state_prep=True):
    """
    Wrapper function for the various types of noiseless, noisy, and harware
    simulations that we want to do for the Z2 gauge with staggered matter
    in (1+1)D. These simulations allow for including and excluding randomized
    compiling (pauli twirling) and richardson extrapolations.
    
    Parameters
    ----------
    nsites : int
        number of sites for the lattice
        (the number of qubits is given by 2 * nsites - 1
         because there are nsites and nsites - 1 gauge links)
    epsilon : float
        the Trotter step size in time (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    qsim : boolean, optional
        Whether to set this simulation up for running on a quantum computer
        or on a qasm/statevector simulator. The default is True
    richardson_level : int (optional)
        Determines how many times the CNOT gate is interleaved.
        # N_CNOTS = richardson_level * 2 - 1
    state_prep : boolean, (optional)
        whether to use a state prep circuit or not. Default is True
        
    Returns
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit that will be simulated
    """

    # get number of qubits from number of sites
    nqubits = 2 * nsites - 1
    
    # create a quantum circuit
    qc = QuantumCircuit(nqubits, nqubits)

    ###############################################################
    # prepare initial state -- eventually move to separate function
    ###############################################################
    # Hademard operator on gauge links puts photon qubits in + state of x basis 
    # (recall that H|0> = |+>)
    # needed for gauge invariance? Vacuum quantum numbers? Something else???
    if state_prep:
        qc.h([2 * i + 1 for i in range(num_sites - 1)])

    # build quantum circuit for running on a quantum computer
    # using different layouts for 2 sites (3 qubits) and 4 sites (7 qubits)
    if (qsim):
        if nsites == 2:
            if state_prep:
                qc.x([0])
                qc.z([1])
            qc = trotter_evolution_2sites(qc, epsilon, mass, ntrotter,
                                          twirl=twirl, richardson_level=richardson_level)
        elif nsites == 4:
            if state_prep:
                qc.x([2, 4])
                qc.z([3])
            qc = trotter_evolution_4sites(qc, epsilon, mass, ntrotter,
                                          twirl=twirl, richardson_level=richardson_level)
        else:
            exception_str = "For quantum hardware,"
            exception_str += "only 2- and 4-site simulation code is available."
            raise Exception(exception_str)
            
    # build quantum circuit for running on a simulator backend
    else:
        if state_prep:
            mid = 2 * (nsites // 2)
            qc.x([mid, mid + 2])
            qc.z(mid + 1)
        qc = trotter_evolution_generic(qc, nsites, epsilon, mass, ntrotter,
                                       twirl=False, save_state_vector=True)

    return qc


def trotter_evolution_generic(qc, nsites: int, epsilon: float,
                              mass: float, ntrotter: int, twirl=False,
                              save_state_vector=True):
    """
    Quantum circuit for a second-order Trotter using the
    most efficient encoding of the Trotterization.
    Use this when running on a simulator backend.
    
    Apply the Suzuki-Trotter evolution operator ntrotter times on the quantum
    circuit:

        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2) U_{matter-gauge}(dt)
        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2)

    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be simulated
    nsites : int
        number of sites for the lattice
        (the number of qubits is given by 2 * nsites - 1
         because there are nsites and nsites - 1 gauge links)
    epsilon : float
        the Trotter step size in time (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    save_state_vector : boolean (optional)
        save the statevectors of the simulation. This is good if we want to
        minimize the computing time. The default is True.
        
    Returns
    ----------
    qc : qiskit.QuantumCircuit
        quantum circuit object after the gates have been applied
    """

    # apply second order Suzuki-Trotter operator ntrotter times
    for step in range(ntrotter):
        
        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, nsites, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, nsites, epsilon)

        # apply the fermion hopping term across the qubits
        qc = Z2gates.apply_fermion_hopping(qc, nsites, epsilon, eta=1.0, twirl=False)

        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, nsites, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, nsites, epsilon)
        
        # save the state vector
        if (save_state_vector):
            qc.save_statevector(label=str(step))
        
    return qc


def trotter_evolution_2sites(qc, epsilon: float, mass: float, ntrotter: int,
                             twirl=False, richardson_level=1):
    """
    Quantum circuit for a second-order Trotter using an encoding of the
    Trotterization suitable for running a 2-site simulation on a quantum computer.
    
    Apply the Suzuki-Trotter evolution operator ntrotter times on the quantum
    circuit:

        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2) U_{matter-gauge}(dt)
        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2)

    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be simulated
    epsilon : float
        the Trotter step size in time (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    richardson_level : int (optional)
        Determines how many times the CNOT gate is interleaved.
        # N_CNOTS = richardson_level * 2 - 1
        
    Returns
    ----------
    qc : qiskit.QuantumCircuit
        quantum circuit object after the gates have been applied
    """

    # apply the Suzuki-Trotter operator ntrotter times
    for step in range(ntrotter):
        
        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, 2, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, 2, epsilon)

        # apply the fermion hopping term across the qubits
        qc = Z2gates.apply_fermion_hopping_2sites(qc, epsilon, eta=1.0,
                                                  twirl=twirl, richardson_level=richardson_level)

        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, 2, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, 2, epsilon)

    return qc


def trotter_evolution_4sites(qc, epsilon: float, mass: float, ntrotter: int,
                             twirl=False, richardson_level=1):
    """
    Quantum circuit for a second-order Trotter using an encoding of the
    Trotterization suitable for running a 4-site simulation on a quantum computer.
    
    Apply the Suzuki-Trotter evolution operator ntrotter times on the quantum
    circuit:

        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2) U_{matter-gauge}(dt)
        U_{pure_gauge}(dt / 2) U_{matter}(dt / 2)

    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be simulated
    epsilon : float
        the Trotter step size in time (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    richardson_level : int (optional)
        Determines how many times the CNOT gate is interleaved.
        # N_CNOTS = richardson_level * 2 - 1
        
    Returns
    ----------
    qc : qiskit.QuantumCircuit
        quantum circuit object after the gates have been applied
    """
    
    # apply the Suzuki-Trotter Operator ntrotter times
    for step in range(ntrotter):
        
        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, 4, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, 4, epsilon)

        # apply the fermion hopping term across the qubits
        qc = Z2gates.apply_fermion_hopping_4sites(qc, epsilon, eta=1.0, twirl=twirl,
                                                  richardson_level=richardson_level)

        # apply the rotations corresponding to the mass operators
        qc = Z2gates.apply_mass_terms(qc, 4, mass, epsilon)

        # apply the rotations corresponding to the gauge dynamics operators
        qc = Z2gates.apply_gauge_terms(qc, 4, epsilon)

    return qc


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
             np.kron(np.identity(4),
                     np.kron(fhop.conjugate().transpose(), np.identity(4)))]

    oper = np.identity(2 ** 7).dot(r1q)
    for op in fhops:
        oper = op @ oper
    oper = r1q @ oper
    oper = Operator(oper)
    opert = Operator(trotter_evolution_4sites(QuantumCircuit(7), dt, mass, 1))
    # print(np.array(oper)[:10, :10], np.array(opert)[:10, :10])
    print(oper.equiv(opert))
    print('=' * 20)
    try:
        trotter_evolution(8, 0.1, 1.0, 1, twirl=False, qsim=True,
                          state_prep=False)
    except (Exception) as e:
        print('====' * 5)
        print(e)
        print('---')
        print('exception working')
    print('=' * 100)
    trotter_evolution(8, 0.1, 1.0, 1, twirl=False, qsim=False,
                      state_prep=False)

