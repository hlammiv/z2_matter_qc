# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:10:00 2022
Last edited on Last edited on Thu Mar 24 2022
@authors: Erik Gustafson, Elizabeth Hardt, Norman Hogan,
            Henry Lamm, Ruth Van de Water, and Mike Wagman
Basic Z2 gates.
Used to build circuits in ./Trotterization.py
"""

import qiskit
import numpy as np
import os
import sys
sys.path.append("..")
from scipy.linalg import expm

#import pauli_twirling.circuit_twirling as circuit_twirling

def apply_mass_terms(qc, nsites, mass, epsilon):
    '''
    Apply Rz(mass * epsilon / 2 (-1)^{i}) rotation
    to the even qubits corresponding to the sites.
    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be modified
    nsites : int
        number of sites for the lattice
        (the number of qubits is given by 2 * nsites - 1
         because there are nsites and nsites - 1 gauge links)
    epsilon : float
        the Trotter step size in time (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    Returns
    -------
    qc : qiskit.QuantumCircuit
        the quantum circuit with the sites rotated
    '''

    # Iterate through lattice sites and apply sign rotation.
    for site in range(nsites):
        qc.rz((-1) ** site * mass * epsilon / 2, 2 * site)
    return qc

def apply_gauge_terms(qc, nsites, epsilon):
    '''
    Apply gauge rotations to the odd qubits corresponding to the links.
    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be modified
    nsites : int
        number of sites for the lattice
        (the number of qubits is given by 2 * nsites - 1
         because there are nsites and nsites - 1 gauge links)
    epsilon : float
        the Trotter step size in time (1/a) units
    Returns
    -------
    qc : qiskit.QuantumCircuit
        the quantum circuit with the gauge links rotated
    '''

    # Iterate through nsites-1 gauge links.
    for link in range(nsites - 1):
        qc.rx(epsilon / 2, 2 * link + 1)

    return qc

#####################################################################
############### Erik will add twirling functionality? ###############
#####################################################################
def apply_fermion_hopping(qc, nsites : int, epsilon : float,
                          eta=1.0, twirl=False):
    '''
    Apply fermion hoppping term across an arbitrary connectivity grid and
    number of lattice sites.
    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        the quantum circuit to be modified
    nsites : int
        number of sites for the lattice
        (the number of qubits is given by 2 * nsites - 1
         because there are nsites and nsites - 1 gauge links)
    epsilon : float
        the Trotter step size in time (1/a) units
    eta : float (optional)
        Lattice anisotropy (may be needed to renormalize the speed of light).
        The default is 1.0.
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
        NOT CURRENTLY IMPLEMENTED.
    Returns
    -------
    qc : qiskit.QuantumCircuit
        the quantum circuit with the fermion hopping gates appended
    '''

    # iterate through the even sites
    for site in range(0, nsites - 1, 2):
        q1, q2, q3 = 2 * site, 2 * site + 1, 2 * site + 2
        qc.z(q1)
        qc.sxdg(q2)
        qc.s(q2)
        qc.cx(q1, q2)
        qc.sx(q1)
        qc.cx(q1, q3)
        qc.rx(epsilon * eta / 4, q1)
        qc.ry(epsilon * eta / 4, q3)
        qc.cx(q1, q3)
        qc.sxdg(q1)
        qc.cx(q1, q2)
        qc.sdg(q2)
        qc.sx(q2)
        qc.z(q1)

    # iterate through the odd sites
    for site in range(1, nsites - 1, 2):
        q1, q2, q3 = 2 * site, 2 * site + 1, 2 * site + 2
        qc.z(q1)
        qc.sxdg(q2)
        qc.s(q2)
        qc.cx(q1, q2)
        qc.sx(q1)
        qc.cx(q1, q3)
        qc.rx(-epsilon * eta / 4, q1)
        qc.ry(-epsilon * eta / 4, q3)
        qc.cx(q1, q3)
        qc.sxdg(q1)
        qc.cx(q1, q2)
        qc.sdg(q2)
        qc.sx(q2)
        qc.z(q1)

    return qc

def apply_fermion_hopping_2sites(qc, epsilon, eta=1.0,
                                 twirl=False, richardson_level=1):
    """
    Apply 4-cnot version of the fermion hopping gate shown in the overleaf
    for a 2-site staggered simulation.
    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        quantum circuit for the simulation
    epsilon : float
        the Trotter step size in time (1/a) units
    eta : float
        Lattice anisotropy (may be needed to renormalize the speed of light).
        The default is 1.0.
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    richardson_level : int (optional)
        Determines how many times the CNOT gate is interleaved.
        # N_CNOTS = richardson_level * 2 - 1
    Returns
    -------
    qc : qiskit.QuantumCircuit
        quantum circuit with with fermion hopping gates applied
    """

    # define the number of times to apply a cnot
    ncnots = richardson_level * 2 - 1

    # apply the unitary to semi-diagonalize the fermion hopping gate
    # qc.z(0)
    qc.sxdg(1)
    qc.s(1)

    # if the twirl flag is use applied a cnot twirl
    if twirl:
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 3, [(0, 1)])
    # otherwise just apply a regular cnot
    else:
        for n in range(ncnots):
            qc.barrier(0, 1)
            qc.cx(0, 1)

    # apply the next single qubit gate
    qc.sx(0)

    # if the twirl flag is used apply a cnot twirl
    if twirl:
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 3, [(0, 2)])
    # otherwise just apply a regular cnot
    else:
        for n in range(ncnots):
            qc.barrier(0, 2)
            qc.cx(0, 2)

    # apply the single qubit phase rotations
    qc.rx(-epsilon * eta / 4, 0)
    qc.ry(-epsilon * eta / 4, 2)

    # now we begin the process of inverting the semi-diagonalizing operator
    # twirl the cnot if flag is raised
    if twirl:
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 3, [(0, 2)])
    # otherwise do nothing
    else:
        for n in range(ncnots):
            qc.barrier(0, 2)
            qc.cx(0, 2)

    # single qubit gate
    qc.sxdg(0)

    # twirl the CNOT from p to gamma if flag raised
    if twirl:
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 3, [(0, 1)])
    # otherwise do nothing
    else:
        for n in range(ncnots):
            qc.barrier(0, 1)
            qc.cx(0, 1)

    # final single qubit gates
    qc.sdg(1)
    qc.sx(1)
    # qc.z(0)

    return qc

#####################################################################
############### Inline comments in code hard to follow ##############
#####################################################################
def apply_fermion_hopping_4sites(qc, epsilon, eta=1.0,
                                 twirl=False, richardson_level=1):
    """
    Apply 4-cnot version of the fermion hopping gate shown in the overleaf
    for a 4-site staggered simulation.
    Parameters
    ----------
    qc : qiskit.QuantumCircuit
        quantum circuit for the simulation
    epsilon : float
        the Trotter step size in time (1/a) units
    eta : float
        Lattice anisotropy (may be needed to renormalize the speed of light).
        The default is 1.0.
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
    richardson_level : int (optional)
        Determines how many times the CNOT gate is interleaved.
        # N_CNOTS = richardson_level * 2 - 1
    Returns
    -------
    qc : qiskit.QuantumCircuit
       quantum circuit with with fermion hopping gates applied
    """

    # define the number of times to apply a cnot
    ncnots = richardson_level * 2 - 1

    # single qubit rotations for the hopping from sites 0-1 and 2-3
    #================#
    # fermion term 1 #
    #================#
    qc.z([0, 4])
    qc.sxdg([1, 5])
    qc.s([1, 5])
    # if twirl flag is raised apply a twirled CNOT HARD CYCLE
    # (simultaneously applied )
    #=========================
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(0, 1), (4, 5)])
    # otherwise just a regular cnot across (2, 1) and (4, 5)
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            qc.barrier(0, 1)
            qc.barrier(4, 5)
            qc.cx(0, 1)
            qc.cx(4, 5)
    # single qubit gate rotation
    #================================
    qc.sx([0, 4])
    # if twirling flag is raised apply the second CNOT gate
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(0, 2), (4, 6)])
    # otherwise just apply a plain old cnot across 2-0 and 4-6
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            qc.barrier(0, 2)
            qc.barrier(4, 6)
            qc.cx(0, 2)
            qc.cx(4, 6)
    #single qubit rotions for the evolution
    #===============================
    qc.rx(epsilon / 2, [0, 4])
    qc.ry(epsilon / 2, [2, 6])
    # let's un do the semi-diagonalization operations across the even pairs

    # if twirling flag is raised apply the second CNOT gate
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(0, 2), (4, 6)])
    # otherwise just apply a plain old cnot across 2-0 and 4-6
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            qc.barrier(0, 2)
            qc.barrier(4, 6)
            qc.cx(0, 2)
            qc.cx(4, 6)
    #===============================
    qc.sxdg([0, 4])
    #=========================
    #=========================
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(0, 1), (4, 5)])
    # otherwise just a regular cnot across (2, 1) and (4, 5)
    else:
        # iterate the cnots
        for n in range(ncnots):
            qc.barrier(0, 1)
            qc.barrier(4, 5)
            qc.cx(0, 1)
            qc.cx(4, 5)
    #================================
    qc.sdg([1, 5])
    qc.sx([1, 5])
    qc.z([0, 4])
    # now we are going to apply the fermion hopping term to the sites 1 - 2
    #===========#
    # Fhop pt2  #
    #===========#
    qc.h(3)
    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:

        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 2)])

        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 4)])
    else:
        for n in range(ncnots):
            qc.barrier(3, 2)
            qc.cx(3, 2)
        for n in range(ncnots):
            qc.barrier(3, 4)
            qc.cx(3, 4)
    # some more single qubit rotations
    #================================
    qc.rx(-epsilon / 2, 3)
    qc.h(2)
    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:

        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 2)])
    else:
        for n in range(ncnots):
            qc.barrier(3, 2)
            qc.cx(3, 2)
    #================================
    qc.h(2)
    qc.s(2)
    qc.z(3)
    qc.h(4)

    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:
        for n in range(ncnots):
            qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 4)])
    else:
        for n in range(ncnots):
            qc.barrier(3, 4)
            qc.cx(3, 4)
    qc.h(4)
    qc.s(4)
    qc.rx(-epsilon / 2, 3)

    # if the twirling flag is raised apply a cnot hard cycle twirl here
    if twirl:
        qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 4)])
        qc = circuit_twirling.twirl_hard_cycle(qc, 7, [(3, 2)])
    # otherwise just leave it as is
    else:
        qc.cx(3, 4)
        qc.cx(3, 2)

    qc.h(3)
    qc.sdg([2, 4])

    return qc
