# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:10:00 2022

@author: Erik Gustafson, Henry Lamm, Ruth Van der Water
Mike Wagman

This script file contains the various gates that are used in the
simulations
"""

import qiskit
import numpy as np
import os
import sys
from scipy.linalg import expm
sys.path.append("..")
import pauli_twirling.circuit_twirling as circuit_twirling
#from ..pauli_twirling.circuit_twirling import paulitwirling

def apply_mass_terms(quantum_circuit, number_of_sites,
                      mass, epsilon):
    '''

    this function applies an Rz(mass * epsilon / 2 (-1)^{i})
    rotation to the qubits corresponding to the sites

    Parameters
    ----------
    quantum_circuit : Qiskit QuantumCircuit Object
        quantum circuit for the simulation
    number_of_sites : Int
        number of sites in the simulation
    mass : float
        The mass of fermion
    epsilon : float
        The suzuki Trotter step size

    Returns
    -------
    quantum_circuit : qiskit.QuantumCircuit Object
    '''
    # Here we iterate through the sites, 1 to N,
    # and then apply the appropriate sign rotation
    for site in range(number_of_sites):
        # the sites correspond to every other qubit
        quantum_circuit.rz((-1) ** site * mass * epsilon / 2,
                            2 * site)
    # Let's return the quantum circuit
    # this is not strictly necessary but is safe
    return quantum_circuit


def apply_gauge_terms(quantum_circuit, number_of_sites, epsilon):
    '''
    Apply the gauge rotations to the qubits corresponding to the links
    on the lattice

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit object
        quantum circuit for the simulation
    number_of_sites : int
        number of sites for the quantum circuit
    epsilon : float
        the time evolution operator

    Returns
    -------
    quantum_circuit

    '''
    # iterate through the odd qubits which correspond to the gauge links
    for link in range(number_of_sites - 1):
        quantum_circuit.rx(epsilon / 2, 2 * link + 1)
    # return the quantum circuit
    return quantum_circuit


def apply_fermion_hopping(quantum_circuit, number_of_sites : int,
                          epsilon : float,
                          eta=1.0, twirl=False):
    '''
    apply a fermion hoppping term across an arbitrary connectivity grid and
    number of lattice sites.

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit
        quantum circuit to add the fermion hopping term to
    number_of_sites : int
        number of physical lattice sites
    epsilon : float
        the Suzuki-Trotter step size
    eta : float, optional
        anisotropy (maybe needed to renormalize the speed of light).
        The default is 1.0.
    twirl : boolean, optional
        whether to twirl or not. NOT CURRENTLY IMPLEMENTED.
        The default is False.

    Returns
    -------
    quantum_circuit : qiskit.QuantumCircuit
        the quantum circuit with the fermion hopping gates appended

    '''
    # iterate through the even sites
    for site in range(0, number_of_sites - 1, 2):
        q1, q2, q3 = 2 * site, 2 * site + 1, 2 * site + 2
        quantum_circuit.z(q1)
        quantum_circuit.sxdg(q2)
        quantum_circuit.s(q2)
        quantum_circuit.cx(q1, q2)
        quantum_circuit.sx(q1)
        quantum_circuit.cx(q1, q3)
        quantum_circuit.rx(epsilon * eta / 4, q1)
        quantum_circuit.ry(epsilon * eta / 4, q3)
        quantum_circuit.cx(q1, q3)
        quantum_circuit.sxdg(q1)
        quantum_circuit.cx(q1, q2)
        quantum_circuit.sdg(q2)
        quantum_circuit.sx(q2)
        quantum_circuit.z(q1)
    # iterate through the odd sites
    for site in range(1, number_of_sites - 1, 2):
        q1, q2, q3 = 2 * site, 2 * site + 1, 2 * site + 2
        quantum_circuit.z(q1)
        quantum_circuit.sxdg(q2)
        quantum_circuit.s(q2)
        quantum_circuit.cx(q1, q2)
        quantum_circuit.sx(q1)
        quantum_circuit.cx(q1, q3)
        quantum_circuit.rx(-epsilon * eta / 4, q1)
        quantum_circuit.ry(-epsilon * eta / 4, q3)
        quantum_circuit.cx(q1, q3)
        quantum_circuit.sxdg(q1)
        quantum_circuit.cx(q1, q2)
        quantum_circuit.sdg(q2)
        quantum_circuit.sx(q2)
        quantum_circuit.z(q1)
    return quantum_circuit


def apply_fermion_hopping_2sites(quantum_circuit, epsilon,
                                  eta=1.0, twirl=False,
                                  richardson_level=1):
    """

    apply the 4 cnot version of the fermion hopping gate shown in the overleaf
    for a two site simulation

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit
        quantum circuit for the simulation
    epsilon : float
        suzuki Trotter step size
    eta : float
        an scale for renormalization of speed of light (?)
    twirl : boolean (optional)
        whether to use twirled or not twirled gates
    richardson_level: int (optional)
        the richardson scaling level

    Returns
    -------
    quantum_circuit : qiskit.QuantumCircuit
        quantum circuit with with fermion hopping gates applied

    """

    # define the number of times to apply a cnot
    ncnots = richardson_level * 2 - 1
    # apply the unitary to semi-diagonalize the fermion hopping gate
    quantum_circuit.z(0)
    quantum_circuit.sxdg(1)
    quantum_circuit.s(1)
    # if the twirl flag is use applied a cnot twirl
    if twirl:
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                3, [(0, 1)])
    # otherwise just apply a regular cnot
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(0, 1)
            quantum_circuit.cx(0, 1)
    # apply the next single qubit gate
    quantum_circuit.sx(0)
    # if the twirl flag is used apply a cnot twirl
    if twirl:
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                            3, [(0, 2)])
    # otherwise just apply a regular cnot
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(0, 2)
            quantum_circuit.cx(0, 2)
    # apply the single qubit phase rotations
    quantum_circuit.rx(epsilon * eta / 4, 0)
    quantum_circuit.ry(epsilon * eta / 4, 2)
    # now we begin the process of inverting the semi-diagonalizing operator
    # twirl the cnot if flag is raised
    if twirl:

        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                3, [(0, 2)])
    # otherwise do nothing
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(0, 2)
            quantum_circuit.cx(0, 2)
    # single qubit gate
    quantum_circuit.sxdg(0)
    # twirl the CNOT from p to gamma if flag raised
    if twirl:
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                            3, [(0, 1)])
    # otherwise do nothing
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(0, 1)
            quantum_circuit.cx(0, 1)
    # final single qubit gates
    quantum_circuit.sdg(1)
    quantum_circuit.sx(1)
    quantum_circuit.z(0)

    return quantum_circuit


def apply_fermion_hopping_4sites(quantum_circuit, epsilon, eta=1.0,
                                 twirl=False, richardson_level=1):
    """
      perform the fermion hopping operation across a 4 site staggered
      fermion lattice

    Parameters
    ----------
    quantum_circuit : qiskit.QuantumCircuit
        the quantum circuit we are going to be applying gates to
    epsilon : float
        the suzuki trotter step size
    eta : int, optional
        scale setting for the renomalization of the speed of light.
        . The default is 1.0.
    twirl : boolean, optional
        whether to pauli twirl (randomized compiling)
        or not. The default is False.
    richardson_level : int, optional
        how much to distort the cnot gates
    Returns
    -------
    quantum_circuit : qiskit.QuantumCircuit
        the quantum_circuit object with the gates applied.

    """
    # this the relation to distort the cnot gates by for a given richardson
    # level
    ncnots = richardson_level * 2 - 1
    # single qubit rotations for the hopping from sites 0-1 and 2-3
    #================#
    # fermion term 1 #
    #================#
    quantum_circuit.z([0, 4])
    quantum_circuit.sxdg([1, 5])
    quantum_circuit.s([1, 5])
    # if twirl flag is raised apply a twirled CNOT HARD CYCLE
    # (simultaneously applied )
    #=========================
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(0, 1),
                                                                    (4, 5)])
    # otherwise just a regular cnot across (2, 1) and (4, 5)
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            quantum_circuit.barrier(0, 1)
            quantum_circuit.barrier(4, 5)
            quantum_circuit.cx(0, 1)
            quantum_circuit.cx(4, 5)
    # single qubit gate rotation
    #================================
    quantum_circuit.sx([0, 4])
    # if twirling flag is raised apply the second CNOT gate
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(0, 2),
                                                                    (4, 6)])
    # otherwise just apply a plain old cnot across 2-0 and 4-6
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            quantum_circuit.barrier(0, 2)
            quantum_circuit.barrier(4, 6)
            quantum_circuit.cx(0, 2)
            quantum_circuit.cx(4, 6)
    #single qubit rotions for the evolution
    #===============================
    quantum_circuit.rx(epsilon / 2, [0, 4])
    quantum_circuit.ry(epsilon / 2, [2, 6])
    # let's un do the semi-diagonalization operations across the even pairs

    # if twirling flag is raised apply the second CNOT gate
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(0, 2),
                                                                    (4, 6)])
    # otherwise just apply a plain old cnot across 2-0 and 4-6
    else:
        # iterate the untwirled cnots
        for n in range(ncnots):
            quantum_circuit.barrier(0, 2)
            quantum_circuit.barrier(4, 6)
            quantum_circuit.cx(0, 2)
            quantum_circuit.cx(4, 6)
    #===============================
    quantum_circuit.sxdg([0, 4])
    #=========================
    #=========================
    if twirl:
        # iterate the twirled cnots
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(0, 1),
                                                                    (4, 5)])
    # otherwise just a regular cnot across (2, 1) and (4, 5)
    else:
        # iterate the cnots
        for n in range(ncnots):
            quantum_circuit.barrier(0, 1)
            quantum_circuit.barrier(4, 5)
            quantum_circuit.cx(0, 1)
            quantum_circuit.cx(4, 5)
    #================================
    quantum_circuit.sdg([1, 5])
    quantum_circuit.sx([1, 5])
    quantum_circuit.z([0, 4])
    # now we are going to apply the fermion hopping term to the sites 1 - 2
    #===========#
    # Fhop pt2  #
    #===========#
    quantum_circuit.h(3)
    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:

        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(3, 2)])

        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(3, 4)])
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(3, 2)
            quantum_circuit.cx(3, 2)
        for n in range(ncnots):
            quantum_circuit.barrier(3, 4)
            quantum_circuit.cx(3, 4)
    # some more single qubit rotations
    #================================
    quantum_circuit.rx(-epsilon / 2, 3)
    quantum_circuit.h(2)
    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:

        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(3, 2)])
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(3, 2)
            quantum_circuit.cx(3, 2)
    #================================
    quantum_circuit.h(2)
    quantum_circuit.s(2)
    quantum_circuit.z(3)
    quantum_circuit.h(4)

    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:
        for n in range(ncnots):
            quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                                7, [(3, 4)])
    else:
        for n in range(ncnots):
            quantum_circuit.barrier(3, 4)
            quantum_circuit.cx(3, 4)
    quantum_circuit.h(4)
    quantum_circuit.s(4)
    quantum_circuit.rx(-epsilon / 2, 3)
    # if the twirling flag is raised apply a cnot hard cycle twirl here
    # otherwise just leave it as is
    if twirl:
        quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                            7, [(3, 4)])
        quantum_circuit = circuit_twirling.twirl_hard_cycle(quantum_circuit,
                                                            7, [(3, 2)])
    else:
        quantum_circuit.cx(3, 4)
        quantum_circuit.cx(3, 2)
    quantum_circuit.h(3)
    quantum_circuit.sdg([2, 4])
    # return the quantum circuit object
    return quantum_circuit


