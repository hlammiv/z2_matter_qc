# functions for creating a circuit for J1 current insertion. 
# author: Mariia Kharchenko harchenkomv@icloud.com, mariiak@fnal.gov 
# note: this code is work in progress 

import qiskit as qis
import numpy as np
import matplotlib.pyplot as plt

def create_init_circuit(nqubits,vector):
    qc = qis.QuantumCircuit(nqubits)
    qc.h(0) #bring ancilla into |+> state
    qc.initialize(vector,[i for i in range(1,nqubits)]) #place-holder until state prep group gives us init code
    return qc


def XzY_vm1_v(nqubits,v):
    # ancilla is the 0th qubit
    qc = qis.QuantumCircuit(nqubits)
    qc.cx(0,v-2) # v-1 space on lattice
    qc.cz(0,v-1) # v-1,v gauge
    qc.cy(0,v)   # v

    return qc

def XzY_vp1_v(nqubits,v):
    # ancilla is the 0th qubit
    qc = qis.QuantumCircuit(nqubits)
    qc.cx(0,v)   # v space on lattice
    qc.cz(0,v+1) # v,v+1 gauge space on lattice
    qc.cy(0,v+2) # v+1 space on lattice

    return qc

def YzX_vm1_v(nqubits,v):
    # ancilla is the 0th qubit
    qc = qis.QuantumCircuit(nqubits)
    qc.cy(0,v-2) # v-1 space on lattice
    qc.cz(0,v-1) # v-1,v (gauge) space on lattice
    qc.cx(0,v)   # v space on lattice

    return qc

def YzX_vp1_v(nqubits,v):
    # ancilla is the 0th qubit
    qc = qis.QuantumCircuit(nqubits)
    qc.cy(0,v)   # v space on lattice
    qc.cz(0,v+1) # v,v+1 gauge space on lattice
    qc.cx(0,v+2) # v+1 space on lattice

    return qc




def to_XZY_vm1_v(nqubits,v_site):
    qc = qis.QuantumCircuit(nqubits)
    qc.h(0) # bring ancilla back to the initial state

    qc.h(v_site-2) # apply hadomard gate to v-1 lattice site to bring it to X basis

    # do nothing with v-1,v site because it is already in computational basis

    qc.sdg(v_site) # step 1 of brining vth site on the lattice to Y-basis
    qc.h(v_site) # step 2 of brining vth site on the lattice to Y-basis

    return qc

def to_XZY_vp1_v(nqubits,v_site):
    qc = qis.QuantumCircuit(nqubits)
    qc.h(0) # bring ancilla back to the initial state

    qc.h(v_site) # apply hadomard gate to vth lattice site to bring it to X basis

    # do nothing with v,v+1 site (link) because it is already in computational basis

    qc.sdg(v_site+2) # step 1 of brining v+1th site on the lattice to Y-basis
    qc.h(v_site+2) # step 2 of brining v+1th site on the lattice to Y-basis

    return qc


def to_YZX_vm1_v(nqubits,v_site):
    qc = qis.QuantumCircuit(nqubits)

    qc.h(0) # bring ancilla back to the initial state

    qc.sdg(v_site-2) # step 1 of brining v-1th site on the lattice to Y-basis
    qc.h(v_site-2) # step 2 of brining v-1th site on the lattice to Y-basis

    # do nothing with v-1,v site because it is already in computational basis

    qc.h(v_site) # apply hadomard gate to Vth lattice site to bring it to X basis

    return qc

def to_YZX_vp1_v(nqubits,v_site):
    qc = qis.QuantumCircuit(nqubits)

    qc.h(0) # bring ancilla back to the initial state

    qc.sdg(v_site) # step 1 of brining vth site on the lattice to Y-basis
    qc.h(v_site) # step 2 of brining vth site on the lattice to Y-basis

    # do nothing with v,v+1 site (link) because it is already in computational basis

    qc.h(v_site+2) # apply hadomard gate to v+1th lattice site to bring it to X basis


    return qc


def apply_J1(nqubits,qc,v,which_order): # v is the site on lattice
    """
    this function applies J1 to a given quantum circuit qc
    which_order:
    variable is used to indicate which part of J1 you want to apply:
      > which_order == 0: XzY for v-1,v
      > which_order == 1: XzY for v,v+1
      > which_order == 2: YzX for v-1,v
      > which_order == 3: YzX for v,v+1
    """



    if which_order == 0:

        qc_XzY_vm1_v = XzY_vm1_v(nqubits,v)
        qc.append(qc_XzY_vm1_v,[i for i in range(nqubits)])
    if which_order == 1:

        qc_XzY_vp1_v = XzY_vp1_v(nqubits,v)
        qc.append(qc_XzY_vp1_v,[i for i in range(nqubits)])

    if which_order == 2:

        qc_YzX_vm1_v = YzX_vm1_v(nqubits,v)
        qc.append(qc_YzX_vm1_v,[i for i in range(nqubits)])

    if which_order == 3:

        qc_YzX_vp1_v = YzX_vp1_v(nqubits,v)
        qc.append(qc_YzX_vp1_v,[i for i in range(nqubits)])

    return qc

def transform_basis(nqubits, qc,v_site,which_basis):
    """
    this function transforms the given quantum circuit qc to an appropriate basis for measuremnet
    which_basis:
    variable is used to indicate which part of J1 you want to apply:
      > which_basis == 0: XzY for v-1,v
      > which_basis == 1: XzY for v,v+1
      > which_basis == 2: YzX for v-1,v
      > which_basis == 3: YzX for v,v+1
    """

    if which_basis == 0:

        basis_qc_XzY_vm1_v = to_XZY_vm1_v(nqubits,v_site)
        qc.append(basis_qc_XzY_vm1_v,[i for i in range(nqubits)])

    if which_basis == 1:

        basis_qc_XzY_vp1_v = to_XZY_vp1_v(nqubits,v_site)
        qc.append(basis_qc_XzY_vp1_v,[i for i in range(nqubits)])

    if which_basis == 2:

        basis_qc_YzX_vm1_v = to_YZX_vm1_v(nqubits,v_site)
        qc.append(basis_qc_YzX_vm1_v,[i for i in range(nqubits)])

    if which_basis == 3:

        basis_qc_YzX_vp1_v = to_YZX_vp1_v(nqubits,v_site)
        qc.append(basis_qc_YzX_vp1_v,[i for i in range(nqubits)])

    return qc


def construct_J1(nqubits,vector,which_qc, v_site):
    #initialize the circuit
    qc = create_init_circuit(nqubits,vector)
    #apply J1:
    qc = apply_J1(nqubits,qc,v_site,which_qc)
    qc = transform_basis(nqubits,qc,v_site,which_qc)
    qc.measure_all()
    return qc


def run_job(nqubits,vector,which_qc,v_site):

    qc = construct_J1(nqubits,vector,which_qc, v_site)
    job = qis.execute(qc, backend=qis.Aer.get_backend('aer_simulator'), shots=10**3,seed_simulator=100)
    counts=job.result().get_counts()
    if (which_qc != 1 and which_qc != 0 and which_qc !=2 and which_qc !=3 ):
        counts = {}

    return counts


