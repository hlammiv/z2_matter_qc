# -*- coding: utf-8 -*-
"""
Created on Tue Tue May 10 2022
Last edited on Tue May 10 2022

@authors: Florian Herren, Ruth Van de Water
          
Prepare gauge-invariant basis states for Z2 gauge theory simulation
Still need to write general function for arbitrary number of states/qubits
"""

## The functions below all return `qiskit.QuantumCircuit` objects.

def prepare_zbasis_state(nq : int, state : list):
    qc = QuantumCircuit(nq,nq)
    for i,bit in enumerate(state):
        if (bit == 1):
            qc.x(nq-i-1)

    return qc

def prepare_states_nsites2(nq : int):
    if nq != 3: 
        print('*** WARNING: number of qubits is 3 for 2-site system ***')
        return
    
    # create dictionary to hold gauge-invariant, Q_net=0 states
    state_dict = {}
    
    # prepare empty (=vacuum) state and add to dictionary
    qc_vac = QuantumCircuit(nq,nq)
    qc_vac.h(1)
    state_dict['vacuum'] = qc_vac
    
    # prepare excited state with one e+ and one e-
    qc_mes = QuantumCircuit(nq,nq)
    for n in range(nq):
        qc_mes.x(n)
    qc_mes.h(1)
    state_dict['meson'] = qc_mes
    
    return state_dict

def prepare_states_nsites4(nq : int):
    if nq != 7: 
        print('*** WARNING: number of qubits is 7 for 4-site system ***')
        return
    
    # create dictionary to hold gauge-invariant, Q_net=0 states
    state_dict = {}
    
    # index of middle qubit
    midpoint = int((nq-1)/2)
    
    # prepare empty (=vacuum) state and add to dictionary
    qc_vac = QuantumCircuit(nq,nq)
    for n in range(1,nq,2):
        qc_vac.h(n)
    state_dict['vacuum'] = qc_vac
    
    # prepare excited state with one e+ and one e-
    # in the middle of lattice
    qc_mes = QuantumCircuit(nq,nq)
    for n in range(midpoint-1,midpoint+2):
        qc_mes.x(n)
    for n in range(1,nq,2):
        qc_mes.h(n)
    state_dict['meson'] = qc_mes
    
    # prepare excited state with e+ and one e-
    # in first two sites
    qc_mes_start = QuantumCircuit(nq,nq)
    for n in range(3):
        qc_mes_start.x(n)
    for n in range(1,nq,2):
        qc_mes_start.h(n)
    state_dict['meson-start'] = qc_mes_start
    
    # prepare excited state with e+ and one e-
    # in last two sites
    qc_mes_end = QuantumCircuit(nq,nq)
    for n in range(nq-1,nq-4,-1):
        qc_mes_end.x(n)
    for n in range(1,nq,2):
        qc_mes_end.h(n)
    state_dict['meson-end'] = qc_mes_end
    
    # prepare excited state with e- in first site
    # and e+ in last site
    qc_mes_stretch = QuantumCircuit(nq,nq)
    qc_mes_stretch.x(0)
    qc_mes_stretch.x(nq-1)
    for n in range(1,nq,2):
        qc_mes_stretch.x(n)
        qc_mes_stretch.h(n)
    state_dict['meson-stretch'] = qc_mes_stretch
    
    # prepare excited state with all sites filled
    qc_2mes = QuantumCircuit(nq,nq)
    for n in range(0,nq,2):
        qc_2mes.x(n)
    for n in range(1,nq,4):
        qc_2mes.x(n)
    for n in range(1,nq,2):
        qc_2mes.h(n)
    state_dict['2meson'] = qc_2mes
    
    return state_dict