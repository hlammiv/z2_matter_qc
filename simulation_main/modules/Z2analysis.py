# -*- coding: utf-8 -*-
"""
Created on Tue Fri Mar 25 2022
Last edited on Tue May 10 2022

@authors: Clement Charles, Florian Herren, Sara Starecheski, Ruth Van de Water
          
Analyze output of Z2 gauge theory simulation
Will expand comments later
"""

import numpy as np
from typing import Optional

# calculate number of qubits given number of lattice sites
def nsites2nqubits(ns : int):
    nq = 2*ns-1
    return int(nq)

# calculate number of lattice sites given number of qubits
def nqubits2nsites(nq : int):
    ns = (nq+1)/2
    return int(ns)

# calculate net charge of z-basis state given as string of 0s and 1s
def net_charge(state : str):
    # initialize to zero
    Qnet = 0
    
    # loop through even sites
    for s in state[0::4]:
        if (s == '1'):
            Qnet += 1
            
    for s in state[2::4]:
        if (s == '1'):
            Qnet += -1
        
    return Qnet

def get_particle_number(counts, nshots: int, n: int, nboots: Optional[int] = None):
    """
    Calculates the percentage of shots for which
    the nth qubit, q[n], has a particle (electron or positron) in it. 
    The site is empty if q[n] = |0> and occupied if q[n] = |1>.
    Parameters
    ----------
    counts : qiskit.result.counts.Counts
        
    nshots : int 
        number of shots in quantum simulation
    n : int
        index of the qubit for which we want the mean occupancy
    nboots : int (optional)
        number of bootstrap ensembles for calculating uncertainty in the mean.
        default is same as nshots.
    Returns
    ----------
    mean : float
        the average occupancy of q[n]
    err : float
        the bootstrap error in the mean
    """

    keys = list(counts.keys())
    values = list(counts.values()) 
    total_counts = sum(values) #total counts!=nshots for mitigated results
    
    # first calculate mean
    mean = 0
    for s in counts:
        p = s[n]
        if p == '1':
             mean = mean + (counts[s]/total_counts)
                
    # next calculate bootstrap error
    if not(nboots):
        nboots = nshots
    prob = [counts[k]/total_counts for k in keys] #use total_counts in denominator so that probabilities sum to 1
    means = []
    for b in range(nboots):
        m = 0
        samples = np.random.choice(keys, size=nshots, p=prob)
        for s in samples:
            p = s[n]
            if p == '1':
                m = m + (1/nshots)
        means.append(m)
    err = np.std(means)
    
    return mean, err
