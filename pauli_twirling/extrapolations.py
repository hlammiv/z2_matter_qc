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

import qiskit"""


"""

This code assumes the master data arrangement to be an observable that is formatted like a
dictionary contained like so:

data  =    {    [key]   :  [ step1[<O>_q0,<O>_q1,<O>_q2],  step2[<O>_q0,<O>_q1,<O>_q2], ... ] }   
             ---------------------             --------------------              -----
            |Multiplicity of cnots|     ->    |trotter steps (avgd)|     ->     |qubit|
             ---------------------             --------------------              -----

 For example, 
 data[1] gives all data with 1:1 cnots in our expected circuit, whereas
 data[3] gives all data with 3:1 cnots in our expected circuit

 data[1][0] gives the observable of each qubit at the first trotter step and
 data[1][2][1] gives the observable measured on qubit 1 for the third trotter step

 data[5][2][1] also gives the observable measured on qubit 1 for the third trotter step,
               but the cnots have a multiplicity of 5


"""


import numpy
import gvar as gv

def noiserescale(target_observable,estimation_observable):
    """
    Uses a noise-estimation circuit to rescale the observable of our target circuit by 
    some factor related to the noise from cnots.

    The estimation circuit is similar in structure to the target circuit we are interested
    in simulating, except there is no initialization of a state (all qubits in the default
    | 0 > state) and all single-qubit rotational gates are made trivial (done by taking 
    epsilon->0).

    For example, if the observable being measured is < Z >, then for a qubit in state | 0 >
    we expect the result to be +1. If instead our measurement of the estimation circuit
    gives 0.97, then we can assume the noise in both the estimation and target circuit should 
    have reduced the signal by ~3%.


    Parameters
    ----------
    target_observable: list of GVars
        The observable at each trotter step of the target circuit
    etimation_observable : list of GVars
        The observable at each trotter step of the estimation circuit
    
    Returns
    -------
    rescaled_observable : list of GVars
        The observable at each trotter step, rescaled
    
    """

    rescaled_observable = target_observable/estimation_observable

    return rescaled_observable


def ZNEfits(observable, NCX, nqubits, nsteps):
    """
    The Richardson Extrapolation (or Zero-Noise Extrapolation (ZNE)) assumes the 
    coherent errors from cnots in a circuit cause an exponential decrease in signal 
    as a function of 'Richardson Level' or 'r'. The Richardson Level is related to the 
    multiplicity of cnots (in our case, r = [1,3,5,7,...]).

                        The signal decay function is described as:

                                    A(r) = A0 * exp[B * r] 

    where A0 is the zero-noise estimation and B is a factor related to the noise model. 
    Since we are interested in the zero noise limit (r=0) only parameter A0 is needed for
    this extrapolation. 


    Parameters
    ----------
    observable: Dictionary containing lists of GVars
        For each multiplicity of cnots (keys), the observable at each trotter step (values)
    NCX : list
        The multiplicities of cnots for each step of the ZNE
    nqubits: int
        The total number of qubits
    nsteps: int
        The total number of trotter steps
    
    Returns
    -------
    fits : Dictionary containing lists of GVars
        For each qubit(key), the trotter steps with their zero-noise estimations A0 (values)
    
    """

    #Create the exponential function and its prior parameters
    exponential = lambda x, p: [p['A'] * gv.exp(-p['B'] * x[i]) for i in range(len(x))]

    #This prior was chosen to keep the fits within what is physically possible
    prior = {'f(A)': gv.BufferDict.uniform('f', 1, -1), 'log(B)': gv.gvar(0, 2000)}

    #Initiate the dictonary of fits
    fits = {}

    #Loop over each qubit, and perform the fit for each trotter step
    for qubit in range(nqubits):
        fits['q'+str(qubit)] = [lsqfit.nonlinear_fit(data=(NCX,[observable[j][step][qubit] for j in NCX]), 
                                prior=prior, fcn=exponential, debug=True).p['A'] for step in range(nsteps)]

    return fits

    



