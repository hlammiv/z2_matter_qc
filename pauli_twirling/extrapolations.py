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
    Uses exponential fits to determine the zero-noise limit


    Parameters
    ----------
    observable: Dictionary containing lists of GVars
        For each multiplicity of cnots (keys), the observable at each trotter step (values)
    NCX : list
        The multiplicities of cnots for each step of the zero-noise extrapolation (ZNE)
    nqubits: int
        The total number of qubits
    nsteps: int
        The total number of trotter steps
    
    Returns
    -------
    fits : Dictionary of lsqfit.nonlinear_fit
        The fits and their parameters for each qubit(key)
    
    """

    #Create the exponential function and its prior parameters
    exponential = lambda x, p: [p['A'] * gv.exp(-p['B'] * x[i]) for i in range(len(x))]
    prior = {'f(A)': gv.BufferDict.uniform('f', 1, -1), 'log(B)': gv.gvar(0, 2000)}

    #Initiate the dictonary of fits
    fits = {}

    #Loop over each trotter step
    for step in range(nsteps):
        #Loop over each qubit
        for qubit in range(nqubits):
            fits['q'+str(qubit)] = lsqfit.nonlinear_fit(data=(NCX,[observable[j][step][qubit] for j in NCX]), 
                                    prior=prior, fcn=exponential, debug=True)

    return fits

    



