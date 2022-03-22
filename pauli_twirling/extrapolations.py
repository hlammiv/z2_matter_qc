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

import numpy
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
 data[1][2][1] gives the observable measured on qubit 1 for the third trottter step


"""


import numpy

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


def ZNEfits(observable, NCX):
    """
    Uses exponential fits to determine the zero-noise limit


    Parameters
    ----------
    observable: Dictionary containing lists of GVars
        For each multiplicity of cnots (keys), the observable at each trotter step (values)
    NCX : list
        The multiplicities of cnots for each step of the zero-noise extrapolation (ZNE)
    
    Returns
    -------
    fits : list of lsqfit.nonlinear_fit
        The fits and their parameters for each multiplicity of cnots
    
    """

    



