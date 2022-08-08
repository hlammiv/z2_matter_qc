# -*- coding: utf-8 -*-
'''
Copyright August 1, 2022
Authors: Norman Hogan
Created: August 1, 2022
Last edited: August 1, 2022 by Norman Hogan
'''


import qiskit
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,execute
from qiskit.compiler import transpile
import numpy as np
from numpy import pi as pi
import random
from qiskit import IBMQ
from qiskit.circuit.library import XGate, RZGate
from qiskit.transpiler import PassManager, InstructionDurations
from qiskit.transpiler.passes import ALAPSchedule, DynamicalDecoupling

"""
This code builds the Pass Manager (PM) for dynamical decoupling (DD). Once the PM is 
created, the input circuits are then ran through it to insert DD. Three DD sequences are 
supported for this function: CPMG, XY4, and EDD [elaborated in arXiv:2207.03670v1]. 

*   CPMG works best as first order protection for decoupling states near the equator of the 
    Bloch sphere (|+> or |->) with uncertainty from the assumption of perfect pulses. 
*   XY4 is universal first order protection, but also with uncertainty from the assumption of 
    perfect pulses. 
*   EDD is the most universal first order protection using only pi pulses with more accurate 
    assumptions of an imperfect pulse with finite width, but is a long sequence that is not 
    compatible with every quantum computer (see table below). 

The sequences are built like so:

            CPMG:    f/2 - X - f - X - f/2
            XY4:     Y - f - X - f - Y - f - X - f
            EDD:     X - f - Y - f - X - f - Y - f - Y - f - X - f - Y - f - X - f

where f represents a constant delay period between each gate in the sequence, which is 
automatically determined by the PM. This function assumes that the circuits passed through 
it are already transpiled, so the circuits returned by the function are then ready to be run 
on the quantum computer.

NOTE: Not all quantum computers can support XY4 or EDD sequences. If a CNOT operation 
time of a given machine is around 280 nanoseconds or less, it will likely be unable to 
implement the EDD sequence. For some machines, if a delay operation time is not a multiple 
of 16, the job will not run.
Here are the different machines and their compatibility with implementing each sequence:

                  CPMG     XY4      EDD  
    ------------------------------------            
    ibmq_lima:    Yes      Yes      Yes |
    ibmq_belem:   Yes      Yes      Yes |
    ibmq_quito:   Yes      Yes      No  | EDD incompatible operations: cx(1,0),cx(3,4)
    ibmq_manila:  Yes      Yes      Yes |
    ibmq_jakarta: Yes      Yes      No  | EDD incompatible operations: cx(1,0),cx(1,2),cx(3,5),cx(5,6)
    ibm_oslo:     Yes      No       No  | EDD incompatible operations: cx(1,2),cx(3,5),cx(5,6)
    ibm_nairobi:  Yes      Yes      No  | EDD incompatible operations: cx(1,0),cx(3,5),cx(5,4)
    ibm_lagos:    Yes      Yes      No  | EDD incompatible operations: cx(1,2),cx(3,5),cx(5,6)
    ibm_perth:    Yes      Yes      No  | EDD incompatible operations: cx(1,2),cx(3,5),cx(5,6)

"""
def DD_PassManager(circuits: list, DD_sequence: str, providerstr: list, backendstr: str):

    """
    Parameters
    ----------
    circuits: list 
        The transpiled circuts to be run through the DD PM
    DD_sequence: string ("CPMG","XY4","EDD" are the only choices)
        The DD sequence to be implemented during CNOTs
    providerstr: list of 3 strings
        The 'hub', 'group', and 'project' attributes to the backend's provider 
    backendstr: string
        The backend for the quantum computer in which the circuits were transpiled for and 
        will be run on.
    
    Returns
    -------
    DD_circuits: list
        The circuits with DD implemented, ready to be run on the specified backend
    
    """

    #Initialize DD circuits:
    DD_circuits = []

    #Initialize the DD sequences (XGate(), RZGate(pi) == Y in terms of basis gates):
    CPMG = [XGate(), XGate()]
    XY4 = [XGate(), RZGate(pi), XGate(), XGate(), RZGate(pi), XGate()]
    EDD = [XGate(), XGate(), RZGate(pi), XGate(), XGate(), RZGate(pi), 
           XGate(), RZGate(pi), XGate(), XGate(), RZGate(pi), XGate()]

    #Initialize delay spacing for each sequence (Must sum up to 1):
    f_CPMG = [1/4,1/2,1/4]
    f_XY4 = [0,0,1/4,1/4,0,1/4,1/4]
    f_EDD = [0,1/8,0,1/8,1/8,0,1/8,0,1/8,1/8,0,1/8,1/8]

    #Get the provider for the backend:
    provider = IBMQ.get_provider(hub=providerstr[0], 
                                group=providerstr[1], 
                                project=providerstr[2])

    #Get the backend:
    backend = provider.get_backend(backendstr)

    #Save the basis gate durations (i.e. single qubit gate and CNOT time) from the backend
    instruct = InstructionDurations.from_backend(backend)

    #Implement CPMG DD sequence if specified:
    if DD_sequence == 'CPMG':

        #Build the pass manager with CPMG DD:
        pm = PassManager([ALAPSchedule(instruct), 
                   DynamicalDecoupling(instruct, CPMG, spacing=f_CPMG)])

    #Implement XY4 DD sequence if specified:
    elif DD_sequence == 'XY4':

        #Build the pass manager with XY4 DD:
        pm = PassManager([ALAPSchedule(instruct), 
                   DynamicalDecoupling(instruct, XY4, spacing=f_XY4)])

    #Implement EDD DD sequence if specified:
    elif DD_sequence == 'EDD':

        #Build the pass manager with EDD DD:
        pm = PassManager([ALAPSchedule(instruct), 
                   DynamicalDecoupling(instruct, EDD, spacing=f_EDD)])



    #Run the PM on the transpiled circuits
    DD_circuits = pm.run(circuits)


    return DD_circuits

    