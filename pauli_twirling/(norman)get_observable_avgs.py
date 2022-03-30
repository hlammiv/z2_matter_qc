import gvar as gv
import numpy as np

# import txt file with measure and keys dictionaries?


"""
Immediate results of job running function assumed to be in the form:

job_results = [list of size richardson_level of
                    [list of length ncircuits of 
                        [lists of length ntrotter of dicts {'state': count, ...}]

where
job_results is a list of length of number of Richardson levels
job_results[0] gives a list of trotter evolutions for multiple circuit runs to be averaged over with richardson_level = 1
job_results[0][0] gives an ntrotter-length list of dictionaries of counts for the first circuit
job_results[0][0][0] gives the dictionary of counts for the first trotter step of the first circuit
...

"""

# The observable has two expecation values for a measurement of either 0 or 1. for example:

#   < Z > = { 0: 1, 1: -1}
#   < S_z > = { 0: 0.5*hbar, 1: -0.5*hbar}
#   etc......

def avg_observable_per_time_step(job_results, ntrotter, nqubits, shots, observable):
    """

    Parameters
    ----------------------------
    job_results = list
        The counts from running the trotterized circuit, with ncircuit amount of circuits per trotter step 
        for each richardson_level
    ntrotter = int
        The number of trotter steps
    nqubits = int
        The number of qubits
    shots = int
        The number of shots run, from initial parameters
    observable = dict
        The observable's expectation value for either the 0 or 1 state


    Returns
    ----------------------------
    obs_list_trotter = list 
        The observable averaged over the amount of circuits run for each trotter step

    """
    all_obs = np.zeros((len(job_results), ntrotter, nqubits))


    #Loop over ncircuits:
    for i in range(len(job_results)):

        one_evolution = job_results[i]
        # one evolution is ntrotter amount of trotterized steps
        # -----------------------------------------------------

        #loop over ntrotter
        for j in range(ntrotter):
            counts_dict = one_evolution[j]

            # add any missing keys and counts in count dictionary to keep state and operator matrices consistently sized:
            # -----------------------------------------------------------------------------------------------------------
            for key in counts_dict.keys():
                for site in range(nqubits):
                    if key[site] == '0':
                        all_obs[i, j, site] += observable[0] * counts_dict[key] / shots
                    else:
                        all_obs[i, j, site] += observable[1] * counts_dict[key] / shots
            # adds measurement of desired observable to the ncircuit row and ntrotter column
            # ------------------------------------------------------------------------------

    obs_list_trotter = gv.gvar(np.mean(all_obs, axis=0), np.std(all_obs, axis=0)/np.sqrt(len(job_results)))
    # averages down columns, which is the observable for all the circuits run for that trotter step
    # ---------------------------------------------------------------------------------------------

    return obs_list_trotter


def get_all_observables(ntrotter: int, richardson_level: int, nqubits:int, job_results_list, shots, measure):
    """
    Objective of function is to collect and organize all data into format easily used by extrapolation.py script
    Calls avg_observable_per_time_step function

    Parameters
    ---------------------------
    ntrotter = int
        The amount of trotter steps
    richardson_level = int
        MAX richardson_level run; if there is data for several richardson levels, this function takes the highest level
        and creates a dictionary with a key for every possible amount of CNOTs up to that level; job_results_list should
        correspond to amount of levels
    nqubits: int
        The number of qubits
    job_results_list = list
        Job results in a in a list containing results from each richardson level run
    shots = int
        Number of shots from job
    measure = dict
        Exected measurement of an observable for a qubit in either the 0 or 1 state


    Returns
    --------------------------
    data in the format:

    data  =    {    [key]   :  [ step1[<O>_q0,<O>_q1,<O>_q2],  step2[<O>_q0,<O>_q1,<O>_q2], ... ] }
             ---------------------             --------------------              -----
            |Multiplicity of cnots|     ->    |trotter steps (avgd)|     ->     |qubit|
             ---------------------             --------------------              -----

    """

    data = {}

    for level in range(1,richardson_level+1):
        ncx = level * 2 - 1
        
        data[ncx] = avg_observable_per_time_step(job_results_list[level-1], ntrotter, 
                                                    nqubits, shots, measure)

    return data
