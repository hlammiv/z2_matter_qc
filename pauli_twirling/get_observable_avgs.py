import gvar as gv
import numpy as np

# import txt file with measure and keys dictionaries?


"""
Immediate results of job running function assumed to be in the form:

job_results = [list of length ncircuits of [lists of length ntrotter of {'state': count, ...}, {'state': count, ...},...],[{'state': count, ...},...]]

where
job_results is a list of length of number of circuits per trotter step
job_results[0] gives a list of length ntrotter, with each entry being a dictionary of counts for that trotter step
job_results[0][0] gives a dictionary of counts for the first trotter step
job_results[1][0] also gives a dictionary of counts for the first trotter step, for the second time running the full trotterized circuit
...
job_results[ncircuits][0] gives the final dictionary of counts for the first trotter step after running the full trotterized circuit ncircuit times

code is also flexible for job_results = [list of length ntrotter of [list of length ncircuits of {dicts}],[{dicts}]] with a few changes

"""

measure_2sites = {'Z0': [1, 1, 1, 1, -1, -1, -1, -1], 'Z1': [1, 1, 1, 1, 1, 1, -1, -1],
                  'Z2': [1, -1, 1, -1, 1, -1, 1, -1]}
keys_2sites = ['110 000', '000 000', '100 000', '010 000', '001 000', '101 000', '111 000', '011 000']


def avg_observable_per_time_step(ntrotter, job_results, shots, keys, observable):
    """

    Parameters
    ----------------------------
    job_results = the counts from running the trotterized circuit, with ncircuit amount of circuits per trotter step
    shots = number of shots run, from initial parameters
    keys = list of possible states for either the 2 site or 4 site case
    observable = the observable operator in the measure dictionary that we wish to find


    Returns
    ----------------------------
    obs_list_trotter = list of avg observable for each trotter step, averaged over the amount of circuits run each trotter step

    """

    all_obs = np.zeros((len(job_results), ntrotter))

    for i in range(len(job_results)):

        one_evolution = job_results[i]
        # one evolution is ntrotter amount of trotterized steps
        # -----------------------------------------------------

        for j in range(len(one_evolution)):
            counts_dict = one_evolution[j]

            # add any missing keys and counts in count dictionary to keep state and operator matrices consistently sized:
            # -----------------------------------------------------------------------------------------------------------
            for k in range(len(keys)):
                if keys[k] not in counts_dict:
                    counts_dict[keys[k]] = 0

            counts = list(counts_dict.values())
            state = np.diag(counts[::-1])
            all_obs[i, j] = (np.trace(state @ np.diag(observable))) / shots
            # adds measurement of desired observable to the ncircuit row and ntrotter column
            # ------------------------------------------------------------------------------

    obs_list_trotter = np.average(all_obs, axis=0)
    # averages down columns, which is the observable for all the circuits run for that trotter step
    # ---------------------------------------------------------------------------------------------

    err = np.sqrt((1 - obs_list_trotter ** 2) / shots)
    for z in obs_list_trotter:
        obs_list_trotter[z] = gv.gvar(z, (err[i] for i in range(len(err))))
        # combines value and error into single gvar list
        # ----------------------------------------------

    return obs_list_trotter


def get_all_observables(ntrotter: int, richardson_level: int, job_results_list, shots=shots, measure=measure_2sites,
                        keys=keys_2sites):
    """
    Objective of function is to collect and organize all data into format easily used by extrapolation script
    Calls avg_observable_per_time_step function

    Parameters
    ---------------------------
    ntrotter = amount of trotter steps
    richardson_level = MAX richardson_level run; if there is data for several richardson levels, this function takes the highest level
                        and creates a dictionary with a key for every possible amount of CNOTs up to that level; job_results_list should
                        correspond to amount of levels
    job_results_list = job results in format discussed above, in a list containing results from each richardson level run
    shots = shots from job
    measure = dictionary of operators for 2 site or 4 site observables
    keys = list of possible states for 2 site or 4 site results


    Returns
    --------------------------
    data in the format:

    data  =    {    [key]   :  [ step1[<O>_q0,<O>_q1,<O>_q2],  step2[<O>_q0,<O>_q1,<O>_q2], ... ] }
             ---------------------             --------------------              -----
            |Multiplicity of cnots|     ->    |trotter steps (avgd)|     ->     |qubit|
             ---------------------             --------------------              -----

    """

    data = {}

    for level in range(1, richardson_level):
        ncx = level * 2 - 1

        qubit_obs_list = []
        qubit_each_step = []
        total_steps = []
        qubit_obs_list.append(avg_observable_per_time_step(ntrotter, job_results[level - 1], shots, keys,
                                                           (measure[key] for key in measure)))
        for step in range(ntrotter):
            for qubit in qubit_obs_list:
                qubits_each_step.append(qubit[step])
            total_steps.append(qubits_each_step)
        data[ncx] = total_steps

    return data



