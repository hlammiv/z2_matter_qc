# -*- coding: utf-8 -*-
"""
Created on Tue Tue May 10 2022
Last edited on Tue May 10 2022

@authors: Ruth Van de Water
          
Time evolve initial states (without noise) and plot resulting final states
Will expand comments later
"""

#############################################
### Process `qiskit.statevector` objects. ###
#############################################

## In the functions below, quantum_state is always a `qiskit.QuantumCircuit` objects.

# Round coefficient of each z-basis to fixed decimals (default=2)
def round_statevector(quantum_state, ndecimals: Optional[int] = 2):
    for key in quantum_state:
        quantum_state[key] = np.round(quantum_state[key],ndecimals)
    return quantum_state

# Replace coefficient of each z-basis state by absolute value
def abs_statevector(quantum_state):
    for key in quantum_state:
        quantum_state[key] = np.abs(quantum_state[key])
    return quantum_state

# Get unique z-basis states from list of statevectors
def get_unique_states(statevecs):
    state_list = []
    for state in statevecs:
        state_list = state_list+[k for k in state]
    state_list = list(set(state_list))
    return state_list

# Remove all dictionary entries for basis state if coefficient 
# is below cutoff value for all statevectors in list
def clean_statevectors(statevectors, cutoff = 1e-5):
    
    state_list = get_unique_states(statevectors) 
    
    # run over all unique basis states
    for s in state_list:
        c_list = []
    
        # run over all statevectors appending to c_list
        for i,statevec in enumerate(statevectors):
            
            # if coefficient > cutoff append False (nonzero value ==> cannot delete key)
            #     " " "      < cutoff, append True (negligible value ==> can delete key)
            try:
                coeff = statevec[s]
                if (coeff < cutoff):
                    c_list.append(True)
                else:
                    c_list.append(False)
            except:
                continue

        # if all c_list entries are True, delete key from all statevectors
        if all(c_list):
            for j,statevec in enumerate(statevectors):
                statevec.pop(s,None)
            
    return

######################################################
### Evolve single state for multiple Trotter steps ###
######################################################

def time_evolve_singlestate(state : dict, name : str, eps_list : list, 
                            mass : float, ntrotter: int, outfile : str, 
                            twirl = False, print_state = True, figsize = (10,5)):
    """
    Time evolve quantum circuit by exp[-i H epsilon]^{ntrotter}, 
    where H is the 1+1d Z2 Hamiltonian, for multiple step sizes epsilon. 
    Then save histogram of final states to png image file.

    Parameters
    ----------
    state : qiskit.QuantumCircuit 
        the initial state vector
    name : str
        state-vector label
    epsilon : list
        entries are Trotter step sizes in (1/a) units 
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation        
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
        
    Output
    ----------
    outfile : str
        save histogram of final states as image file
        named outfile.png
    """
    
    statevectors = []
    legend = []
    probabilities_zbasis = []
    print_state = True

    print('initial\nstate    statevector\n---------------------------------------------')
    # loop over list of epsilon values
    for epsilon in eps_list:
    
        # reset initial state each time
        initial_state = state.copy()
        psi_dict = Statevector(initial_state).to_dict() # human-readable dictionary
        
        # compute and print net charge
        Qnet = list(set([net_charge(key) for key in psi_dict]))
        if (len(Qnet) > 1):
            print('*** WARNING: state does not have definite net charge ***')
        if (print_state):
            print('vacuum  ', round_statevector(psi_dict))
            print(f'Qnet = {Qnet[0]}')
            print_state = False
            
        # infer number of lattice sites from state vector
        nqubits = initial_state.num_qubits
        nsites = nqubits2nsites(nqubits)
        
        # time evolve initial_state and record final state
        final_state, probabilities, total_prob = trotter_evolve(initial_state, nsites,
                                                                epsilon, mass, ntrotter, twirl=twirl)
        statevectors.append(final_state)
        probabilities_zbasis.append(probabilities)  # useful for checking unitarity
       
        # histogram labels
        legend.append("eps.="+str(round(epsilon,2))+"; sum of probabilities="+str(round(total_prob,2)))

    # remove dictionary entries for basis state if coefficient is zero for all statevectors in list
    clean_statevectors(probabilities_zbasis)      
    state_list = get_unique_states(probabilities_zbasis) 
    print('\nall states present after time evolution '+str(state_list))

    plot_histogram(probabilities_zbasis,title='State vector after evolving '+name+' state by single Trotter step',legend=legend,bar_labels=False,color=colors[:len(statevectors)],figsize=figsize,filename=outfile)
    
######################################################
### Evolve multiple states for single Trotter step ###
######################################################

def time_evolve_multistates(states : dict, epsilon : float, mass : float, 
                            ntrotter: int, outfile : str, twirl = False, 
                            print_state = True, figsize = (10,5)):
    """
    Time evolve initial states by exp[-i H epsilon]^{ntrotter}, 
    where H is the 1+1d Z2 Hamiltonian,
    and save histogram of final states to png image file.

    Parameters
    ----------
    states : dict 
        dictionary values are the initial states (qiskit.QuantumCircuit objects)
        and keys are the states' labels (strings)
    epsilon : float
        the Trotter step size in (1/a) units
    mass : float
        mass of the fermion in units of a (a=lattice spacing)
    ntrotter : int
        the number of Trotter steps for the simulation        
    twirl : boolean (optional)
        whether to implement this circuit with randomized
        compiling. The default is False.
        
    Output
    ----------
    outfile : str
        save histogram of final states as image file
        named outfile.png
    """    
    
    statevectors_zbasis = []
    probabilities_zbasis = []
    legend = []
    
    # loop over dictionary of initial states
    for name in states:
        initial_state = states[name] 
        psi_dict = Statevector(initial_state).to_dict() # human-readable dictionary  
        
        # compute and print net charge
        Qnet = [net_charge(key) for key in psi_dict][0] 
        print(state, psi_dict, Qnet)
        
        # infer number of lattice sites from state vector
        nqubits = initial_state.num_qubits
        nsites = nqubits2nsites(nqubits)
        
        # time evolve initial_state and record final state
        final_state, probabilities, total_prob = trotter_evolve(initial_state, nsites, epsilon, 
                                                                mass, ntrotter, twirl=twirl)
        statevectors_zbasis.append(final_state)     
        probabilities_zbasis.append(probabilities)  # useful for checking unitarity
    
        # histogram labels
        legend.append("psi0="+name+"; sum of probabilities="+str(round(total_prob,2)))
    
    plot_histogram(probabilities_zbasis,title="State vector after single Trotter step with epsilon="+str(epsilon),legend=legend,bar_labels=False,color=colors[:len(statevectors_zbasis)],figsize=figsize,filename=outfile)
    
    

    

