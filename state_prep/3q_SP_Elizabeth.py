#import statements needed as always 
import numpy as np
from qiskit import *

# this is just required overall
#define Pauli Matricies 
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])
Had = 1./np.sqrt(2)*np.array([[1,1],[1,-1]])

#helper functions 
#Hamiltonian index function
def index(e1,p1,y1,e2,p2,y2,m,g):
    ''' This function takes in arguements e1,p1,y1,e2,p2,y2,m,g. 
    e1,p1,y1,e2,p2,y2 correspond to the 2 dimensional vectors relating 
    to the 0 and 1 states. m,g are the paramaters corresponding to the
    size of the mass and coupling terms of the hamiltionian. This function
    is a helper function inside the 3q state prep.'''
    first  = (m/2)*np.dot(np.dot(e1,Z),e2)
    second = (-m/2)*np.dot(np.dot(p1,Z),p2)
    third = 1*np.dot(np.dot(y1,X),y2)/2
    forth = (g/4)*np.dot(np.dot(y1,Z),y2)*np.dot(np.dot(e1,X),e2)*np.dot(np.dot(p1,X),p2)
    fifth = (g/4)*np.dot(np.dot(y1,Z),y2)*np.dot(np.dot(e1,Y),e2)*np.dot(np.dot(p1,Y),p2)
    return first+second+third+forth+fifth


#state array function
def state_array(q):
    '''This function takes in the number of qubits, q
    and outputs an array containing the qubit ordering in 2d states. 
    This function is a helper function inside state prep.'''
    #define the number of qubits, currently we are using 3
    qubits = q
    #define the total number of possibilities 
    length = 2**qubits
    #deine a list to hold the bianary representation
    stat = []
    #define a function to give that bianary representation from an integer
    getbinary = lambda x, n: format(x, 'b').zfill(n)
    #loop through the integers to get the numbers in bianary
    #aka, the possible states 
    for i in range(length):
        stat.append(getbinary(i, qubits))

    #define a list to hold the vector representation of the states 
    state = []
    #loop over all states 
    for i in range(len(stat)):
        #break up states by each value 
        fornow = list(stat[i])
        #define a list to hold the vector form 
        fornow2 = []
        #loop over string values, if add vector 1, if 0, add vector 0
        for j in range(len(fornow)):
            if fornow[j] == '1':
                fornow2.append([0,1])
            if fornow[j] == '0':
                fornow2.append([1,0])
        #add the list of vectors to the full list with all states 
        state.append(fornow2)
    #convert the list of lists of lists into a numpy array 
    state = np.array(state)
    return state


#make Hamiltonian
def Hamiltonian_maker(state,m,g):
    ''' this function takes in the coupling terms and the array
    of states and outputs the hamiltionian. This function is a 
    helper function inside 3q state prep'''
#loop over all possible states
    Hamil= []
    #this is also for the specified hamiltonian, may need to change depending on new hamiltonian 
    for i in range(len(state)):
        e1 = np.array(state[i])[0,:]
        y1 = np.array(state[i])[1,:]
        p1 = np.array(state[i])[2,:]
        Ham = []
    #for each state loop over all possible states to get all terms
        for j in range(len(state)):
            e2 = np.array(state[j])[0,:]
            y2 = np.array(state[j])[1,:]
            p2 = np.array(state[j])[2,:]
            Ham.append(index(e1,p1,y1,e2,p2,y2,m,g))
        Hamil.append(Ham)
    #convert list into numpy array
    Hamiltonian = np.array(Hamil)
    return Hamiltonian


#sort eigenvalues and vectors - not entierly needed but just put here anyways
def arr(evals,evec):
    ''' this function takes in eignvalues and eigenvectors and arranges them. 
    this function is a helper function in state prep'''
    item = np.argsort(evals)
    vectors = []
    values = []
    for i in range(len(evals)):
        values.append(evals[item[i]])
        vectors.append(evec[item[i]])
    return vectors, values


#gauge invarient matricies 
def Guage0(state):
    ''' this function takes in the state array and defines one of the gauge invarient matricies.
    it is used as a helper function state prep '''
#loop over all possible states
    Gauge= []
    #this is also for the specified hamiltonian, may need to change depending on new hamiltonian 
    for i in range(len(state)):
        e1 = np.array(state[i])[0,:]
        y1 = np.array(state[i])[1,:]
        p1 = np.array(state[i])[2,:]
        G = []
    #for each state loop over all possible states to get all terms
        for j in range(len(state)):
            e2 = np.array(state[j])[0,:]
            y2 = np.array(state[j])[1,:]
            p2 = np.array(state[j])[2,:]
            G.append(y1@X@y2*p1@p2*e1@Z@e2)
        Gauge.append(G)
    #convert list into numpy array
    Final = np.array(Gauge)
    return Final

def Guage2(state):
    ''' this function takes in the state array and defines one of the gauge invarient matricies.
    it is used as a helper function state prep '''
#loop over all possible states
    Gauge= []
    #this is also for the specified hamiltonian, may need to change depending on new hamiltonian 
    for i in range(len(state)):
        e1 = np.array(state[i])[0,:]
        y1 = np.array(state[i])[1,:]
        p1 = np.array(state[i])[2,:]
        G = []
    #for each state loop over all possible states to get all terms
        for j in range(len(state)):
            e2 = np.array(state[j])[0,:]
            y2 = np.array(state[j])[1,:]
            p2 = np.array(state[j])[2,:]
            G.append(y1@X@y2*p1@-Z@p2*e1@e2)
        Gauge.append(G)
    #convert list into numpy array
    Final = np.array(Gauge)
    return Final


#Gramm-Schmidt
def gs(A):
    '''this function applies gramm-schmidt to a given set of row vectors. 
    it is used as a helper function in state prep.'''
    A = 1.0*A
    #delete the all zero columns in the original matrix
    A = np.delete(A, np.argwhere(np.all(A[:,...] < 1e-6, axis=0)), axis=0)
    (m,n) = A.shape
    O = np.zeros((m,n), dtype=A.dtype)
    o_i = 0
    ndel = m-1
    for i in range(m):
        #print(O)
        q = A[i] # i-th row of A
        for j in range(o_i):
            q = q - np.dot(np.conj(O[j]), q) * O[j]
        if (np.linalg.norm(q)/np.linalg.norm(A[i]) < 1e-6):
            O = np.delete(O, obj=ndel, axis=0)
            ndel -= 1
        else:
            q = q / np.linalg.norm(q)
            O[o_i] = q
            o_i += 1
    return O


#Energies function
def efind(psi, H): 
    '''this function finds the eigenenergies from a given set of
    eigenvectors and hamiltonian corresponding to those vectors'''
    #psi is the vectors 
    #H is the hamiltonian
    E = []
    (m,n) = psi.shape
    for i in range(m):
        #E.append((np.dot(np.dot(psi[i,:],H),psi[i,:]))/(np.dot(psi[i,:],psi[i,:])))
        E.append((np.conj(psi[i,:])@(H@psi[i,:]))/np.dot(psi[i,:],psi[i,:]))
    E = np.array(E)
    return E


#change basis
def zero_plus_one(state):
    '''this function defines the zero-plus-one state to change
    the invarient vectors basis. it is a helper function in state prep'''
    #this is also for the specified hamiltonian, may need to change depending on new hamiltonian 
    for i in range(len(state)):
        e = np.array(state[0])[0,:]
        y = np.array(state[0])[1,:]
        p = np.array(state[0])[2,:]
        zpo = np.kron(np.kron(e, Had@y), X@p)
    return zpo

def one_minus_zero(state):
    '''this function defines the one-minus-zero state to change
    the invarient vectors basis. it is a helper function in state prep'''
    #this is also for the specified hamiltonian, may need to change depending on new hamiltonian 
    for i in range(len(state)):
        e = np.array(state[0])[0,:]
        y = np.array(state[0])[1,:]
        p = np.array(state[0])[2,:]
        omz = np.kron(np.kron(X@e, (Had@X)@y), p)
    return omz


#state prep
def three_q_SP(m,g,n):
    '''This function takes in m,g which are the coefficents to the
    Hamiltionian and n which is the invarentvector to calculate the theta from 
    and outputs 2 values: theta and the energy of that gauge 
    invarient vector of the Hamiltonian'''
    
    #get the states
    state = state_array(3)
    #get the gauge invarient matrices 
    G0 = Guage0(state)
    G2 = Guage2(state)
    #get the projectors 
    P0 = (np.eye(8) + G0)/2
    P2 = (np.eye(8) + G2)/2
    P = P0@P2
    P_non = (np.eye(8) - P)
    
    #get the hamiltonian 
    Hamiltonian = Hamiltonian_maker(state,m,g).real
    #print("Hamiltonian with m = ", m, 'g = ', g)
    #print(Hamiltonian)
    
    #calculate eigen values and vectors
    val, vec = np.linalg.eig(Hamiltonian)
    #sort them
    svec, sval = arr(val,vec)
    vecsorted = np.array(svec)
    valsorted = np.array(sval)
    
    #apply prjector to eigenvectors 
    #apply as P@vec because each vector is column of vec
    PV = P@vecsorted
    PV2 = P_non@vecsorted
    
    #apply gramm-schmidt 
    Invarient = gs(np.transpose(PV))
    NonInvarient_all = gs(np.transpose(PV2))

    #get energies just cuz
    E = efind(Invarient, Hamiltonian).real
    #print("The Engergies are: ")
    #print(E)
    
    #get the change of basis vectors 
    V_base = zero_plus_one(state)
    P_base = one_minus_zero(state)
    
    #get alpha and beta just cuz
    alpha = []
    beta = []
    #get theta
    theta = []
    for i in range(len(Invarient)):
        alpha.append(np.real(np.dot(V_base,Invarient[i])))
        beta.append(np.real(np.dot(P_base,Invarient[i])))
        theta.append(2*np.arccos(alpha[i]))

    return theta[n],E[n]

#circuit of state prep 
def three_q_prep_circuit(theta):
    '''This function takes in one value of theta corresponding to one
    gauge invaient eigenvector. It outputs the quantum circuit that will
    be the prepared state.'''
    q = QuantumRegister(3)

    qc = QuantumCircuit(q)

    qc.rx(theta,1)
    qc.h(1)
    qc.cnot(1,0)
    qc.cnot(1,2)
    qc.x(2)
    qc.h(1)
    
    return qc

#aternate option for state prep if many vectors wanted 
def alternate_three_q_SP(m,g):
    '''This function takes in m,g which are the coefficents to the
    Hamiltionian. It outputs theta and the energy of that gauge invarient vectors
    of the Hamiltonian. It is differs from three_q_SP in that it outputs the values for all
    of the invarient vectors not just one prechoosen one.'''
    
    #get the states
    state = state_array(3)
    #get the gauge invarient matrices 
    G0 = Guage0(state)
    G2 = Guage2(state)
    #get the projectors 
    P0 = (np.eye(8) + G0)/2
    P2 = (np.eye(8) + G2)/2
    P = P0@P2
    P_non = (np.eye(8) - P)
    
    #get the hamiltonian 
    Hamiltonian = Hamiltonian_maker(state,m,g).real
    #print("Hamiltonian with m = ", m, 'g = ', g)
    #print(Hamiltonian)
    
    #calculate eigen values and vectors
    val, vec = np.linalg.eig(Hamiltonian)
    #sort them
    svec, sval = arr(val,vec)
    vecsorted = np.array(svec)
    valsorted = np.array(sval)
    
    #apply prjector to eigenvectors 
    #apply as P@vec because each vector is column of vec
    PV = P@vecsorted
    PV2 = P_non@vecsorted
    
    #apply gramm-schmidt 
    Invarient = gs(np.transpose(PV))
    NonInvarient_all = gs(np.transpose(PV2))

    #get energies just cuz
    E = efind(Invarient, Hamiltonian).real
    #print("The Engergies are: ")
    #print(E)
    
    #get the change of basis vectors 
    V_base = zero_plus_one(state)
    P_base = one_minus_zero(state)
    
    #get alpha and beta just cuz
    alpha = []
    beta = []
    #get theta
    theta = []
    for i in range(len(Invarient)):
        alpha.append(np.real(np.dot(V_base,Invarient[i])))
        beta.append(np.real(np.dot(P_base,Invarient[i])))
        theta.append(2*np.arccos(alpha[i]))

    return theta,E