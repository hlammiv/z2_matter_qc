# -*- coding: utf-8 -*-
"""
Created on Tue Fri Jun 23 2023
Last edited on Fri Jun 23 2023
          
Prepare gauge-invariant basis states for Z2 gauge theory simulation with 3 Qubits
Code written by Immanol Benitez and edited by Florian Herren
"""
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister


#Compute the projector
def projall():
    id2 = np.identity(2)
    ox = np.array([[0,1],[1,0]])
    oy = np.array([[0+0j,0-1j],[0+1j,0+0j]])
    oz = np.array([[1,0],[0,-1]])
    h_gate = 1/np.sqrt(2)*(ox+oz)
    projector_z0_x1 = (1/2)*np.kron(oz,np.kron(ox,id2))+(1/2)*np.identity(8)
    projector_x1_z2 = -(1/2)*np.kron(id2,np.kron(ox,oz))+(1/2)*np.identity(8)
    projector_all = projector_x1_z2@projector_z0_x1
    return projector_all

#finding eigenvalues and eigenvectors of a given matrix eg. Hamiltonian generated above. a 2d list, 0iter - 
def eigens(matrix):
    evals,evec = np.linalg.eig(matrix)
    return evals,np.transpose(evec)

#arrangement of eigenvectors with respect to eigenvalues from smallest to greatest
def arr(evals,evec):
    item = np.argsort(evals)
    vectors = []
    for i in range(len(evals)):
        vectors.append(evec[item[i]])
    return vectors

# Code to prepare the Hamiltonian
def Ham(mass,photon,interaction):
    id2 = np.identity(2)
    ox = np.array([[0,1],[1,0]])
    oy = np.array([[0+0j,0-1j],[0+1j,0+0j]])
    oz = np.array([[1,0],[0,-1]])
    h_gate = 1/np.sqrt(2)*(ox+oz)

    mass_term = mass*(1/2)*(np.kron(oz,np.kron(id2,id2))-np.kron(id2,np.kron(id2,oz)))
    photon_term = photon*(1/2)*np.kron(id2,np.kron(ox,id2))
    interaction_term = interaction*(1/4*np.kron(ox,np.kron(oz,ox))+(1/4)*np.kron(oy,np.kron(oz,oy)))
    return mass_term + photon_term + interaction_term
    
def state_circuit(mass, interaction, E_n=0):
    hamil = Ham(mass, 1, interaction)
    eigenvec = eigens(hamil)[1]
    eigenval = eigens(hamil)[0]
    arrangedvecs = arr(eigenval, eigenvec)
    
    projector_all = projall()
    
    projected = projector_all @ arrangedvecs
    gsvecs = gs(np.transpose(projected))

    vec_0 = gsvecs[0]
    vec_1 = gsvecs[1]


    ox = np.array([[0,1],[1,0]])
    oz = np.array([[1,0],[0,-1]])
    h_gate = 1/np.sqrt(2)*(ox+oz)
    #1-0 state
    minus10 = np.kron(np.kron(ox@[1,0],h_gate@ox@[1,0]),[1,0])

    # 0+1 state
    plus01 = np.kron(np.kron([1,0],h_gate@[1,0]),ox@[1,0])

    a_0 = np.dot(plus01, vec_0)
    a_1 = np.dot(plus01, vec_1)

    t_0 = 2 * np.arccos(a_0)
    t_1 = 2 * np.arccos(a_1)
    thetas = [t_0.real, t_1.real]

    q = QuantumRegister(3)
    qc = QuantumCircuit(q)
    qc.rx(thetas[E_n], 1)
    qc.h(1)
    qc.cnot(1, 0)
    qc.cnot(1, 2)
    qc.x(2)
    qc.h(1)

    return (qc)
    
def gs(A):
    A = np.transpose(A)
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
    
    
def energy_calculation(ham):   # Calculating the energy of the hamiltonian by using psi@Hamiltonian@psi equation
    arrangedevecs = arr(eigens(np.array(ham,float))[0],np.transpose(eigens(np.array(ham,float))[1]))
    projector_all = projall()
    projevecs = projector_all@arrangedevecs
    gevecs = gs(projevecs)@swap@cp
    Energ= []
    for i in range(len(gevecs)):
        vec = gevecs[i]
        Energ.append(np.dot((vec),np.matmul(ham,np.transpose(vec)))/np.dot(vec,vec))  
    return Energ
    
def magic_function(mass, interaction, E_n=0):
    hamil = Ham(mass, 1, interaction)
    eigenvec = eigens(hamil)[1]
    eigenval = eigens(hamil)[0]
    arrangedvecs = arr(eigenval, eigenvec)
    
    
    projector_all = projall()
    
    projected = projector_all @ arrangedvecs
    gsvecs = gs(projected) # gs(np.transpose(projected))

    vec_0 = gsvecs[0]
    vec_1 = gsvecs[1]

    a_0 = np.dot(plus01, vec_0)
    a_1 = np.dot(plus01, vec_1)

    t_0 = 2 * np.arccos(a_0).real
    t_1 = 2 * np.arccos(a_1).real

    thetas = [t_0,t_1]

    energies = energy_calculation(hamil)

    return thetas[E_n], energies[E_n].real
