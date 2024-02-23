import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt

def test(x):
    return x+1

def compute_schmidt_states_new(result, time_index):
    global_state = result.states[time_index]
    density_matrix = qt.ptrace(global_state, [0]) # Calculate the density matrix at the specified time
    eigenvalues, eigenstates = density_matrix.eigenstates() # Compute the eigenstates and eigenvalues of the density matrix
    eigenstates = [np.array(state) for state in eigenstates]
    schmidt_states_s = []
    schmidt_values = []
    i=0
    for state, eigenvalue in zip(eigenstates, eigenvalues):
        schmidt_values.append(eigenvalue)
        if eigenvalue < 10e-14:
            # If the eigenvalue is zero, set the Schmidt state to a zero vector
            schmidt_states_s.append(np.zeros_like(state))
        else:
            #print(f"state {state}")
            i=i+1
            N=abs(np.vdot(state.conjugate(),state))
            schmidt_states_s.append(state/np.sqrt(N)) # Normalize

    # Sort the Schmidt states by eigenvalue in descending order
    schmidt_states_s, schmidt_values = zip(*sorted(zip(schmidt_states_s, schmidt_values), key=lambda x: -x[1]))
    d=np.size(global_state)
    d1 = np.size(schmidt_states_s[0])
    d2=d//d1
    #compute the schmidt states of the environement.
    schmidt_states_e = []
    I = np.eye(d2)
    for j in range(i):
        state = schmidt_states_s[j]
        P_a_state = np.kron(np.outer(state,state.conjugate().T),I)     #def projector |><|xId, outer transposes teh second one
        temp = np.dot(P_a_state,global_state) #vdot is conjugate on first one.
        #print(np.sum(temp))
        N = abs(np.vdot(global_state,temp))
        #print(N)
        schmidt_states_e.append(temp/np.sqrt(N))
    #print(schmidt_values)
    return schmidt_states_s,schmidt_states_e,schmidt_values


def compute_schmidt_full(result,idx,s=1):
    ss, se, sv = compute_schmidt_states_new(result, idx)
    if s==1:
        a = ss[0] #schmidt 1 on system 1
        a = np.squeeze(a)
        b = se[0] #schmidt 1 on system 2
        b=np.squeeze(b)
        g = np.outer(a, b).flatten()
        g=np.squeeze(g)
    elif s==2:
        a = ss[1] #schmidt 2 on system 1
        a = np.squeeze(a)
        b = se[1] #schmidt 2 on system 2
        b=np.squeeze(b)
        g = np.outer(a, b).flatten()
        g=np.squeeze(g)
    return g



def compute_schmidt_states(result, time_index, subsystem_index=0, trigger=0):
    density_matrix = qt.ptrace(result.states[time_index], [subsystem_index]) # Calculate the density matrix at the specified time
    eigenvalues, eigenstates = density_matrix.eigenstates() # Compute the eigenstates and eigenvalues of the density matrix
    eigenstates = [np.array(state) for state in eigenstates]
    schmidt_states = []
    schmidt_values = []
    for state, eigenvalue in zip(eigenstates, eigenvalues):
        schmidt_values.append(eigenvalue)
        if eigenvalue == 0:
            # If the eigenvalue is zero, set the Schmidt state to a zero vector
            schmidt_states.append(np.zeros_like(state))
        else:
            schmidt_states.append(state / np.linalg.norm(state)) # Normalize
            #schmidt_states.append(state / np.sqrt(np.array(eigenvalue)))  # Normalize
    # Sort the Schmidt states by eigenvalue in descending order
    schmidt_states, schmidt_values = zip(*sorted(zip(schmidt_states, schmidt_values), key=lambda x: -x[1]))
    if trigger == 0:
        return schmidt_states
    else:
        return schmidt_values

