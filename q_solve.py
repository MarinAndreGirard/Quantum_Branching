import numpy
import qutip as qt
import math
import matplotlib.pyplot as plt

def test(x):
    a = numpy.linspace(0,x,10)
    return a 

def generate_result(d1 = 10,d2 = 200,w = 0.3, E_spacing = 1.0, Int_strength = 0.03, tmax= 10, ind_nb = 100,log=0):
# Define the Hilbert space dimensions for the two quantum systems
    dim_system_1 = d1  # Dimension of system 1. make both dimensions even
    dim_system_2 = d2  # Dimension of system 2 (changed to 10)
    dim_total = dim_system_1 * dim_system_2  # Total Hilbert space dimension
    
    # Create basis states for system 1 and system 2
    basis_system_1 = [qt.basis(dim_system_1, i) for i in range(dim_system_1)]
    basis_system_2 = [qt.basis(dim_system_2, i) for i in range(dim_system_2)]
    ket_0 = qt.basis(dim_system_1, round(d1/4))  # |0> state
    ket_1 = qt.basis(dim_system_1, round(d1*3/4))  # |2> state, int(dim_system_1/2)
    
    # Define random Hermitian matrices as Hamiltonians for system 1 and system 2
    H_system_1 = qt.qeye(dim_system_1) #qt.rand_herm(dim_system_1)  # Random Hermitian matrix for system 1
    energy_spacing = E_spacing  # Adjust as needed
    diagonal_elements = numpy.arange(0, dim_system_1) * energy_spacing
    H_q = qt.Qobj(numpy.diag(diagonal_elements)) # Create a diagonal matrix with increasing diagonal elements
    H_system_2 = qt.rand_herm(dim_system_2)  # Random Hermitian matrix for system 2

    # Define initial states for system 1 and system 2
    initial_state_system_1 = (math.sqrt(w)*ket_0 + math.sqrt(1-w)*ket_1).unit()
    #initial_state_system_2 = qt.rand_ket(dim_system_2)
    ev ,es = H_system_2.eigenstates()
    initial_state_system_2 = es[round(d2/2)]#qt.ket(es[round(d2/2)])  # You can change this to any basis state. note does not appear on schmidt spectral decomp, because is not an eigenstate of the rd matrix    
    #define initial state of full system
    state = qt.tensor(initial_state_system_1, initial_state_system_2)
    
    ### BASIC SIMU
    interaction_strength = Int_strength  # Adjust as needed
    #H_interaction = interaction_strength * (qt.rand_herm(dim_total)) #note here the interaction mat is fully random. 
    H_interaction = interaction_strength * qt.tensor(H_q, qt.rand_herm(dim_system_2))  
    #not a full rd matrix will make things non-zero in terms in which it was not nonzero in the sstem 1 basis
    
    # Expand the system Hamiltonians to the full Hilbert space dimensions
    H_system_1_ext = qt.tensor(H_system_1, qt.qeye(dim_system_2))
    H_system_2_ext = qt.tensor(qt.qeye(dim_system_1), H_system_2)
    
    # Define the total Hamiltonian
    H_total = H_system_1_ext + H_system_2_ext + H_interaction
    
    # Define the time settings for the simulation
    
    tlist = numpy.linspace(0, tmax, ind_nb)  # Adjust the time range and step count as needed

    import numpy as np

    if log == 0:
        tlist = numpy.linspace(0, tmax, ind_nb)  # Linear spacing
    elif log == 1:
        tlist = numpy.logspace(numpy.log10(1), numpy.log10(tmax+1), ind_nb)-1  # Logarithmic spacing
    else:
        raise ValueError("Invalid value for 'log'. It should be either 0 or 1.")
    
    # Perform time evolution of the combined system
    result = qt.mesolve(H_total, state, tlist, [], [])
    
    return result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2
