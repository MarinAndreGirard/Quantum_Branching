import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from q_solve import generate_result
from VN import plot_VN
from VN import compute_VN_time
from Schmidt_solve import compute_schmidt_states

######This is a checkup to see if you program works

def checkup():

    #defining all quantities used in the simulation
    d1, d2 = 10, 200
    w = 0.3# smallest weight of the 2 eingestates, set btw 0 and 1
    E_spacing = 1.0
    Int_strength = 0.03
    # Define the time settings for the simulation
    tmax= 10
    ind_nb = 100
    #using all above to obtain QM qtts
    result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result(d1,d2,w, E_spacing, Int_strength, tmax, ind_nb,1)
    
    # Eigenstates and eigenenergies 
    eigenenergies_system_2, eigenstates_system_2 = H_system_2.eigenstates() 
    eigenenergies_system_total, eigenstates_system_total = H_total.eigenstates() 
    eigenenergies_system_1, eigenstates_system_1 = H_q.eigenstates() 
    # Create a 2x2 grid of plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    v2 = compute_VN_time(result,tlist)
    
    # Plot 1: VN entropy over time
    axes[0, 0].plot(tlist, v2)
    axes[0, 0].set_title("VN entropy over time")
    axes[0, 0].set_xlabel("time index")
    axes[0, 0].set_ylabel("VN entropy")
    
    # Plot 2: Distribution of the Environment state over the environment energy eigenstates
    time_index = 0  # Define the time index
    s0 = compute_schmidt_states(result, time_index, 1)[0]
    schmidt_coefficients0 = [abs(np.vdot(s0, eigenstate)) ** 2 for eigenstate in eigenstates_system_2]
    axes[0, 1].plot(eigenenergies_system_2, schmidt_coefficients0, marker='o', label=f'Energy {eigenenergies_system_2}')
    axes[0, 1].set_title("Distribution of the Environment state over the environment energy eigenstates")
    axes[0, 1].set_xlabel("Energy Eigenstates")
    axes[0, 1].set_ylabel("Schmidt Coefficients")
    
    
    
    # Plot 3: Distribution of the System state over the system energy eigenstates
    s0 = compute_schmidt_states(result, time_index, 0)[0]
    schmidt_coefficients0 = [abs(np.vdot(s0, eigenstate)) ** 2 for eigenstate in eigenstates_system_1]
    axes[1, 0].plot(eigenenergies_system_1, schmidt_coefficients0, marker='o', label=f'Energy {eigenenergies_system_1}')
    axes[1, 0].set_title("Distribution of the System state over the system energy eigenstates")
    axes[1, 0].set_xlabel("Energy Eigenstates")
    axes[1, 0].set_ylabel("Schmidt Coefficients")
    
    # Plot 4: Distribution of the state over the energy eigenstates
    s0 = result.states[0]
    schmidt_coefficients0 = [abs(np.vdot(s0, eigenstate)) ** 2 for eigenstate in eigenstates_system_total]
    axes[1, 1].plot(eigenenergies_system_total, schmidt_coefficients0, marker='o', label=f'Energy {eigenenergies_system_total}')
    axes[1, 1].set_title("Distribution of the state over the energy eigenstates")
    axes[1, 1].set_xlabel("Energy Eigenstates")
    axes[1, 1].set_ylabel("Schmidt Coefficients")
    
    plt.tight_layout()  # Adjust spacing between subplots
    plt.show()