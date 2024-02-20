import numpy as np
import matplotlib.pyplot as plt
from q_solve import generate_result
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from PIL import Image
from overlap import compute_schmidt_full


def compute_schmidt_2(result,idx,s=1):
    if s==1:
        a = compute_schmidt_states(result, idx, 0, 0)[0] #schmidt 1 on system 1
        a = np.squeeze(a)
        b = compute_schmidt_states(result, idx, 1, 0)[0] #schmidt 2 on system 1
        b=np.squeeze(b)
        g = np.outer(a, b).flatten()
        g=np.squeeze(g)
    elif s==2:
        a = compute_schmidt_states(result, idx, 0, 0)[1] #schmidt 1 on system 1
        a = np.squeeze(a)
        b = compute_schmidt_states(result, idx, 1, 0)[1] #schmidt 2 on system 1
        b=np.squeeze(b)
        g = np.outer(a, b).flatten()
        g=np.squeeze(g)
    return g


# Function to update the plot for each frame of the animation

def update_plot(frames):
    # Clear previous plot
    plt.clf()
    
    print(frames)
    state = compute_schmidt_2(result,frames,1)
    state2 = compute_schmidt_2(result,frames,2)
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    plt.plot(eigenenergies_total,energy_coeff);
    plt.plot(eigenenergies_total,energy_coeff2);
    plt.title(f"Plot of the probability that Schmidt1 and 2 are in the energy eigenstates for EI={EI})")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")

def update_plot1(frames):
    # Clear previous plot
    plt.clf()
    
    print(frames)
    state = compute_schmidt_2(result,frames,1)
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    plt.plot(eigenenergies_total,energy_coeff);
    plt.title(f"Plot of the probability that Schmidt1 and 2 are in the energy eigenstates for EI={EI})")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")

def update_plot2(frames):
    # Clear previous plot
    plt.clf()
    
    print(frames)
    state2 = compute_schmidt_2(result,frames,2)
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    plt.plot(eigenenergies_total,energy_coeff2);
    plt.title(f"Plot of the probability that Schmidt1 and 2 are in the energy eigenstates")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")

def make_gif_distribs1s2(d1,d2,w, E_spacing, Int_strength, tmax, ind_nb):#H_total,result,EI
    result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result(d1,d2,w, E_spacing, Int_strength, tmax, ind_nb,1)    
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 
    # Create a figure
    fig = plt.figure(figsize=(10, 2))
    
    # Create the animation
    ani = FuncAnimation(fig, update_plot, frames=100, interval=100)
    
    # Save the animation as a GIF
    ani.save('distrib_schmidt1_2_over_energy_spectrum.gif', writer='pillow')

def make_gif_distribs1(H_total,result,EI):
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 
    result = result
    EI = EI
    # Create a figure
    fig = plt.figure(figsize=(10, 2))
    
    # Create the animation
    ani = FuncAnimation(fig, update_plot1, frames=100, interval=100)
    
    # Save the animation as a GIF
    ani.save('distrib_schmidt1_over_energy_spectrum.gif', writer='pillow')

def make_gif_distribs2(H_total,result,EI):
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    
    result = result
    EI = EI
    # Create a figure
    fig = plt.figure(figsize=(10, 2))
    
    # Create the animation
    ani = FuncAnimation(fig, update_plot2, frames=100, interval=100)
    
    # Save the animation as a GIF
    ani.save('distrib_schmidt2_over_energy_spectrum.gif', writer='pillow')


