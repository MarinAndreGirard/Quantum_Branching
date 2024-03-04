import numpy as np
import matplotlib.pyplot as plt
from q_solve import generate_result
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from PIL import Image
#from overlap import compute_schmidt_full
from Schmidt_solve import compute_schmidt_full

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
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 

    # Create a figure
    fig = plt.figure(figsize=(10, 2))
    
    # Create the animation
    ani = FuncAnimation(fig, update_plot, frames=100, interval=100)
    
    # Save the animation as a GIF
    ani.save('distrib_schmidt1_2_over_energy_spectrum.gif', writer='pillow')

def make_gif_distribs1(H_total,result,EI):
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




def update_plot_new(frames,result,eigenstates_total,eigenenergies_total,EI,w,env):
    # Clear previous plot
    plt.clf()
    frames = frames+1
    state = compute_schmidt_full(result,frames,1)
    state2 = compute_schmidt_full(result,frames,2)
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    plt.plot(eigenenergies_total, energy_coeff)
    plt.plot(eigenenergies_total, energy_coeff2)
    plt.title(f"Plot of the probability that Schmidt1 and 2 are in the energy eigenstates for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.legend(["Schmidt1", "Schmidt2"])
    plt.ylim(0, 0.35)
    # Calculate the mean
    mean1 = np.sum(np.array(energy_coeff) * np.array(eigenenergies_total))
    mean2 = np.sum(np.array(energy_coeff2) * np.array(eigenenergies_total))
    st1_tst1 = np.mean((np.array(energy_coeff) * np.array(eigenenergies_total)-mean1)**2)
    st1_tst2 = np.mean((np.array(energy_coeff2) * np.array(eigenenergies_total)-mean2)**2)
    std1 = np.std(np.array(energy_coeff) * np.array(eigenenergies_total))
    std2 = np.std(np.array(energy_coeff2) * np.array(eigenenergies_total))
    # Add a vertical line at the mean for energy_coeff
    plt.axvline(x=mean1, color='b', linestyle='--')
    # Add a vertical line at the mean for energy_coeff2
    plt.axvline(x=mean2, color='r', linestyle='--')
    # Add a vertical line at the mean plus one standard deviation for energy_coeff
    plt.axvline(x=mean1 + st1_tst1, color='g', linestyle='--')
    # Add a vertical line at the mean minus one standard deviation for energy_coeff
    plt.axvline(x=mean1 - st1_tst1, color='g', linestyle='--')
    # Add a vertical line at the mean plus one standard deviation for energy_coeff2
    plt.axvline(x=mean2 + st1_tst2, color='c', linestyle='--')
    # Add a vertical line at the mean minus one standard deviation for energy_coeff2
    plt.axvline(x=mean2 - st1_tst2, color='c', linestyle='--')

    #plt.legend()



def make_gif_distribs1s2_new(EI,w,result,eigenstates_total,eigenenergies_total,env,d1,d2,E_spacing,tmax,ind_nb):
    print(env)
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_new,fargs=(result,eigenstates_total,eigenenergies_total,EI,w,env), frames=99, interval=100)
    # Save the animation as a GIF
    if env!=[0]:
        e=env[0:2]
        ani.save(f'Gifs/distrib_schmidt1_2_over_energy_spectrum_EI_{EI}_w_{w}_env_{e}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    else:
        ani.save(f'Gifs/distrib_schmidt1_2_over_energy_spectrum_EI_{EI}_w_{w}_env_NA_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')        
    plt.close()
'''
def plot_strd_time(EI,w,result,eigenstates_total,eigenenergies_total,env,d1,d2,E_spacing,tmax,ind_nb):
    for i in range(len(tlist)-1):

'''