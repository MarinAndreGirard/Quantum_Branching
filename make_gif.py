import numpy as np
import matplotlib.pyplot as plt
from q_solve import generate_result
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from PIL import Image
#from overlap import compute_schmidt_full
from Schmidt_solve import compute_schmidt_full
from Schmidt_solve import compute_schmidt_states_new

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
    plt.legend(["Schmidt1", "Schmidt2", "Mean1", "Mean2", "Mean1 + Std1", "Mean1 - Std1", "Mean2 + Std2", "Mean2 - Std2"])
    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)


def update_plot_new_zoomed(frames,result,eigenstates_total,eigenenergies_total,EI,w,env):
    # Clear previous plot
    plt.clf()
    frames = frames+1
    state = compute_schmidt_full(result,frames,1)
    state2 = compute_schmidt_full(result,frames,2)
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    plt.plot(eigenenergies_total, energy_coeff)
    plt.plot(eigenenergies_total, energy_coeff2)
    plt.title(f"Zoomed plot of the probability that Schmidt1 and 2 are in the energy eigenstates for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.legend(["Schmidt1", "Schmidt2"])
    plt.ylim(0, 0.02)
    plt.xlim(0 , 2)

    # Calculate the me
    mean1 = np.sum(np.array(energy_coeff) * np.array(eigenenergies_total))
    mean2 = np.sum(np.array(energy_coeff2) * np.array(eigenenergies_total))
    st1_tst1 = np.mean((np.array(energy_coeff) * np.array(eigenenergies_total)-mean1)**2)
    st1_tst2 = np.mean((np.array(energy_coeff2) * np.array(eigenenergies_total)-mean2)**2)
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

    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

def update_plot_pointer(frames,result,eigenstates_int,eigenenergies_int,EI,w,env):
    # Clear previous plot
    plt.clf()
    frames = frames+1
    
    state = compute_schmidt_full(result,frames,1) 
    state2 = state = compute_schmidt_full(result,frames,2)
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_int]
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_int]
    plt.plot(eigenenergies_int, energy_coeff)
    plt.plot(eigenenergies_int, energy_coeff2)
    plt.title(f"Plot of the probability that the schmidts  are in the energy eigenstates of H_int for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.legend(["Schmidt1", "Schmidt2"])
    plt.ylim(0, 0.25)
    #plt.xlim(5, 8)

    # Calculate the mean
    mean1 = np.sum(np.array(energy_coeff) * np.array(eigenenergies_int))
    mean2 = np.sum(np.array(energy_coeff2) * np.array(eigenenergies_int))
    st1_tst1 = np.mean((np.array(energy_coeff) * np.array(eigenenergies_int)-mean1)**2)
    st1_tst2 = np.mean((np.array(energy_coeff2) * np.array(eigenenergies_int)-mean2)**2)
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

    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

def update_plot_pointers1(frames,result,eigenstates_int,eigenenergies_int,EI,w,env):
    # Clear previous plot
    plt.clf()
    frames = frames+1
    state = compute_schmidt_full(result,frames,1) 
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_int]
    plt.plot(eigenenergies_int, energy_coeff)
    plt.title(f"Plot of the probability that the schmidt state 1 is in the energy eigenstates of H_int for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.legend(["Schmidt 1"])
    plt.ylim(0, 0.25)

    #plt.xlim(5, 8)

    # Calculate the mean
    mean1 = np.sum(np.array(energy_coeff) * np.array(eigenenergies_int))
    st1_tst1 = np.mean((np.array(energy_coeff) * np.array(eigenenergies_int)-mean1)**2)
    # Add a vertical line at the mean for energy_coeff
    plt.axvline(x=mean1, color='b', linestyle='--')
    # Add a vertical line at the mean plus one standard deviation for energy_coeff
    plt.axvline(x=mean1 + st1_tst1, color='g', linestyle='--')
    # Add a vertical line at the mean minus one standard deviation for energy_coeff
    plt.axvline(x=mean1 - st1_tst1, color='g', linestyle='--')
    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

def update_plot_pointers2(frames,result,eigenstates_int,eigenenergies_int,EI,w,env):
    # Clear previous plot
    plt.clf()
    frames = frames+1
    
    state2 = state = compute_schmidt_full(result,frames,2)
    energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_int]
    plt.plot(eigenenergies_int, energy_coeff2)
    plt.title(f"Plot of the probability that the schmidt state 2  are in the energy eigenstates of H_int for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.legend(["Schmidt 2"])
    plt.ylim(0, 0.25)
    #plt.xlim(5, 8)

    # Calculate the mean
    mean2 = np.sum(np.array(energy_coeff2) * np.array(eigenenergies_int))
    st1_tst2 = np.mean((np.array(energy_coeff2) * np.array(eigenenergies_int)-mean2)**2)
    # Add a vertical line at the mean for energy_coeff2
    plt.axvline(x=mean2, color='r', linestyle='--')
    # Add a vertical line at the mean plus one standard deviation for energy_coeff2
    plt.axvline(x=mean2 + st1_tst2, color='c', linestyle='--')
    # Add a vertical line at the mean minus one standard deviation for energy_coeff2
    plt.axvline(x=mean2 - st1_tst2, color='c', linestyle='--')

    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)


def update_plot_interference(frames,result,eigenstates_int,eigenenergies_int,EI,w,env):
    # Clear previous plot
    frames = frames+1
    ss, se, sv = compute_schmidt_states_new(result, frames)
    s_val_0 = sv[0]
    s_val_1 = sv[1]

    state = result.states[frames].full()
    state1 = compute_schmidt_full(result,frames,1)
    state2 = compute_schmidt_full(result,frames,2)
    interference = [abs(np.sqrt(s_val_0*s_val_1)*(2*np.real(np.vdot(eig,state1))*np.real(np.vdot(eig,state2))+2*np.imag(np.vdot(eig,state1))*np.imag(np.vdot(eig,state2)))) for eig in eigenstates_int]
    plt.clf()
    plt.plot(eigenenergies_int, interference)
    
    plt.title(f"Plot of the interference btw s1 and s2 in the energy eigenbasis (H_I) for EI={EI} and w={w} and env={env}")
    plt.xlabel("Eigenenergies of H")
    plt.ylabel("Interference")
    plt.ylim(0, 0.10)
    #plt.xlim(5, 8)

    # Calculate the mean
    mean1 = np.sum(np.array(interference) * np.array(eigenenergies_int))
    st1_tst1 = np.mean((np.array(interference) * np.array(eigenenergies_int)-mean1)**2)
    # Add a vertical line at the mean for energy_coeff
    plt.axvline(x=mean1, color='b', linestyle='--')
    # Add a vertical line at the mean plus one standard deviation for energy_coeff
    plt.axvline(x=mean1 + st1_tst1, color='g', linestyle='--')
    # Add a vertical line at the mean minus one standard deviation for energy_coeff
    plt.axvline(x=mean1 - st1_tst1, color='g', linestyle='--')
    # Add clock
    plt.text(0.95, 0.95, f"Frame: {frames}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)


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


def make_gif_distribs1s2_new_zoomed(EI,w,result,eigenstates_total,eigenenergies_total,env,d1,d2,E_spacing,tmax,ind_nb):
    
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_new_zoomed,fargs=(result,eigenstates_total,eigenenergies_total,EI,w,env), frames=99, interval=100)
    # Save the animation as a GIF
    if env!=[0]:
        e=env[0:2]
        ani.save(f'Gifs/zoomed_distrib_schmidt1_2_over_energy_spectrum_EI_{EI}_w_{w}_env_{e}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    else:
        ani.save(f'Gifs/zoomed_distrib_schmidt1_2_over_energy_spectrum_EI_{EI}_w_{w}_env_NA_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')        
    plt.close()


def make_gif_distrib_pointer(EI,w,result,eigenstates_int,eigenenergies_int,env,d1,d2,E_spacing,tmax,ind_nb):
    
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_pointer,fargs=(result,eigenstates_int,eigenenergies_int,EI,w,env), frames=99, interval=100)
    # Save the animation as a GIF
    if env!=[0]:
        e=env[0:2]
        ani.save(f'Gifs/zoomed_distrib_pointer_over_energy_spectrum_EI_{EI}_w_{w}_env_{e}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    else:
        ani.save(f'Gifs/zoomed_distrib_pointer_over_energy_spectrum_EI_{EI}_w_{w}_env_NA_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')        
    plt.close()


def make_gif_distrib_pointer_s1(EI,w,result,eigenstates_int,eigenenergies_int,q1,d1,d2,E_spacing,tmax,ind_nb):
    
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_pointers1,fargs=(result,eigenstates_int,eigenenergies_int,EI,w,[0]), frames=99, interval=100)
    # Save the animation as a GIF
    ani.save(f'Gifs/Pointer_gifs/distrib_pointers1_over_energy_spectrum_EI_{EI}_w_{w}_q1_{q1}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    plt.close()
    

def make_gif_distrib_pointer_s2(EI,w,result,eigenstates_int,eigenenergies_int,q1,d1,d2,E_spacing,tmax,ind_nb):
    
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_pointers2,fargs=(result,eigenstates_int,eigenenergies_int,EI,w,[0]), frames=99, interval=100)
    # Save the animation as a GIF
    ani.save(f'Gifs/Pointer_gifs/distrib_pointers2_over_energy_spectrum_EI_{EI}_w_{w}_q1_{q1}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    plt.close()


def make_gif_distrib_interf(EI,w,result,eigenstates_int,eigenenergies_int,q1,d1,d2,E_spacing,tmax,ind_nb):
    
    # Create a figure
    fig = plt.figure(figsize=(10, 5))

    # Create the animation
    ani = FuncAnimation(fig, update_plot_interference,fargs=(result,eigenstates_int,eigenenergies_int,EI,w,[0]), frames=99, interval=100)
    # Save the animation as a GIF
    ani.save(f'Gifs/Pointer_gifs/distrib_interference_over_energy_spectrum_EI_{EI}_w_{w}_q1_{q1}_d1_{d1}_d2_{d2}_Espace_{E_spacing}_tmax_{tmax}_ind_nb_{ind_nb}.gif', writer='pillow')
    plt.close()
