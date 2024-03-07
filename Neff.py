import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_states
from Schmidt_solve import compute_schmidt_full
from q_solve import generate_result
from q_solve import generate_result_envi_superpo


#computign the neffschmidt/neffttot
#delta = Neff_schmidt/Neff_tot

def plot_Neff_schmidt(H_total,tlist, result,log=True):
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    coef1 = []
    coef2 = []
    coef3 = []


    for idx in range(len(tlist)-1):
        g=compute_schmidt_full(result,idx,s=1)
        h=compute_schmidt_full(result,idx,s=2)
        
        p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        coef1.append(1/np.sum(p2))
        coef2.append(1/np.sum(p2_))
        
    plt.figure(figsize=(10, 2))
    if log:
        plt.xscale("log")
    else:
        plt.xscale("linear")
    plt.yscale("log")
    plt.plot(tlist[0:len(tlist)-1],coef1)
    plt.plot(tlist[0:len(tlist)-1],coef2)
    
def plot_Neff_schmidt_Neff_tot(H_total, result,tlist,EI,log=True):
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    p=[(abs(np.vdot(result.states[1], eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
    coef1 = []
    coef2 = []
    coeft = 1/np.sum(p)

    for idx in range(len(tlist)-1):
        g=compute_schmidt_full(result,idx+1,s=1)
        h=compute_schmidt_full(result,idx+1,s=2)
        
        p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        coef1.append(1/np.sum(p2))
        coef2.append(1/np.sum(p2_))
        
    c1=coef1/coeft
    c2=coef2/coeft
    plt.figure(figsize=(10, 2))
    if log:
        plt.xscale("log")
    else:
        plt.xscale("linear")
    plt.yscale("log")
    plt.plot(tlist[0:len(tlist)-1], c1)
    plt.plot(tlist[0:len(tlist)-1], c2)
    plt.title(f"Graph of Neff_schmidt/Neff_tot over time for EI={EI}")
    plt.xlabel("Time")
    plt.ylabel("Neff_schmidt/Neff_tot")
    plt.legend(['Schmidt 1', 'Schmidt 2'])
    plt.savefig(f'Graphs/Neff_schmidt_Neff_tot_{EI}.png')



def plot_Neff_schmidt_Neff_tot_characterize(d1=10,d2=200,w=[0.1,0.2,0.3,0.4], E_spacing=1.0, EI=[0.03,0.05,0.07,0.09],tmax=10, ind_nb =100,env=[0]): #
    delta_list_1 = []
    delta_list_2 = []
    if env==[0]:
        for wi in w:  
            for EIi in EI:
                print(f'wi_{wi}_EI_{EIi}')           
                result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result(d1,d2,wi, E_spacing, EIi, tmax, ind_nb,1)
                eigenenergies_total, eigenstates_total = H_total.eigenstates()
                print("generated")
                p=[(abs(np.vdot(result.states[1], eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                coef1 = []
                coef2 = []
                coeft = 1/np.sum(p)

                for idx in range(len(tlist)-1):
                    g=compute_schmidt_full(result,idx+1,s=1)
                    h=compute_schmidt_full(result,idx+1,s=2)
                    
                    p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                    p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                    coef1.append(1/np.sum(p2))
                    coef2.append(1/np.sum(p2_))
                print("Schmidt calculated")
                c1=coef1/coeft
                c2=coef2/coeft
                delta_list_1.append(c1)
                delta_list_2.append(c2)

    else:
        for EIi in EI:
            for wi in w: 
                print(f'wi{wi}_EI_{EIi}') 
                result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result_envi_superpo(d1,d2,wi, E_spacing, EIi, tmax, ind_nb,1,env)
                eigenenergies_total, eigenstates_total = H_total.eigenstates() 
                print("generated")
                p=[(abs(np.vdot(result.states[1], eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                coef1 = []
                coef2 = []
                coeft = 1/np.sum(p)
                
                for idx in range(len(tlist)-1):
                    g=compute_schmidt_full(result,idx+1,s=1)
                    h=compute_schmidt_full(result,idx+1,s=2)
                    
                    p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                    p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
                    coef1.append(1/np.sum(p2))
                    coef2.append(1/np.sum(p2_))
                print("Schmidt calculated")
                c1=coef1/coeft
                c2=coef2/coeft
                delta_list_1.append(c1)
                delta_list_2.append(c2)


    fig, axs = plt.subplots(len(w), len(EI), figsize=(10, 2*len(w)), sharex=True, sharey=True)

    for i in range(len(w)):
        for j in range(len(EI)):
            axs[i, j].plot(tlist[0:len(tlist) - 1], delta_list_1[i*len(EI) + j])
            axs[i, j].plot(tlist[0:len(tlist) - 1], delta_list_2[i*len(EI) + j])
            axs[i, j].set_title(f"Graph of Neff_schmidt/Neff_tot over time for w={w[i]}, EI={EI[j]}, env={env}", fontsize=5)
            axs[i, j].set_xlabel("Time", fontsize=6)
            axs[i, j].set_ylabel("Neff_schmidt/Neff_tot", fontsize=8)
            axs[i, j].legend(['Schmidt 1', 'Schmidt 2'])
            axs[i, j].set_yscale("log")  # Set y-axis to logarithmic scale
            axs[i, j].set_xscale("log")  # Set y-axis to logarithmic scale


    plt.tight_layout()
    plt.savefig(f'Graphs/Neff_schmidt_Neff_characterization_EI_{EI},w_{w},env_{env}.png')
    plt.show()

