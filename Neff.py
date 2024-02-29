import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_states
from Schmidt_solve import compute_schmidt_full


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
    coef3 = []
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




