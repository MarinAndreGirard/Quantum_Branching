import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_states



#computign the neffschmidt/neffttot
#delta = Neff_schmidt/Neff_tot

def plot_Neff_schmidt(H_total, result):
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    coef1 = []
    coef2 = []
    coef3 = []

    for idx in range(len(tlist)):
        a = compute_schmidt_states(result, idx, 0, 0)[0] #schmidt 1 on system 1
        b = compute_schmidt_states(result, idx, 1, 0)[0] #schmidt 2 on system 1
        c = compute_schmidt_states(result, idx, 0, 0)[1] #schmidt 1 on system 1
        d = compute_schmidt_states(result, idx, 1, 0)[1] #schmidt 2 on system 1
        e = compute_schmidt_states(result, idx, 0, 0)[2] #schmidt 1 on system 1
        f = compute_schmidt_states(result, idx, 1, 0)[2] #schmidt 2 on system 1
        g = np.tensordot(a, b, 0)
        h = np.tensordot(c, d, 0)
        i = np.tensordot(e, f, 0)
        p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2__=[(abs(np.vdot(i, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        coef1.append(1/np.sum(p2))
        coef2.append(1/np.sum(p2_))
        coef3.append(1/np.sum(p2__))
    
    plt.figure(figsize=(10, 2))
    plt.xscale("log")
    plt.plot(coef1)
    plt.plot(coef2)
    plt.plot(coef3)

def plot_Neff_schmidt_Neff_tot(H_total, result,tlist,EI):
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    p=[(abs(np.vdot(result.states[1], eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
    coef1 = []
    coef2 = []
    coef3 = []
    coeft = 1/np.sum(p)

    for idx in range(len(tlist)):
        a = compute_schmidt_states(result, idx, 0, 0)[0] #schmidt 1 on system 1
        b = compute_schmidt_states(result, idx, 1, 0)[0] #schmidt 2 on system 1
        c = compute_schmidt_states(result, idx, 0, 0)[1] #schmidt 1 on system 1
        d = compute_schmidt_states(result, idx, 1, 0)[1] #schmidt 2 on system 1
        e = compute_schmidt_states(result, idx, 0, 0)[2] #schmidt 1 on system 1
        f = compute_schmidt_states(result, idx, 1, 0)[2] #schmidt 2 on system 1
        g = np.tensordot(a, b, 0)
        h = np.tensordot(c, d, 0)
        i = np.tensordot(e, f, 0)
        p2=[(abs(np.vdot(g, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2_=[(abs(np.vdot(h, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        p2__=[(abs(np.vdot(i, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
        coef1.append(1/np.sum(p2))
        coef2.append(1/np.sum(p2_))
        coef3.append(1/np.sum(p2__))
    
    c1=coef1/coeft
    c2=coef2/coeft
    c3=coef3/coeft
    plt.figure(figsize=(10, 2))
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(c1)
    plt.plot(c2)
    plt.plot(c3)
    plt.title(f"Graph of Neff_schmidt/Neff_tot over time for EI={EI}")
    plt.xlabel("Time")
    plt.ylabel("Neff_schmidt/Neff_tot")
    plt.legend(['Schmidt 1', 'Schmidt 2','Schmidt 3'])
