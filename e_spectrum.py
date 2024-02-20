import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt


def Neff(H_total,result):
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 
    state = result.states[0]
    p2=[(abs(np.vdot(state, eigenstate)) ** 2) ** 2 for eigenstate in eigenstates_total]
    Neff = 1/np.sum(p2)
    return Neff
    
def plot_e_spectrum(H_total, result,EI):
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 
    state = result.states[0]
    energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
    c = np.count_nonzero(energy_coeff)
    N = Neff(H_total,result)
    print(f"Neff_total is {N}")
    num_bins=100
    min_energy=min(eigenenergies_total)
    max_energy=max(eigenenergies_total)
    plt.figure(figsize=(10, 2))
    plt.plot(eigenenergies_total,energy_coeff);
    plt.title(f"Plot of the probability that the global state be in an energy eigenstate for EI={EI})")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Probabilities")
    plt.show()
    #looks like a nice wigner semicircle, this is the thing, who's shape changes as the interaction energy increases. lets check that
    plt.figure(figsize=(10, 2))
    plt.hist(eigenenergies_total, bins=num_bins, range=(min_energy, max_energy), edgecolor='black');
    plt.title(f"Distribution of the spectrum of H_total for EI={EI})")
    plt.xlabel("Eigenenergies of H_total")
    plt.ylabel("Count")
    plt.show()
    
    #PERFECT wigner semi cirlce, vs weird cowboy hat