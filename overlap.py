import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_states
from Schmidt_solve import compute_schmidt_full
from q_solve import generate_result

#Defining the overlap of states in probability space
# a measure of similarity of weights in the energie eigenbasis

def get_p_s2(state,eigenstates_total):
    p=[abs(np.vdot(state, eigenstate)) ** 2  for eigenstate in eigenstates_total]
    return p

#def get_p(state,eigenstates_total):
#    p=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
#    return p

#def get_p_2(state,eigenstates_total):
#    p=[(abs(np.vdot(state, eigenstate)) ** 2) **2 for eigenstate in eigenstates_total]
#    return p

def p_overlap(state1,state2,eigenstates_total):
    sqrt_p1 = np.sqrt(get_p_s2(state1,eigenstates_total))
    sqrt_p2 = np.sqrt(get_p_s2(state2,eigenstates_total))
    overlap = np.dot(sqrt_p1, sqrt_p2)
    return overlap

'''
def compute_schmidt_full(result, idx,s=1):
    if s==1:
        a = compute_schmidt_states(result, idx, 0, 0)[0] #schmidt 1 on system 1
        b = compute_schmidt_states(result, idx, 1, 0)[0] #schmidt 2 on system 1
        g = np.tensordot(a, b, 0)
    elif s ==2:
        c = compute_schmidt_states(result, idx, 0, 0)[1] #schmidt 1 on system 1
        d = compute_schmidt_states(result, idx, 1, 0)[1] #schmidt 2 on system 1
        g = np.tensordot(c, d, 0)
    elif s==3:
        e = compute_schmidt_states(result, idx, 0, 0)[2] #schmidt 1 on system 1
        f = compute_schmidt_states(result, idx, 1, 0)[2] #schmidt 2 on system 1
        g = np.tensordot(e, f, 0)
    else: 
        print("wrong input value")
    return g
'''

def plot_p_overlap_graph(tlist,result,H_total,w,EI,log_scale=False):
    o01 = []
    o02 = []
    o12 = []
    eigenenergies_total, eigenstates_total = H_total.eigenstates()
    for idx in range(len(tlist)-2):
        s1=compute_schmidt_full(result,idx+1,1)
        s2=compute_schmidt_full(result,idx+1,2)
        global_state = result.states[idx+1]
        #s3=compute_schmidt_full(result,idx,3)
        o01.append(p_overlap(global_state,s1,eigenstates_total))
        o02.append(p_overlap(global_state,s2,eigenstates_total))
        o12.append(p_overlap(s1,s2,eigenstates_total))

    plt.figure(figsize=(7, 4))
    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
    else:
        plt.xscale("linear")
        plt.yscale("linear")
    plt.plot(tlist[0:len(tlist)-2],o01)
    plt.plot(tlist[0:len(tlist)-2],o02)
    plt.plot(tlist[0:len(tlist)-2],o12)
    plt.title("Graph of the overlap of global v schmidt 1 v schmidt 2")
    plt.xlabel("Time")
    plt.ylabel(r"$\sqrt{\vec{pi}} \cdot \sqrt{\vec{pj}}$")
    plt.legend(['overlap Global-Schmidt 1', 'overlap Global-Schmidt 2', 'overlap Schmidt 1-Schmidt 2'])
    plt.savefig(f'Graphs/overlap_EI_{EI}_w_{w}.png')


def plot_p_overlap_graph_characterize(d1=10,d2=200,w=[0.1,0.2,0.3,0.4], E_spacing=1.0, EI=[0.03,0.05,0.07,0.09],tmax=30, ind_nb =100,env=[0]): #
    
    o01_list = []
    o02_list = []
    o12_list = []
    for wi in w:  
        for EIi in EI:
            print(f'wi_{wi}_EI_{EIi}')           
            result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result(d1,d2,wi, E_spacing, EIi, tmax, ind_nb,1)
            eigenenergies_total, eigenstates_total = H_total.eigenstates()
            print("generated")
            o01 = []
            o02 = []
            o12 = []
            for idx in range(len(tlist)-1):
                s1=compute_schmidt_full(result,idx+1,1)
                s2=compute_schmidt_full(result,idx+1,2)
                global_state = result.states[idx+1]
                #s3=compute_schmidt_full(result,idx,3)
                o01.append(p_overlap(global_state,s1,eigenstates_total))
                o02.append(p_overlap(global_state,s2,eigenstates_total))
                o12.append(p_overlap(s1,s2,eigenstates_total))
            o01_list.append(o01)
            o02_list.append(o02)
            o12_list.append(o12)
        
    fig, axs = plt.subplots(len(w), len(EI), figsize=(10, 2*len(w)), sharex=True, sharey=True)
    #plt.title(f"Plot of the means and standard dev of the distributions of Schmidt 1 and 2 w={w}, EI={EI}, env={env}", fontsize=10)
    for i in range(len(w)):
        for j in range(len(EI)):
            axs[i, j].plot(tlist[0:len(tlist) - 1], o01_list[i*len(EI) + j])
            axs[i, j].plot(tlist[0:len(tlist) - 1], o02_list[i*len(EI) + j])
            axs[i, j].plot(tlist[0:len(tlist) - 1], o12_list[i*len(EI) + j])
            #axs[i, j].set_title(f"Plot of the means and standard dev of the distributions of Schmidt 1 and 2 w={w[i]}, EI={EI[j]}, env={env}", fontsize=5)
            axs[i, j].set_xlabel("Time", fontsize=6)
            axs[i, j].set_ylabel("Overlap", fontsize=8)
    
    plt.legend(["o01", "o02", "012"])
    plt.tight_layout()
    plt.savefig(f'Graphs/overlap_characterization_EI_{EI},w_{w},env_{env}.png')
    plt.show()
