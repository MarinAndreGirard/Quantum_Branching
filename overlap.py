import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_states
from Schmidt_solve import compute_schmidt_full

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
