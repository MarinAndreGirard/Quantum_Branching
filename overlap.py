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
            mean = get_mean_rd_overlap(w[i],EI[j])
            axs[i, j].axhline(y=mean, color='red', linestyle='--')
            
    plt.legend(["o01", "o02", "012"])
    plt.tight_layout()
    plt.savefig(f'Graphs/overlap_characterization_EI_{EI},w_{w},env_{env}.png')
    plt.show()


#This is the piece of code I used to find the mean overlap between newly initialized eigenstates. But note that this is for a very specific set of parameters
def get_mean_rd_overlap(w = 0.3,Int_strength = 0.052):
    d1, d2 = 10, 200
    E_spacing = 1.0
    
    # Create basis states for system 1 and system 2
    basis_system_1 = [qt.basis(d1, i) for i in range(d1)]
    basis_system_2 = [qt.basis(d2, i) for i in range(d2)]
    ket_0 = qt.basis(d1, 3)  # |0> state
    ket_1 = qt.basis(d1, 7)  # |2> state, int(dim_system_1/2)
        
    # Define random Hermitian matrices as Hamiltonians for system 1 and system 2
    H_system_1 = qt.qeye(d1) #qt.rand_herm(dim_system_1)  # Random Hermitian matrix for system 1
    energy_spacing = E_spacing  # Adjust as needed
    diagonal_elements = np.arange(0, d1) * 1.0
    H_q = qt.Qobj(np.diag(diagonal_elements)) # Create a diagonal matrix with increasing diagonal elements
    H_system_2_1 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_2 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_3 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_4 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_5 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_6 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_7 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_8 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_9 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    H_system_2_10 = qt.rand_herm(d2,1)  # Random Hermitian matrix for system 2
    # Define initial states for system 1 and system 2
    initial_state_system_1 = (math.sqrt(w)*ket_0 + math.sqrt(1-w)*ket_1).unit()
    #initial_state_system_2 = qt.rand_ket(dim_system_2)
    ev1 ,es1 = H_system_2_1.eigenstates()
    ev2 ,es2 = H_system_2_2.eigenstates()
    ev3 ,es3 = H_system_2_3.eigenstates()
    ev4 ,es4 = H_system_2_4.eigenstates()
    ev5 ,es5 = H_system_2_5.eigenstates()
    ev6 ,es6 = H_system_2_6.eigenstates()
    ev7 ,es7 = H_system_2_7.eigenstates()
    ev8 ,es8 = H_system_2_8.eigenstates()
    ev9 ,es9 = H_system_2_9.eigenstates()
    ev10,es10 = H_system_2_10.eigenstates()

    initial_state_system_2_1 = es1[round(d2/2)]
    initial_state_system_2_2 = es2[round(d2/2)]
    initial_state_system_2_3 = es3[round(d2/2)]
    initial_state_system_2_4 = es4[round(d2/2)]
    initial_state_system_2_5 = es5[round(d2/2)]
    initial_state_system_2_6 = es6[round(d2/2)]
    initial_state_system_2_7 = es7[round(d2/2)]
    initial_state_system_2_8 = es8[round(d2/2)]
    initial_state_system_2_9 = es9[round(d2/2)]
    initial_state_system_2_10 = es10[round(d2/2)]
    #define initial state of full system
    states=[]
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_1))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_2))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_3))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_4))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_5))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_6))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_7))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_8))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_9))
    states.append(qt.tensor(initial_state_system_1, initial_state_system_2_10))

    interaction_strength = Int_strength  # Adjust as needed
    H_interaction = interaction_strength * qt.tensor(H_q, qt.rand_herm(d2,1))  
        
    H_system_1_ext = qt.tensor(H_system_1, qt.qeye(d2))
    H_system_2_ext = 0.75*qt.tensor(qt.qeye(d1), H_system_2_1)
    H_total = H_system_1_ext + H_system_2_ext + H_interaction

    eigenenergies_total, eigenstates_total = H_total.eigenstates() 

    st=[]
    for s in states:
        st.append(s.full().squeeze())

    state_0=st[0]
    p_0=[abs(np.vdot(state_0, eigenstate)) for eigenstate in eigenstates_total]
    st.pop(0)

    overlap_list=[]
    for s in st:
        p = [abs(np.vdot(s, eigenstate)) for eigenstate in eigenstates_total]
        overlap_list.append(np.dot(p_0, p))

    mean_overlap = np.mean(overlap_list)
    return mean_overlap
