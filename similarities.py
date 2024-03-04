import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from q_solve import generate_result
from Schmidt_solve import compute_schmidt_full
from Schmidt_solve import compute_schmidt_states_new
from scipy.spatial import distance
from VN import plot_VN

def probs_schmidt_in_energy_eigenstates(result,eigenenergies_total, eigenstates_total,tlist,EI,w):
    prob_list1 = []
    prob_list2 = []
    for idx in range(len(tlist)-1):
        state = compute_schmidt_full(result,idx+1,1)
        state2 = compute_schmidt_full(result,idx+1,2)
        energy_coeff=[abs(np.vdot(state, eigenstate)) ** 2 for eigenstate in eigenstates_total]
        energy_coeff2=[abs(np.vdot(state2, eigenstate)) ** 2 for eigenstate in eigenstates_total]
        prob_list1.append(energy_coeff)
        prob_list2.append(energy_coeff2)
    return prob_list1, prob_list2


def cos_similarity_btw_s1_s2_plot(s1_list,s2_list,tlist):
    similarities = []
    for idx in range(len(tlist)-1):
        d = 1-distance.cosine(s1_list[idx], s2_list[idx])
        similarities.append(d)
    plt.figure(figsize=(10, 6))
    plt.plot(tlist[0:len(tlist)-1], similarities)
    plt.xlabel('Time')
    plt.ylabel('cosine similarity')
    plt.title('Evolution of similarity between schmidt 1 and schmidt 2')


def metric_similarity_btw_s1_s2_plot(s1_list,s2_list,tlist):
    similarities = []
    for idx in range(len(tlist)-1):
        d=0
        for i in range(len(s1_list[idx])):
            d = d + s1_list[idx][i]+s2_list[idx][i]-abs(s1_list[idx][i]-s2_list[idx][i])
        similarities.append(d/2)
    plt.figure(figsize=(10, 6))
    plt.plot(tlist[0:len(tlist)-1], similarities)
    plt.xlabel('Time')
    plt.ylabel('similarity metric')
    plt.title('Evolution of similarity between schmidt 1 and schmidt 2')

def time_cos_similarity_plot_new(result,tlist,H_total,EI,w,step=1):
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 
    p1, p2 = probs_schmidt_in_energy_eigenstates(result,eigenenergies_total, eigenstates_total,tlist,EI,w)
    p1 = np.sqrt(p1)
    p2 = np.sqrt(p2)
    similarities = []
    similarities2 = []
    for i in range(len(tlist)-step-1):
        d = 1-distance.cosine(p1[i], p1[i+step])
        similarities.append(d)
        d1= 1-distance.cosine(p2[i], p2[i+step])
        similarities2.append(d1)

    # Plot the similarity over time
    plt.figure(figsize=(10, 6))
    plt.plot(tlist[0:len(tlist)-step-1], similarities)
    plt.plot(tlist[0:len(tlist)-step-1], similarities2)
    plt.xlabel('Time')
    plt.ylabel('Similarity')
    plt.title(f'Evolution of Similarity between environment schmidt 0 with itself a little before. Step size {step}')
    plt.legend(["Env Schmidt 0", "Env Schmidt 1"])


    plt.show()

def time_cos_similarity_plot(result,tlist):
    #need to that the schmidt states become stable over time
    #to do that, we get a list of ss[0] vectors over time and check that as timer goes, they change less and less.
    se0_list=[]
    se1_list=[]    
    for i in range(len(tlist)-1):
        ss, se, sv = compute_schmidt_states_new(result, i+1)
        se1 = np.abs(se[0]).flatten()
        #print(f"i{i}se{se}")
        se0_list.append(se1)
        se2 = np.abs(se[1]).flatten()
        se1_list.append(se2)
    # Calculate the similarity between ss[0] vectors at different time indices
        step = 1
    similarities = []
    similarities2 = []
    for i in range(len(tlist)-step-1):
        d = 1-distance.cosine(se0_list[i], se0_list[i+step])
        similarities.append(d)
        d1= 1-distance.cosine(se1_list[i], se1_list[i+step])
        similarities2.append(d1)

    # Plot the similarity over time
    plt.figure(figsize=(10, 6))
    plt.plot(tlist[0:len(tlist)-step-1], similarities)
    plt.plot(tlist[0:len(tlist)-step-1], similarities2)
    plt.xlabel('Time')
    plt.ylabel('Similarity')
    plt.title(f'Evolution of Similarity between environment schmidt 0 with itself a little before. Step size {step}')
    plt.legend(["Env Schmidt 0", "Env Schmidt 1"])


    plt.show()


def similarity_btw_s1_s2_plot_compare_w(d1,d2,w, E_spacing, EI,tmax=10,ind_nb=100):
    #defining all quantities used in the simulation
    d1, d2 = 10, 200
    E_spacing = 1.0
    Int_strength = 0.03

    w_vec=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    for w in w_vec:
        result, tlist, H_q, H_system_2, H_system_1_ext, H_system_2_ext, H_interaction, H_total, ket_0, ket_1, initial_state_system_2 = generate_result(d1,d2,w, E_spacing, Int_strength, tmax, ind_nb,0)
        eigenenergies_total, eigenstates_total = H_total.eigenstates()
        plot_VN(result,tlist)
        s1_list,s2_list=probs_schmidt_in_energy_eigenstates(result,eigenenergies_total,eigenstates_total,tlist,EI,w)
        similarities = []
        for idx in range(len(tlist)-1):
            d = 1-distance.cosine(s1_list[idx], s2_list[idx])
            similarities.append(d)
            
        plt.plot(tlist[0:len(tlist)-1], similarities)
        plt.xlabel('Time')
        plt.ylabel('cosine similarity')
        plt.title('Evolution of similarity between schmidt 1 and schmidt 2')
        plt.legend(['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'])
    plt.show()
