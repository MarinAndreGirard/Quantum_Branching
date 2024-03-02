import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from Schmidt_solve import compute_schmidt_full
from Schmidt_solve import compute_schmidt_states_new
from scipy.spatial import distance

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

def time_cos_similarity_plot(result,tlist):
    #need to that the schmidt states become stable over time
    #to do that, we get a list of ss[0] vectors over time and check that as timer goes, they change less and less.
    se0_list=[]
    for i in range(len(tlist)-1):
        ss, se, sv = compute_schmidt_states_new(result, i+1)
        se = np.abs(se[0]).flatten()
        se0_list.append(se)
    # Calculate the similarity between ss[0] vectors at different time indices
        step = 1
    similarities = []
    for i in range(len(tlist)-step-1):
        d = 1-distance.cosine(se0_list[i], se0_list[i+step])
        similarities.append(d)

    # Plot the similarity over time
    plt.figure(figsize=(10, 6))
    plt.plot(tlist[0:len(tlist)-step-1], similarities)
    plt.xlabel('Time')
    plt.ylabel('Similarity')
    plt.title(f'Evolution of Similarity between environment schmidt 0 with itself a little before. Step size {step}')


    plt.show()