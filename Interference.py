import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from Schmidt_solve import compute_schmidt_states_new
from Schmidt_solve import compute_schmidt_full

def interference_plot(result,H_total,tlist,start_num=1000,end_num=1010):
    #set things up
    prob_list=[]
    prob_list2=[]
    prob_list3=[]
    eigenenergies_total, eigenstates_total = H_total.eigenstates() 

    #get the schmidt values over time
    s_val_0=[]
    s_val_1=[]
    t_ind = range(len(tlist))
    for idx in t_ind:
        ss, se, sv = compute_schmidt_states_new(result, idx)
        s_val_0.append(sv[0])
        s_val_1.append(sv[1])


    # Create an array of integers between start_num and end_num
    arr = np.arange(start_num, end_num + 1)
    arr2 = np.arange(start_num, end_num)
    for i in arr:
        eig = eigenstates_total[i]
        prob=[]
        prob2=[]
        prob3=[]
        for idx in range(len(tlist)-1):
            prob.append(abs(np.vdot(compute_schmidt_full(result,idx+1,1), eig)) ** 2)
            prob2.append(abs(np.vdot(compute_schmidt_full(result,idx+1,2), eig)) ** 2)
            prob3.append(abs(np.vdot(result.states[idx].full(), eig)) ** 2)
            
        prob_list.append(prob)
        prob_list2.append(prob2)
        prob_list3.append(prob3)
        print(i)

    '''
    plt.figure(figsize=(10, 6))
    for prob in prob_list:
        plt.plot(tlist[0:99], prob)
    plt.title('Probability for s1 to be in 10 different total energy eigenstate over time')
    plt.xlabel('Time')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.legend([f'Prob in eigenstate {i}' for i in arr])
    plt.show()

    plt.figure(figsize=(10, 6))
    for prob in prob_list2:
        plt.plot(tlist[0:99], prob)
    plt.title('Probability for s2 to be in 10 different total energy eigenstate over time')
    plt.xlabel('Time')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.legend([f'Prob in eigenstate {i}' for i in arr])
    plt.show()

    plt.figure(figsize=(10, 6))
    for prob in prob_list3:
        plt.plot(tlist[0:99], prob)
    plt.title('Probability that the global state be in 10 different energy eigenstate over time')
    plt.xlabel('Time')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.legend([f'Prob in eigenstate {i}' for i in arr])
    plt.show()
    '''
    k=0
    b = 0
    for i in arr2:
        if prob_list[i-start_num+1][0]>b:
            b=prob_list[i-start_num+1][0]
            k=i-start_num+1
        #if prob_list3[i-start_num+1][2]>prob_list3[i-start_num][2]:
        #    k=i-start_num+1

    i = k
    j=start_num+k

    eig = eigenstates_total[j]
    interf=[]
    for idx in range(len(tlist)-1):
        interf.append(np.sqrt(s_val_0[idx]*s_val_1[idx])*(2*np.real(np.vdot(eig,compute_schmidt_full(result,idx+1,1)))*np.real(np.vdot(eig,compute_schmidt_full(result,idx+1,2)))+2*np.imag(np.vdot(eig,compute_schmidt_full(result,idx+1,1)))*np.imag(np.vdot(eig,compute_schmidt_full(result,idx+1,2)))))
        #interf.append(np.sqrt(s_val_0[idx]*s_val_1[idx])*(np.vdot(eig,compute_schmidt_full(result,idx+1,1))*(np.vdot(eig,compute_schmidt_full(result,idx+1,2)).conjugate()) + np.vdot(eig,compute_schmidt_full(result,idx+1,2))*(np.vdot(eig,compute_schmidt_full(result,idx+1,1)).conjugate())))

    weighted = np.multiply(prob_list[i], s_val_0[0:len(tlist)-1])+np.multiply(prob_list2[i], s_val_1[0:len(tlist)-1])
    weighted_plus_interf = weighted + interf
    plt.figure(figsize=(10, 6))
    plt.xscale('linear')
    plt.plot(tlist[0:len(tlist)-1], weighted)
    plt.plot(tlist[0:len(tlist)-1], weighted_plus_interf)
    plt.plot(tlist[0:len(tlist)-1], prob_list3[i])
    plt.plot(tlist[0:len(tlist)-1],np.multiply(prob_list[i], s_val_0[0:len(tlist)-1]))
    plt.plot(tlist[0:len(tlist)-1],np.multiply(prob_list2[i], s_val_1[0:len(tlist)-1]))
    plt.legend(['Weighted sum of schmidts without interference','Weighted sum of schmidts with interference','global state','s1','s2'])
    plt.title(f'Probability that various states be in energy eigenstates {j} over time')
    plt.xlabel('Time')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.show()


#What this graph tells me is that somehow the sum of the probabilities of the 2 worlds being in  
#an energy eigenstate is not conserved over time.

# I need to do a more through examination of this for more states
#I need to verify theoretically that this is what we expect.
#I need to make sure thjat the multiplication by schmidt weight mak es sense.

