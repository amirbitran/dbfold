#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:24:08 2020

@author: amirbitran
"""


import hmmlearn
from hmmlearn.hmm import MultinomialHMM   #this one is used, do not delete!

import numpy as np

import copy as cp
import matplotlib.pyplot as plt
import dbfold.utils as utils
import dbfold.load_data as load_data
import dbfold.compute_PMF as compute_PMF
import natsort

import scipy

import joblib

def all_folding_rates(Arrhenius_temps, transitions_of_interest, data, clusters, PDB_files,  protein, unique_tops, top_free_energies, eq_temperatures,  Ntrials,  min_trans):
    """    
    First, computes Arrhenius parameters, as well as bootstrap distributions for them 
    Next, computes folding rates at all temperatures without errors, just using the raw Arrhenius data
    But everything is saved, so it's easy to get a distribution of folding rates at any desired temperature
    
    Note: transitions_of_interest should be a list of UNFOLDING transitions, as in Arrhenius fit
    """
    activation_energies, prefactors, Bootstrapped_EAs, Bootstrapped_intercepts= bootstrap(Arrhenius_temps, transitions_of_interest, data, PDB_files, protein, legend_loc='upper right', Ntrials=Ntrials, min_trans= min_trans)    
    G=compute_PMF.cluster_free_energies( clusters, unique_tops, top_free_energies, eq_temperatures)
    folding_rates, unfolding_rates, unused=infer_folding_rates(clusters,activation_energies, prefactors, G, eq_temperatures)
    #pickle.dump([folding_rates, unfolding_rates, activation_energies, prefactors, Bootstrapped_EAs, Bootstrapped_intercepts,G, deltaG_trap, temperatures, clusters, substructures], open( '{}/{}'.format(PMF_directory, filename), 'wb' ))
    folding_info = {'folding rates': folding_rates,
                    'unfolding rates': unfolding_rates,
                    'activation energies': activation_energies,
                    'prefactors': prefactors,
                    'Bootstrapped EAs': Bootstrapped_EAs,
                    'Bootstrapped intercepts': Bootstrapped_intercepts,
                    'cluster free energies': G,
                    'eq temperatures': eq_temperatures}
    print('Folding rate calculation complete!')
    return folding_info

def Arrhenius_fit(temperatures, transitions_of_interest, data, PDB_files, protein, min_trans=5,  labels = [],colors=['r', 'b', 'g', 'y', 'c', 'm', 'r', 'b', 'g', 'y', 'c', 'm'], ax = None, legend_loc='upper right',  legend = True, fontsize = 35,  relabel_trans = None, labelsize = 30, legend_fontsize = 35,show_rsquared = True, temp_norm = 1):
    """
    Sweeps through temperatures and corresponding trajectories in temperatures and trajectories_list
    For each temperature, calculates transition rates between clusters
    Then for each transition rate, applies a linear fit to calculate Arrhenius prefactor and activation energy
    
    
    you can set transitions_of_interest to 'All', in which case you extract all transitions that occur during unfolding simultions
    
    You also have the option of entering signlets into transitions_of_interest, for instnace, [(3,), (,4)]
    What this then does is, for each singlet in the list, you compute the rate of getting to that cluster, 
    not worrying too much about what happened before (i.e. implicitly setting everything before to 0)
    Make sure you have the comma in each parentheses: For instnace, (3,) instead of (3)
    
    
    Also, by default, a transition of interest is only counted at a given temperatures if it is observed at least 5 times
    YOu can change this minimum number using keyword min_trans
    
    Finally, you have the option to relabel_trans in case the cluster numbering is unintuitive
    By default, this variable is set to None, so no transitions are relabeled
    But you can also feed it a list of tuples, in which case, the ith transition is renamed based on that ith tuple
    If that transition is originally called 1->3 for instance. but the ith tuple is (1,2), then it will be relabled as 1->2
    
    Oh also, if temp_norm is not 1, we express x axis as temp_norm/T, where temp_norm is presuambly the melting temp
    """
    n_clusters=int(np.nanmax(data))+1


    #mean_FPTs=np.zeros(( n_clusters, n_clusters, len(temperatures)))
    mean_transition_rates=np.zeros(( n_clusters, n_clusters, len(temperatures)))  #same as in Transition_rates function, with one page for each temperature
    Ns=np.zeros(( n_clusters, len(temperatures)))  #how many trajectories go into calculating each rate?    
    

  
    #Loop through all temps and obtain transition rates
    
    dics = []  #Added 1/18/19, a list of dictionaries. Each dic tells you, for a given temperature, how many counts of each transition occur
    
    for t, temp in enumerate(temperatures):
        print("Computing temperature {}".format(temp))
        dics.append(get_unfolding_transitions(data, PDB_files, temp=temp, verbose=False)) #added 1/18/19

        transmat, times, Ns[:,t]=markov_transmat(temp, data, PDB_files)
        rates=rates_from_transmat(transmat, times[1]-times[0])
        mean_transition_rates[:,:,t]=rates 

        #else:  
        #    mean_FPTs[:,:,t], mean_transition_rates[:,:,t], Ns[:,:,t] = Transition_rates(temp, transitions_of_interest, data, PDB_files, gmix, Filter)
    activation_energies=np.zeros((n_clusters, n_clusters ))
    prefactors=np.zeros((n_clusters, n_clusters ))
    
    
    
    if ax == None:
        fig, ax = plt.subplots() 
    
    if temp_norm ==1:
        ax.set_xlabel('1/Temperature', fontsize=fontsize)
    else:
        ax.set_xlabel('$T_M$/Temperature', fontsize=fontsize)
    ax.set_ylabel('Log of inverse MFPT', fontsize=fontsize)
    ax.set_title('{}'.format(protein), fontsize=fontsize, y=1.03)
    ax.tick_params(axis='both', which='major', labelsize=labelsize, pad=8) 

    #now, we loop through all transitions and compute their Arrhenius plots 
    zzz=0   
    for i in range(n_clusters):
        for j in range(n_clusters):      
            current_rates=mean_transition_rates[i,j,:]
            
            #do we care to include this transition in fit??
            if transitions_of_interest=='All':
                if j>i: fit=True
                else: fit=False
            else:
                if (i,j) in transitions_of_interest: 
                    fit=True
                    n = transitions_of_interest.index((i,j))
                else:
                    fit=False            
            for k in range(len(current_rates)):
                if current_rates[k]==0:
                    current_rates[k]=np.nan
            
            if fit:  #added 1/18/19
                eliminate = np.array([z for z in range(len(current_rates)) if (i,j) in dics[z].keys() and dics[z][(i,j)]<min_trans]) #added 1/18/19
                #print(eliminate)
                if len(eliminate)>0: current_rates[eliminate] =np.nan
            
            mean_transition_rates[i,j,:] = current_rates #update the mean_transition_rates matrix to reflect any new Nans that have been added due to insufficient number of transitions
            
            non_Nan=np.where(np.logical_and(np.logical_not(np.isnan(current_rates)) , np.logical_not(np.isinf(current_rates))))[0] #entires that are neither NaN nor infinity
            current_rates=current_rates[non_Nan]            
            current_temperatures=[temperatures[z] for z in non_Nan] 

            #print(non_Nan)
            if len(current_rates)>1 and fit:  #only do an Ahrrenious fit if you have data over multiple temperatures
                 #coeffs = np.polyfit([1/T for T in current_temperatures], np.log(current_rates), 1, full = False)
                 slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([1/T for T in current_temperatures], np.log(current_rates))
                 #print((slope, intercepts))
                 coeffs = [slope, intercept]
                 activation_energies[i,j]=-coeffs[0]
                 prefactors[i,j]=np.exp(coeffs[1])
                 print(r_value**2)
                 
                 
                 

                 if relabel_trans!=None:
                     iii = relabel_trans[n][0]
                     jjj = relabel_trans[n][1]
                 else:
                     iii = i
                     jjj = j
                
                
                 if len(labels)==0:
                     label='{} -> {}'.format(iii,jjj)
                 else:
                     label = labels[n]
                 ax.scatter([temp_norm/T for T in current_temperatures], np.log(current_rates), label = label, s=50, color=colors[n])
                 #ax.plot([1/T for T in current_temperatures], coeffs[1]+[coeffs[0]*1/T for T in current_temperatures], label='{} -> {} fit'.format(iii,jjj), color=colors[n])
                 ax.plot([temp_norm/T for T in current_temperatures], coeffs[1]+[coeffs[0]*1/T for T in current_temperatures], color=colors[n])
                 ax.annotate('$E_A = {}$ $k_B T$'.format( round(-coeffs[0], 1) ), (1/current_temperatures[1]+0.02,  0.4+coeffs[1]+coeffs[0]*1/current_temperatures[1] ), fontsize=fontsize )
                 if show_rsquared: ax.annotate('$R^2$ = {}'.format( round( r_value**2, 3)), (1/current_temperatures[1]+0.02,  -0.05+coeffs[1]+coeffs[0]*1/current_temperatures[1] ), fontsize=fontsize )
                 zzz+=1
                 
    if legend: ax.legend(fontsize=legend_fontsize, loc=legend_loc) 

            
    return activation_energies, prefactors, mean_transition_rates, Ns


def Bernie_elimination(Assignments, PDB_files, Sanders_thresh, plot=True, cumulative=False, keep_zero_times=True):
    """
    Eliminate assignments that do not represent more than some fraction of the population given by Sanders_thresh
    e.g. Sanders_thresh may be The 1 Percent (i.e.Sanders_thresh=0.01)
    
    Can also do this in a cumulative fashion, so that we sort unique assignemnts by representation, and keep however many we need 
    to cumulatively represent fraction of the population 1 - Sanders_thresh
    
    if keep_zero_times=True, then we will never eliminate assignments that occur at the beginning of a simulation
    """     
    if Sanders_thresh==0 and cumulative:
        Sanders_thresh=-0.1  #to avoid rounding errors when summing the cumulative fraction

    unique=list(set(Assignments))
    times=utils.get_times(PDB_files)
    zero_points=list(set([a for t,a in enumerate(Assignments) if times[t]==0]))  #points that occur at beginning of simulations
    rep=[]
    for u in unique:
        rep.append(len([i for i in Assignments if i==u])/len(Assignments))
    
    
    #We whipe out clusters that represent too small a fraction of the population. We do this in one of two fashions
    
    Bernie_dic={} #this dic will tell you what fraction of the population is represented by a given cluster, OR 1 minus the cumulative fraction accounted for by clusters up to that one when the clusters are sorted in order
    if cumulative:
        sorted_indices=np.argsort(-np.array(rep))
        sorted_points=[unique[s] for s in sorted_indices]
        sorted_rep=[rep[s] for s in sorted_indices]
        cum_rep=[np.sum(sorted_rep[0:i+1]) for i in range(len(unique))]
        
        for i, u in enumerate(sorted_points):
            Bernie_dic[u] = 1 - cum_rep[i]
    else:
        for i,u in enumerate(unique):
            Bernie_dic[u]=rep[i]
    
    #Determine which points we accept or reject
    if keep_zero_times: #Accept points with sufficient representation, AND points that occur at beginnings of simulations
        Accepted=[u for u in unique if Bernie_dic[u]>=Sanders_thresh or u in zero_points]
    else: #Only accept points if they have sufficient representation  
         Accepted=[u for u in unique if Bernie_dic[u]>=Sanders_thresh]
    
    BankstoBust=[u for u in unique if u not in Accepted]  #represent less than Sanders_thresh percent of the population

    #We assign a numerical label to all structures that are accepted
    sums=[-sum([int(i) for i in str]) for str in Accepted] #(negative of) number of formed substructures in each config
    sort=natsort.natsorted(zip(sums, Accepted))
    dic={}  #a dictionary to readily map assignments to labels 
    key=[]
    for i, pair in enumerate(sort):
        struct=pair[1]
        dic[struct]=i
        key.append(struct)    
    #####    
    for u in BankstoBust: dic[u]='BustBigBank'   #for now, rejected structures will not be given a label
    
    labels=[dic[a] for a in Assignments] 
    rejected_indices = np.array([i for i in range(len(labels)) if  labels[i]=='BustBigBank'])
    accepted_indices = np.array([i for i in range(len(labels)) if  labels[i]!='BustBigBank'])

    
    #Now we fix rejected structures
    
    for r in rejected_indices:

        prev_accept=accepted_indices[accepted_indices<r][-1]
        future_accepts=accepted_indices[accepted_indices>r]
        
        if len(future_accepts)==0: #no structures got accepted after this one so we go with the previous one
            prev_accept_struct = Assignments[prev_accept]
            labels[r]=dic[prev_accept_struct]
        else:
            next_accept=future_accepts[0]
            prev_accept_struct=Assignments[prev_accept]
            next_accept_struct=Assignments[next_accept]
            curr_struct=Assignments[r]
            
            #How different is this structure fromt the previous and the next accepted structure
            d_to_prev = np.sum([ np.abs( int(curr_struct[i])- int(prev_accept_struct[i])) for i in range(len(curr_struct))])
            d_to_next = np.sum([ np.abs( int(curr_struct[i])- int(next_accept_struct[i])) for i in range(len(curr_struct))])
            
            if d_to_prev<=d_to_next:
                labels[r]=dic[prev_accept_struct]
            else:
                labels[r]=dic[next_accept_struct] 
    
    if plot:
        utils.histogram(labels, key)

    return labels, key



def bootstrap(temperatures, transitions_of_interest, data, PDB_files, protein, legend_loc='upper right', Ntrials=2000, min_trans = 1):
    """
    A bootstrap for errors in Arrenhius fits
    At each temperature, we sample with replacement a set of N trajectories, where N is the original number of trajectories at that temperature
    We then use all these resampled trajectories at all temperatures to fit an Arrenhius plot, and store the resulting intercept and negative slope (activation energy)
    """
    n_clusters=np.max(data)+1
    All_times=[] 
    All_labels=[] #each page corresponds to the set of trajectories corresponding to a given temperature
    
    #First, we do an Arrenhius fit on the original data
    activation_energies, prefactors, orig_transition_rates, Ns=Arrhenius_fit(temperatures, transitions_of_interest, data, PDB_files, protein,  legend_loc='lower left', min_trans = min_trans)
    
    #print(np.shape(orig_transition_rates))
    #Now, we construct arrays that contain the labels and time values for all the trajectories
    
    for temp in temperatures: 
        times, labels=predict_unfolding_at_temperature(temp, data, PDB_files)
        All_times.append(times)
        All_labels.append(labels)
    
    #4/20 commented out the following:
    #All_labels=np.transpose(All_labels, (1,2,0)) #reshapes so that page corresponds to temperture
    All_times=np.array(All_times)
    
    Bootstrapped_EAs=np.full((n_clusters, n_clusters, Ntrials), np.nan) #will store all activation energies over all iterations of bootstrap. Each page is one iteration. Within a page, entry (i,j) correspodns to activation energy for i->j transition in that iteration
    Bootstrapped_intercepts=np.full((n_clusters, n_clusters,Ntrials),np.nan) #Same, but for the pre-factors
    
    
    mean_transition_rates=np.zeros(( n_clusters, n_clusters, len(temperatures), Ntrials)) 
    print('Running {} bootstrap trials to compute error on unfolding rates...'.format(Ntrials))
    for n in range(Ntrials): #loop throuhg bootstrap iterations
        if np.mod(n,100)==0: print('{} trials completed'.format(n))
        #A matrix that will store the transition rates for the current bootstrap iteration

        for t,temp in enumerate(temperatures):
            labels=All_labels[t][np.random.randint(np.shape(All_labels[t])[0], size=np.shape(All_labels[t])[0]), :]  #randomly choose, with replacement, some set of trajectories at the current temperature..the number of samples you choose is the number you originally had at that temp
            
            #obtain all transition rates for current samples
            transmat, unused_outupt, Ns=markov_transmat(temp, data, PDB_files, sim_labels=labels)
            rates=rates_from_transmat(transmat, times[1]-times[0])
            mean_transition_rates[:,:,t, n]=rates 

        #now compute Arrenhius parameters for current sample trajectories
        for i in range(n_clusters):
            for j in range(n_clusters):      
                current_rates=mean_transition_rates[i,j,:,n]
                for k in range(len(current_rates)): 
                    if current_rates[k]==0:
                        current_rates[k]=np.nan
                    elif np.isinf(current_rates[k]):
                        current_rates[k]=np.nan
                    elif np.isnan(orig_transition_rates[i,j,k]): #added 1/19: Does not do bootstrapping on transitions that were very rare in original Arrhenius plot 
                        current_rates[k]=np.nan
                non_Nan=np.where(np.logical_not(np.isnan(current_rates)) )[0]
                current_rates=current_rates[non_Nan]
                current_temperatures=[temperatures[i] for i in non_Nan ] 
                 
                if len(current_rates)>1:  #only do an Ahrrenious fit if you have data over multiple temperatures
                     coeffs=np.polyfit([1/T for T in current_temperatures], np.log(current_rates), 1)
                     Bootstrapped_EAs[i,j, n]=-coeffs[0]
                     Bootstrapped_intercepts[i,j, n]=coeffs[1]
    print('Bootstrap complete')
    return  activation_energies, prefactors, Bootstrapped_EAs, Bootstrapped_intercepts


def cluster(T_A, unique_labels, distance, S, key):
    """
    The parameter is:
        T_A: The kinetic distance threshold...two topological configurations will be clustered together if their characteristic exchange time is less than this
    """
    ##### Do loop clustering ####
    print('Doing loop clustering...')
    clusters, unused, mean_intercluster, mean_intracluster=utils.loopCluster(T_A, unique_labels, distance, verbose=False)
    
    X=[]
    for s in S:
        clust=[c for c in range(len(clusters)) if s in clusters[c]][0]
        X.append(clust)
    X=np.array(X)
    
    Clusters = []
    for b in clusters: 
        curr_clus = []
        for a in b:
            state = key[a]
            curr_clus.append(state)
        curr_clus = utils.barcodes_to_labels(curr_clus)
        Clusters.append(curr_clus)
    return X, Clusters


def compute_cluster_survival(i,times, sim_labels, min_length = 50, fit = 'Single'):
    """"
    Simply computes probability of surviving in cluster i
    
    sim_labels should be your data (cluster occupancy over time) in an array where each row is an independent trajectory and each column is a time
    
    min_length tells you we're only gonna keep sub_trajectories whose length is at least that length 
    
    fit can be 'Single' or 'Double' for single or double exponential fit, respectively
    
    """
    sub_trajectories=[]   #List of all trajectories, truncated at early times until the first instance of cluster i
    for k, trajec in enumerate(sim_labels):
        traj=cp.deepcopy(trajec)
        if i in traj:
            traj=traj[np.min(np.where( traj==i ) [0]):-1 ]  #Truncate the beginning of the  trajectory before cluster i is reached for the first time
         
            if len(traj)>=min_length:
                sub_trajectories.append(traj)
    
    #Find shortest common length among all trajectories we've kept
    N=len(sub_trajectories)  #how many trajectories go into rate calculation for this temperature?
    if N>0:
        largest_common_length=min(len(T) for T in sub_trajectories) #largest length that is shared by all subtrajectories. Ex. if the subtrajectories are [[1,1,2,2] , [1,1,2]], largest_common_length is 3 
        sub_trajectories=[T[0:largest_common_length] for T in sub_trajectories]

        sub_trajectories=np.array(sub_trajectories) #we can convert to np array since all subtrajectories now have the same length

        survival_probability=[ len(np.where(sub_trajectories[:,t]==i)[0])/np.shape(sub_trajectories)[0] for t in range(np.shape(sub_trajectories)[1])  ]
        
        times=times[0:len(survival_probability)]     
        
        
        if survival_probability[1]!=0 and survival_probability[-1]!=1:       #if survival probability immediately drops from 1 to 0, we cannnot solve for lambda MLE.  Similarly, we can't do this if transition never happens (so survival probability is always 1). Otherwise we can 
            if fit == 'Single':
                params,b=scipy.optimize.curve_fit(utils.exp_decay, times, survival_probability, p0=0.0000001 )

            else:
                params, b = scipy.optimize.curve_fit(utils.double_exp_decay, times, survival_probability, p0=[0.0000001, 0.0000001, 1, 1], bounds = ([0,0, -np.inf, -np.inf], [np.inf,np.inf,np.inf,np.inf]) )
            
        else:
            params=[np.nan]
    
    else: 
        params=[np.nan]
        survival_probability=[]
    
    #if lamb_MLE<0: lamb_MLE=np.nan
        
    return params, survival_probability, times, N




def get_folding_distribution(Bootstrapped_EAs, Bootstrapped_intercepts, free_energies, temperatures, temp, plot=None):
    """
    For a given temperature temp, computes a list of log folding rates 
    For now, we assume that there is no uncertainty regarding the free energy differences, although this may not be a great assumption
    Later, we can add in the uncertainty and treat it as Gaussian (also an approximation), so that we can do a convolution
    
    plot is a tuple that tells you which transition you care to plot
    """
    t=np.where(np.array(temperatures)==temp)[0][0]
    G=free_energies[t,:]
    Log_unfolding_distribution=get_unfolding_distribution(Bootstrapped_EAs, Bootstrapped_intercepts, temp)
    Ntrials=np.shape(Bootstrapped_EAs)[2]
    deltaG=np.zeros(np.shape(Log_unfolding_distribution[:,:,0]))
    for i in range(np.shape(deltaG)[0]):
        for j in range(np.shape(deltaG)[1]):
            deltaG[i,j]=G[j]-G[i]
    deltaG=np.repeat(deltaG[:, :, np.newaxis], Ntrials, axis=2) 
    Log_folding_distribution=Log_unfolding_distribution+deltaG
    Log_folding_distribution=np.transpose(Log_folding_distribution, (1,0,2))
    
    if plot!=None:
        plt.figure()
        plt.hist(Log_folding_distribution[plot[0], plot[1], :], bins=int(Ntrials/10), normed=1)
        plt.xlabel('Log folding rate', fontsize=25)
        plt.ylabel('Probability', fontsize=25)
        plt.title('Log folding rate distribution', fontsize=25)
        plt.tick_params(axis='both', which='major', labelsize=22, pad=2) 
        
    return Log_folding_distribution

def get_unfolding_distribution(Bootstrapped_EAs, Bootstrapped_intercepts, temp, plot=None):
    """
    At some temperature temp, get a list of log unfolding rates by computing the unfolding rate for each activation energy/intercept pair from the bootstrapped trials
    
    """
    Ntrials=np.shape(Bootstrapped_EAs)[2]
    
    Log_unfolding_distribution=np.zeros(np.shape(Bootstrapped_EAs))  #each page will be the unfolding rates from a given bootstrap trial...enry (i,j) in that page will be the rate for transition i->j in that trial
    
    
    for n in range(Ntrials):
        slopes=-Bootstrapped_EAs[:,:,n]
        intercepts=Bootstrapped_intercepts[:,:,n]
        Log_unfolding_distribution[:,:,n]=intercepts+slopes/temp
    if plot!=None:
        plt.figure()
        plt.hist(Log_unfolding_distribution[plot[0], plot[1], :], bins=int(Ntrials/10), normed=1)
        plt.xlabel('Log unfolding rate', fontsize=25)
        plt.ylabel('Probability', fontsize=25)
        plt.title('Log unfolding rate distribution', fontsize=25)
        plt.tick_params(axis='both', which='major', labelsize=22, pad=2) 
            
    return Log_unfolding_distribution    



def get_unfolding_pathways(labels, PDB_files, temperatures='*', verbose=True, goal_state = -1):
    """
    What are all the unfolding pathways that are taken?
    We ignore loops, so 0->1->2->1->3 becomes 0->1->3
    Sorts unfolding pathways in order of decreasing flux through them.
    Can choose to specify only one temperature, for instnace, '0.925_', (that underscore is important) or a list of them
    
    We only keep pathways that reach a desired unfolded state
    By default, this is the fully unfodled state, and so we indicate this with goal_state = -1
    But if you care about pathways that reach some other state, just specify this other state with goal_state keyword
    """
    
    PDB_files, labels=utils.get_trajectory(labels, PDB_files, temperatures)
    
    times=np.array(utils.get_times(PDB_files))
    zero_times=np.where(times==0)[0]
    pathways=[]
    for n, t in enumerate(zero_times[:-1]):
        traj=labels[t:zero_times[n+1]]
        pathway=[]
        for s,p in enumerate(traj):
            if s==0: #always include the initial time
                pathway.append(p)
            else:
                if p!=traj[s-1]: #a change has occured in trajectory
                    if p in pathway:
                        ind = pathway.index(p)
                        pathway=pathway[0:ind+1]
                    else:
                        pathway.append(p)
                
        if goal_state ==-1 and p==np.max(labels):
            pathways.append(pathway)
        elif goal_state in pathway:
            pathways.append(pathway)  #only keep track of pathways that fully unfold
            
    unique_pathways=[]
    fracs=[]  #what fraction of flux goes through each unique pathway?
    for p in pathways:
        if p not in unique_pathways:
            unique_pathways.append(p)
            fracs.append(len([P for P in pathways if P==p])/len(pathways))
    #sort in order of decreasing flux
    order=np.argsort(-np.array(fracs))
    fracs=np.array(fracs)[order]
    unique_pathways=[unique_pathways[o] for o in order]    

    if verbose:
        print ('Pathway \t \t fraction of flux')
        for p, path in enumerate(unique_pathways):
            if fracs[p]>=0.01: print('{} \t \t {} \n'.format(path, np.round(fracs[p], 3)))
    return pathways, unique_pathways, fracs 


def get_unfolding_transitions(labels, PDB_files, temp='*', verbose=True, min_to_show=0):
    """
    FInds unfolding transitions (ex. (0,1)) and ranks by abundance
    Can choose to specify only one temperature or multiple using wildcard notation
    
    min_to_show tells you the minimum number of instances a transition soudl occur if you want it to be printed by this function
    """
    if type(temp)==int:temp=float(temp)
    if type(temp)!=str: #convert to a format the Get_trajectory can work with
        temp=str(temp)
        while len(temp)<5: temp = '{}0'.format(temp)
        temp='{}_'.format(temp) #need this underscore at the end for Get_trajectory to be happy 
    
    PDB_files, labels=utils.get_trajectory(labels, PDB_files, temp)
    
    times=np.array(utils.get_times(PDB_files))
    zero_times=np.where(times==0)[0]
    transitions=[]
    
    for n, tau in enumerate(zero_times[:-1]):
        traj=labels[tau:zero_times[n+1]]
        for t in range(len(traj)-1):
            if traj[t+1]>traj[t]: #an unfolding transition has occured
                transitions.append((traj[t], traj[t+1]))
    
    unique_transitions=[]
    counts=[]  #
    
    dic = {}
    for p in transitions:
        if p not in unique_transitions:
            unique_transitions.append(p)
            counts.append(len([P for P in transitions if P==p]))
            dic[p] =  counts[-1]
            
    #sort in order of decreasing frequency
    order=np.argsort(-np.array(counts))
    counts=np.array(counts)[order]
    unique_transitions=[unique_transitions[o] for o in order]    

    if verbose:
        print ('Transition \t \t number ofinstances')
        for p, path in enumerate(unique_transitions):
             if counts[p]>min_to_show: print('{} \t \t {} \n'.format(path, counts[p]))   
    
    
    return dic



def infer_folding_rates(clusters, activation_energies, prefactors,  G, temperatures):
    """
    Takes Arrenius parameters and uses detailed balance to compute folding rates
    """
    print('Inferring unknown folding rates from detailed balance...')
    Nclusters = len(clusters)
    folding_rates=np.nan*np.zeros((Nclusters, Nclusters, len(temperatures) ))
    unfolding_rates = np.nan*np.zeros((Nclusters, Nclusters, len(temperatures)))

    for b in range(Nclusters):
        for a in range(Nclusters):
            unfolding_rates[a, b,:] = prefactors[a,b]*np.exp(-activation_energies[a,b]/temperatures)
            for t, temp in enumerate(temperatures):
                if -np.log(unfolding_rates[a,b,t]) < (G[t,b] - G[t,a]):  #barrier height is lower than free energy difference...typically this implies Arrhenius approximation is failing
                    unfolding_rates[a,b,t] = np.exp(-( G[t,b] - G[t,a]) ) #Then we use the barrier height
            folding_rates[b,a,:]= unfolding_rates[a,b,:] * np.exp(G[:,b] - G[:,a])  #detailed balance! 
            
    return folding_rates, unfolding_rates, temperatures

def markov_transmat(temp, data, PDB_files,  sim_labels=[]):
    """
    Builds a transition matrix whose entry i, j is the probability per timestep of transitioning from i to j
    Computed simply by counting the number of times i goes to j, normlaized by the total number of times i goes
    anywhere
    
    If Intermediates_matter=True (default), we treat this like an ordinary Markov transition matrix calculation
    
    
    You also have the otpion of eliminating intermediates in some transitio
    For instnace, if a trajectory proceeds as [1,1,2,2,3,3] and you input Elim_int=(1,3),
    then the trajectory will be converted into [1,1,1,1,3,3]..
    
    That is you don't care about everything that happens between the first instnace of i and 
    the first instance of j--the in-between points get simply converted to i
    
    You also have the option of directly inputting the trajectories yourself via the sim_labels argument, althoguh the 
    default is not to (sim_labels=None)
    This is useful if you are bootstrapping
    
    
    THIS CODE IS DISGUISTING!!! CLEAN IT UP PLEASE!!!
    """
    n_clusters=int(np.max(data))+1  #assumes data is labels
    
    if len(sim_labels)==0:
        times, sim_labels= predict_unfolding_at_temperature(temp, data, PDB_files )
    else:
        times=[]

    #if Elim_int!=None:
    #    sim_labels=Eliminate_intermediates(sim_labels, Elim_int[0], Elim_int[1])

    Ns=np.zeros((n_clusters))
    counts=np.zeros((n_clusters, n_clusters))
    norm=np.zeros((n_clusters, n_clusters))

            
    for traj in sim_labels:
        visited=[]
        for t, l in enumerate(traj):
            if t<len(traj)-1:
                norm[int(l),:]+=1
                counts[int(l), int(traj[t+1])]+=1
                if l not in visited: visited.append(l)
        for v in visited: Ns[int(v)]+=1
    transmat=np.divide(counts, norm)    
    return transmat, times, Ns

    

def plot_cluster_survival(temp, clusters, data, PDB_files, protein,  fit = 'Single', ax = None,  fontsize = 20, labelsize = 20, legend_loc='upper right', min_length = 50, legend = True,colors = ['blue', 'orange', 'green', 'red', 'magenta', 'black','yellow', 'cyan' ], Labels = []):
    """
    Plots probability over time of surviving in each cluser within a provided list clusters
    clusters is a list of cluster numbers for which you want to plot survival probs
    fit can be 'Single' or 'Double' for single or double exponential fit, respectively
    
    """
    if ax == None:
        fig, ax = plt.subplots()

    for nn, i in enumerate(clusters):
        times, labels=predict_unfolding_at_temperature(temp, data, PDB_files,)
        params, survival_probability, t, N=compute_cluster_survival(i,times, labels, min_length = min_length, fit = fit) 
        #print(N)
        
        if fit == 'Single':
            lamb_MLE=params[0]
            print('$\lambda_1$ = {}'.format(lamb_MLE))
            ax.scatter(t, survival_probability, color = colors[nn] )
            if len(Labels)==0:
                label = 'Cluster {}'.format(i)
            else:
                label = Labels[nn]
            ax.plot(t, np.exp(-lamb_MLE*np.array(t)), label=label, color = colors[nn] )
        else:
            [lamb1, lamb2, a1, a2] = params
            print('$\lambda_1$ = {} \n $\lambda_2$ = {}'.format(lamb1, lamb2))
            ax.scatter(t, survival_probability, color = colors[nn] )
            if len(Labels)==0:
                label = 'Cluster {}'.format(i)
            else:
                label = Labels[nn]
            ax.plot(t, a1*np.exp(-lamb1*np.array(t)) + a2*np.exp(-lamb2*np.array(t)), label=label, color = colors[nn] )
            
        
        #plt.title('Survival probability for cluster {}'.format(i), fontsize=32, y=1.03)
        #plt.annotate('T={}'.format(temp), (len(t)/2, 0.6), fontsize=30)
        #exp=np.int(np.floor(np.log10(lamb_MLE)))
        #factor=np.round(lamb_MLE/10**exp, 1)
            #plt.annotate('$\lambda$ = {}*10^{}'.format(factor, exp), (t[int(3*len(t)/4)], np.exp(-lamb_MLE*t[int(3*len(t)/4)] ) ), fontsize=30)
    #ax.yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel('MC Step', fontsize=fontsize)
    ax.set_ylabel('Survival probability', fontsize = fontsize)
    if legend: ax.legend(fontsize=fontsize, loc=legend_loc)
    ax.tick_params(axis='both', which='major', labelsize=labelsize, pad=2) 



def plot_distance_map(distance, key):
#distance = distance*10**(6)  #Originally, the distance was in units of millions of MC steps, but we just want it in raw MC steps
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    labels = []
    
    for item in key:
        string = ''
        b=0
        for bit in item:
            if bit=='1':
                string = '{}{}'.format(string, alphabet[b])
            b+=1
        
        if string == '':
            string = '$\emptyset$'
        labels.append(string)
            
    N = len(key)
     
    for i in range(N):
        for j in range(N):
            d = distance[i,j]
            if d>10000:
                distance[i,j] = np.inf

    
    plt.figure()
    plt.imshow(distance, cmap = 'Reds')
    
    
   
    
    for i in range(N):
        for j in range(N):
            d = distance[i,j]
            if d!=np.inf:
                d = int(d)
                plt.text(i,j, d, horizontalalignment='right') #lol...in matplotlib's world, right actually means left
    
    
    plt.xticks(np.arange(0,N,1),labels=labels, rotation=45)
    plt.yticks(np.arange(0,N,1),labels=labels)
    plt.tick_params(axis=u'both', which=u'both',length=0, labelsize = 20)
    
    plt.xlabel('Topological configuration', fontsize = 20)
    
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=20)
    cb.set_label('Kinetic distance $T_K^{i,j} \ (*10^6$ MC steps)', rotation=270, fontsize = 20, labelpad = 40)
    
    plt.tight_layout()
    





def plot_folding_rates(folding_info_path, titlename, full_protein_transitions, legend_loc='upper right', ylabel = True, ymin=None, labels = None,ax = None, linestyle = '-', errorbars = True,fontsize = 20, colors=['r', 'y', 'g', 'b', 'k', 'm', 'c'], labelsize = 20, legend_fontsize = 20, title = True, temp_norm = 1,PMF_cutoff=10, relabel_trans = None, legend = True):
    """
    Plots folding rate vs temperature for a single protein full_protein_name whose info is in folding_info_path
    You can specify a list of transitions you wish to plot, for instnace [(1,0), (2,1)]
    By default only plots points where free energy difference for relevant transition is less than PMF_cutoff, which is 10 kbT by default


    Finally, you have some options for labeling.
    By default, things are labeled like i->j transition, where (i,j) are the tuples provided in full_protein_transitions
    BUt you have the option to relabel_trans in case the cluster numbering is unintuitive
    By default, this variable is set to None, so no transitions are relabeled
    But you can also feed it a list of tuples, in which case, the ith transition is renamed based on that ith tuple
    If that transition is originally called 3->1 for instance. but the ith tuple is (3,2), then it will be relabled as 3->2
    
    
    Alternatively, you can label that transition whatever you want if you set the varaible labels to some list of 
    desired labels (each of which is string), rather than keeping labels as None (default)
    
    Note, you can set errorbars to False if you don't want errorbars
    
    We only plot a folding rate at a temperature T if the magnitude of the folding free energy at that 
    temperature for the current transition is less than PMF_cutoff

    """
    folding_info = joblib.load(folding_info_path)
    temperatures = folding_info['eq temperatures']
    Bootstrapped_EAs=folding_info['Bootstrapped EAs']
    Bootstrapped_intercepts = folding_info['Bootstrapped intercepts']
    full_folding_rates=folding_info['folding rates']
    free_energies=folding_info['cluster free energies']
    

    if ax==None: fig, ax = plt.subplots()
    
    
    
    for kk, full_protein_transition in enumerate(full_protein_transitions):

        uncertainties=np.zeros((len(temperatures),))
        for j,temp in enumerate(temperatures):
            Log_folding_distribution=get_folding_distribution(Bootstrapped_EAs, Bootstrapped_intercepts, free_energies, temperatures, temp, plot=None)
            uncertainties[j]=np.std(Log_folding_distribution[full_protein_transition[0], full_protein_transition[1],:])
        #indices where we trust the PMF since free energy difference is less than 10 in full protein
        Full_trustable_PMF=np.where(np.abs(free_energies[:,full_protein_transition[1]]-free_energies[:,full_protein_transition[0]])<PMF_cutoff)[0]
        #temperatures=np.array(temperatures)/temp_norm
        
        lower_errors=full_folding_rates[full_protein_transition[0], full_protein_transition[1],:]-full_folding_rates[full_protein_transition[0], full_protein_transition[1],:]*np.exp(-uncertainties)
        upper_errors=full_folding_rates[full_protein_transition[0], full_protein_transition[1],:]*np.exp(uncertainties)-full_folding_rates[full_protein_transition[0], full_protein_transition[1],:]
        
        if relabel_trans !=None:
            iii = relabel_trans[kk][0]
            jjj = relabel_trans[kk][1]
        else:
            iii = full_protein_transition[0]
            jjj = full_protein_transition[1]
            
        if labels == None:
            label = '{} -> {} transition'.format(iii, jjj)
        else:
            label = labels[kk]

        if errorbars: 
            ax.errorbar(temperatures[Full_trustable_PMF]/temp_norm, full_folding_rates[full_protein_transition[0], full_protein_transition[1],Full_trustable_PMF], yerr=[lower_errors[Full_trustable_PMF], upper_errors[Full_trustable_PMF]], color=colors[kk],label=label, marker='*') #changed on 2/10/18 to label according to substructures
            ax.plot(temperatures[Full_trustable_PMF]/temp_norm, full_folding_rates[full_protein_transition[0], full_protein_transition[1],Full_trustable_PMF], color=colors[kk], linestyle = linestyle)
        else:
            ax.plot(temperatures[Full_trustable_PMF]/temp_norm, full_folding_rates[full_protein_transition[0], full_protein_transition[1],Full_trustable_PMF], color=colors[kk], linestyle = linestyle, label = label)
    
    ax.set_yscale('log')
    if legend: ax.legend(fontsize=legend_fontsize, loc=legend_loc)
    if temp_norm != 1:
        ax.set_xlabel('Temperature/$T_M$', fontsize=fontsize)
    else:
        ax.set_xlabel('Temperature', fontsize=fontsize)
    if ylabel: ax.set_ylabel('Inverse MFPT (MC steps$^-1$)', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=labelsize, pad=2)
    if title: ax.set_title(titlename, fontsize=fontsize)
    if ymin!=None:
        ax.set_ylim(ymin=ymin)     


def plot_unfolding_trajectory(traj, labels, PDB_files):
    f, l = utils.get_trajectory(labels, PDB_files, traj)
    #if filter_osc==True: l=filter_oscillations(l, thresh=thresh)
    plt.figure()
    plt.plot(utils.get_times(f),l)
    
    unique=list(set(labels))
    if np.nan in unique: unique.remove(np.nan)
    #plt.ylim((np.min(unique)-0.1, np.max(unique)+0.1))


def predict_unfolding_at_temperature(temp, data, PDB_files):
    """
    Function to predict lables for all trajectoires at a given temperature
    Note: The assumption is that at a given temperature, all snapshots are at the same times
    
    Filter should be 'First commit' or 'Last commit' or 'Filter osc' as described in ClusterPCA
    You can also enter None (or anything else besides the options above) in whcih case no filtering is applied
    """
    
    temp=str(temp)
    
    if len(temp)==1:
        temp='{}.'.format(temp)
    while len(temp)<5:  #add zeros so that the temperature is of the form 0.80
        temp='{}0'.format(temp)
    f, trajectories = utils.get_trajectory(data, PDB_files, '{}_'.format(temp) )
    #need to figure out how long are all the trajectories.
    #to figure this out, iterate through the first files until you see a change

    
    go=True
    i=0
    traj_nums=[]
    while go:
        file=f[i]
        file=file.split('{}_'.format(temp))
        suffix=file[1]
        traj_num=suffix.split('.')[0]
        traj_nums.append(traj_num)
        if traj_nums[i]!=traj_nums[i-1]:
            go=False
        else:
            i+=1
    
    traj_len=i
    n_trajectories=int(len(f)/traj_len)
    sim_labels=np.zeros((n_trajectories, traj_len))
    
    times=utils.get_times(f[0:traj_len])
    for n in range(n_trajectories):  
        traj=trajectories[n*traj_len:n*traj_len+traj_len]  
        sim_labels[n,:]=traj
        
    return times, sim_labels
        

def rates_from_transmat(transmat, dt):
    """
    Given a transition matrix (probability per printed step of making transitions) and a timestep (ex. 1 million MC steps),
    return a rates matrix whose entry i,j is a probablity per continuous time (in the timestep->0) of transitioning from i to j
    
    To derive this relationship, we note that 
    
        lambda_{i to j} = [P(i to j during delta t)/P(i to anywhere but itself during delta t) ]*sum_k lamda_{i to k} (1)
    
    We further note that
    
        transmat[i,i] =P survival during delta t = exp(-sum_k lambda_{i to k}*delta t) (2)
    
    We solve equation (2) for sum_k lambda_{i to k} and plug that into equation (1)
    
    
    """
    n_clusters=np.shape(transmat)[0]
    rates=np.nan*np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        sum_lambdai=1/dt*np.log(1/transmat[i,i])
        for j in range(n_clusters):
            #if i<j:
                rates[i,j]=(sum_lambdai*transmat[i,j])/(1-np.exp(-sum_lambdai*dt))
    #print(time.clock()-t1)
    return rates
                


def runHMM(score_path, f,s, m,starting_states, plot = True):
    """
    Parameters are as follows:
        score_path: Path to unfolding simulation substructure scores
        f: thresh for defining whether a substructure is formed. Should be same value as when computing substrucutre PMF
        s: The Sanders thresh, the minimum representation a topological configuration should have among unfolding snapshots for it to be kept
        m: Misassignment probability used in HMM
        starting_state: topological configuraiton of starting state
    
    By the end, all trajectories are fit to HMM, and we can generate a plot of pairwise kinetic distances between states, if we wish
    """    
    print('Loading scores...')
    scores, PDB_files, Substructures=load_data.load_scores(score_path, f, convert_to_binary=True )
    
    
    labels, key=Bernie_elimination(scores, PDB_files, s, plot=False)
    
    
    alphabet = 'abcdefghijklmnopqrstuvwxyz' #passed Kindergarden
    NN = len(key[0]) 
    
    start_states=[]
    
    for starting_state in starting_states:
        ss = ''
        for n in range(NN):
            if alphabet[n] in starting_state:
                ss= '{}{}'.format(ss, (int(1)))
            else:
                ss= '{}{}'.format(ss, (int(0)))
        
        start_states.append(np.where(np.array(key)==ss)[0][0])
        
    
    start_states = np.array(start_states)
    
    times=np.array(utils.get_times(PDB_files))
    zero_times=np.where(times==0)[0]
    
    #print('Zero times: {}'.format(zero_times))
    
    lengths=[zero_times[t+1]-zero_times[t] for t in range(len(zero_times)-1)]
    lengths.append(200) #need to include the last trajectory
    lengths=np.array(lengths)
    
    
    
    
    for i, t in enumerate(zero_times):
        if t==zero_times[-1]:
            labels[t:]=utils.filter_nans(labels[t:])
        else:
            nextt=zero_times[i+1]
            labels[t:nextt]=utils.filter_nans(labels[t:nextt])
    
    unique_labels=np.unique(labels)
    Nlabels=len(unique_labels)
    
    
    P0_val=1/len(start_states)
    
    start_prob=np.zeros(Nlabels)    
    start_prob[start_states]=P0_val
    
    print('The start_prob array is {} \n'.format(start_prob))
    
    #Prepare the emission and starting probabilities 
    
    emit=((m)/(Nlabels-1) )*np.ones(( Nlabels, Nlabels))
    emit[np.diag_indices(Nlabels)]=1 - m
    
    #start_prob[starting_state]=1
    
    
    #Train the HMM !!!
    
    print('Training HMM...')
    labels=np.array(labels).reshape(-1,1)
    
    HMM = hmmlearn.hmm.MultinomialHMM(n_components = Nlabels, params='t', init_params='t')
    HMM.emissionprob_=emit
    HMM.startprob_=start_prob
    HMM.fit(labels, lengths=lengths)
    S=HMM.predict(labels, lengths=lengths)
    score=HMM.score(labels, lengths=lengths)
    transmat=HMM.transmat_
    
    
    #Construct our distance matrix and do loop clustering
    
    print('Computing kinetic distances...')
    distance=np.nan*np.zeros((Nlabels,Nlabels))
    for i in range(Nlabels):
        for j in range(Nlabels):
            P_i_given_ij=len(np.where(labels==i)[0])/(len(np.where(labels==i)[0])+len(np.where(labels==j)[0]))
            P_j_given_ij=1-P_i_given_ij
            
            distance[i,j]=  P_i_given_ij/transmat[i,j] + P_j_given_ij/transmat[j,i]
            #distance[i,j]=  (transmat[i,j]+transmat[i,i]) *P_i_given_ij/transmat[i,j] + (transmat[j,i]+transmat[j,j]) *P_j_given_ij/transmat[j,i]
            
    distance[np.diag_indices(Nlabels)]=0  #should be zero, but sometimes isn't due to oscillations
    
    if plot: plot_distance_map(distance, key)
    
    return unique_labels, distance, S, key, PDB_files

#    
#
#def Compute_LL(S, clusters, key, PDB_files):
#    """
#    Computes a log likelihood for a given clustering as follows:
#    For every cluster (ex. [101, 111]), assigns emission probabilities in accordance with how often
#    a given configuration appears within the cluster
#    
#    Also computes a transition matrix based on how often you see transitions between each pair of clusters
#    The total log likelihood is the sum of the log probabilities of all observed transitions and emissions
#    
#    Vector S tells you your structure assignment at each timepoint, in numbers (which are converted to configurations via the key)
#    
#    Not using this at the moment
#    """
#    N=len(key)
#    M=len(clusters)
#
#    if type(clusters[0][0])==str:  #converts clusters to a list of  lists of numbers, ex. ([[0,1], [2], [3,4]])
#        C=[]
#        for clust in clusters:
#            c=[]
#            for s in clust: 
#                i = key.index(s)
#                c.append(i)
#            C.append(c)
#    else:
#        C=clusters
#    X=[]  #recreates X, which tells you your cluster assignemt at each timepoint
#    for s in S:
#        clust=[c for c in range(len(C)) if s in C[c]][0]
#        X.append(clust)  
#    X=np.array(X)
#    Emissions = np.zeros(N)  #emissions will be a row vector whose ith entry corresponds to the probability of observing state i given that you are in the cluster to whcih i has been assigned
#    
#    for c, clust in enumerate(C):
#        norm = len(np.where(X==c)[0]) #how many times is cluster observed?
#        for s in clust:
#            count = len(np.where(S==s)[0])
#            Emissions[s]=count/norm
#    
#    E_score = np.sum([np.log(Emissions[s]) for s in S]) #score contribution due to emissions
#
#    counts=np.zeros((M, M))
#    norm=np.zeros((M, M))
#
#    #reshape trajectories so you have one list per trajectory that starts at time 0
#    times=np.array(utils.get_times(PDB_files))
#    zero_times=np.where(times==0)[0]
#    sim_labels=[]
#    
#    for i, t in enumerate(zero_times):
#        if i!=len(zero_times)-1:
#            sim_labels.append(X[t:zero_times[i+1]])
#        else:
#            sim_labels.append(X[t:])
#    
#            
#    for traj in sim_labels:
#        for t, l in enumerate(traj):
#            if t<len(traj)-1:
#                norm[int(l),:]+=1
#                counts[int(l), int(traj[t+1])]+=1
#    transmat=np.divide(counts, norm)
#    
#    T_score=0  #score due to transitions
#    
#    for traj in sim_labels:
#        for t in range(len(traj)-1):
#            transprob = transmat[X[t], X[t+1]]
#            T_score+=(np.log(transprob))
#    return T_score + E_score     