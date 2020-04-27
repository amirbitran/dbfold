#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 12:47:46 2018
TODO: Convert this all to a class!

@author: amirbitran
"""

import numpy as np
import matplotlib.pyplot as plt
import joblib



class Kinetic_model:
    
    def __init__(self, protein_name, lengths, rates=[], connectivity=[], Temperatures = [],Clusters = [],folding_rate_paths=[], tiebreak = 'Less folded', omit_transitions = []):
        """
        Initialize a co-translational folding kinetic model for a given 
        protein_name with folding/unfolding rates at different lengths given by
        list/array lengths
        
        There are one of two ways to initialize the rates:
            
        1. One possibility is you input a list of rate matrices...this is a 
        list of array-types where each array is the rates matrix for a 
        particular length of the protein. Each array can potentially have 
        multiple pages due to multiple temperatures
        Note: If you choose this option, you also need to provide a value for 
        connectivity, wich is a list of lists where element i,j tells you, if 
        the protein is at cluster j at the ith length and you synthesize to the
        i+1st length, to whcih cluster does the protein now go?
        Note, in this case, you must also input the clusters! And the temperatures!
        
        2. Alternatively, you can input a list of paths, each of which 
        contains a file with folding info. The info in these files 
        is then used to construct the rates matrix.
        Note: If you choose this option, you can either specify connectivity
        yourself as in option 1, or leave that variable blank and the
        algorith will infer connectivity by finding pairs of clusters at 
        different length with maximal structural similarity
        
        The optional argument tiebreak involves computing connectivity
        If cluster i from length L is equally similar to clusters j and k
        at length L+1, then if tiebreak = 'Less folded', you connect to 
        whichever of j or k is less folded
        
        Otherwise, connect to the more folded one
        
        omit_transitions is a list of duplets for transitions you want to keep
        out of kinetic model, for instnace if folding rate prediction is poor
        For instnace, you may set omit_transitions = [(3,0)],
        in which case transitions from 3 to 0, and vice versa are both omitted 
        
        
        """
        self.protein_name=protein_name
        self.lengths=lengths
        min_len=min(lengths)
        Nclusters=[] #A list, the ith element tells you how many clusters the ith length has
        
        #Clusters=[]  #A list of lists, each of which contains the clusters for a given length
        
        if len(connectivity)==0: self.Infer_connectivity=True
        
        if len(rates)!=0:  #we use the provided list of folding rates -- just convert to arrays first
            for r in range(len(rates)):
                rates[r]=np.atleast_3d(np.array(rates[r]))  #we make the arrays 3 dimensional (even if therr is only one page per array) to be consistent with the possibitliy that there can be multiple temperatures
            #self.rates = rates
            if len(connectivity)==0 and len(rates)>1:
                print('Error! Connectivity must be specified if you provide a custom rate matrix with multiple lengths')
            #Temperatures=[]
            #self.Clusters = Clusters
            for n in range(len(rates)):
                Nclusters.append(len(Clusters[n]))
            
            
        elif len(folding_rate_paths)!=0: #we construct the rate matrix using information in the directories
            #print('moo')
            Temperatures = [] #A list of lists, the ith element tells you the temperatures at the ith length
            for n, path in enumerate(folding_rate_paths):
                folding_info = joblib.load(path)
                folding_rates = folding_info['folding rates']
                unfolding_rates = folding_info['unfolding rates']
                clusters = folding_info['clusters']
                temperatures = folding_info['eq temperatures']

                    
                
                transmat = folding_rates.transpose((1,0,2))*lengths[n]/min_len +unfolding_rates.transpose((1,0,2))*lengths[n]/min_len  #all rates are reweighed based on the shortest protein, essentially converting to monte carlo "sweeps"
                Temperatures.append(temperatures)
                
                if len(omit_transitions)>0:
                    for duple in omit_transitions:
                        i = duple[0]
                        j = duple[1]
                        transmat[i,j]=0
                        transmat[j,i]=0
                
                #Note, we transpose the folding_rates and unfolding_rates matrix because in these matrices, folding_rates[i,j] means rate of transitioning tfom i to j
                #But in a typical transition matrix, element [i,j] actually means the rate of going from j to i
                #Enforce rows need to sum to 0
                for t in range(np.shape(transmat)[2]):
                     for i in range(np.shape(transmat)[0]):
                         transmat[i,i,t] = -np.sum(transmat[:, i, t])
                rates.append(transmat)
                Nclusters.append(len(clusters))                  
                
                Clusters.append(clusters)
        else:
            print('Error! Need to either specify a rate matrix, or input a list of directories with folding rate information')
    
        self.folding_rate_paths=folding_rate_paths
        self.clusters = Clusters
        self.Nclusters = Nclusters
        self.tiebreak = tiebreak
        if self.Infer_connectivity: 
            self.compute_connectivity(self.tiebreak)
        else:
            self.connectivity=connectivity
        self.rates=rates
        self.Temperatures=Temperatures
        self.Cleared_extraneous=False
        

    def compute_connectivity(self, tiebreak):
        connectivity=[]
        for n, direc in enumerate(self.folding_rate_paths): 
            if n>0: 
                conn=[]
                for i, cl in enumerate(old_clusters):
                    mean_distances=[]  #mean distance from cluster i to each of the clusters j
                    for j, cu in enumerate(self.clusters[n]): #now, loop through all states in each cluster
                        dis=[]
                        for k in cl:
                            for l in cu:
                                dis.append( np.sum([ np.abs(int(k[vv]) - int(l[vv])) for vv in range(len(k))  ]))
                        mean_distances.append(np.mean(dis)) 
                    winners = np.where(mean_distances ==np.min(mean_distances))[0] #to which clusters in new length is cluster i from old length closest to
                    
                    if tiebreak =='Less folded':
                        conn.append(winners[-1]) #if there is a tie, connect to the less folded cluster, to be conservative
                    else:
                        conn.append(winners[0]) 
                connectivity.append(conn)
            old_clusters=self.clusters[n]
        self.connectivity=connectivity
    
    def delete_extraneous_clusters(self):
        """
        Deletes clusters that are completely isolated, and adjusts transition matrix accordingly
        Learns this based on transition matrix at each length from the temperature whose index is given 
        by temp_index
        
        
        
        By default, temp_index is 0, so it learns the isolated cluster from the lowest temperature 
        
        Prior to 2/11/19: It learned at the lowest temperature at which transition matrix
        has no Nans, since Nans are an indicator that the PMF was not trusted at given temperature,
        so folding rates there aren't meaningful...but then I realized thi sisn't true since
        Nans are implemented after the fact by Compare_folding_Rates etc..
        moreover, there was an issue with CMK that caused me to just default to temp_idnex = 4 (T = 0.5)
        I did not use 0 since for some reason, MarR had Nans at temp index of 0...not sure why
        """
        rates=[]
        clusters=[]
        temp_index = 4 
        for n in range(len(self.lengths)):
            
            #Deleted the following on 2/12/19 due to a weird case with CMk...see notes
            #temp_index = 0 #find lowest temperature at whcih the transition matrix has no Nans
            #while np.any(  np.isnan(self.rates[n][:,:,temp_index])):
            #    temp_index+=1
            
            TT = self.rates[n][:,:,temp_index]
            keep=[]
            for r in range(len(TT)):
                if np.max(np.abs(TT[:,r]))>0 or np.max(np.abs(TT[r,:]))>0 : #not completely isolated
                    keep.append(r)
                else:
                    print('Length {}: Deleting cluster {}'.format(self.lengths[n],self.clusters[n][r]))
            keep = np.array(keep)
            #keep only rows and columns for connected clusters...first keep rows, then keep columns
            zz = self.rates[n][keep, :, :]
            zz = zz[:, keep, :]
            rates.append(zz)
            self.Nclusters[n] = np.shape(zz)[0]
            CC = [clust for c, clust in enumerate(self.clusters[n]) if c in keep ]
            clusters.append(CC)
            print('\t')
            
        self.rates = rates
        self.clusters=clusters
        self.Cleared_extraneous=True
            
        if self.Infer_connectivity: 
            self.compute_connectivity(self.tiebreak)
   
    
    def diagonalize(self):
        """
        Diagonlize transition matrices

        """
        self.eigenvectors=[]  #list of matrices containing eigenvectors for respective lengths
        self.eigenvalues=[]
        
        if not self.Cleared_extraneous: self.delete_extraneous_clusters()  
        
        for n, transmat in enumerate(self.rates):
            v=np.zeros(np.shape(transmat))
            w = np.zeros((np.shape(transmat)[0], np.shape(transmat)[2]))
            
            for t in range(np.shape(transmat)[2]): #diagonalize each temperature (page) separately
                TT = transmat[:,:,t]
                
                #min_element = np.min(np.abs(TT[np.where(TT)]))
                #w[:,t], v[:,:,t] = np.linalg.eig(TT/min_element)
                if not np.any(np.isnan(TT)): #diagnoalize if there are no Nans, so all rates have been computed
                    w[:,t], v[:,:,t] = np.linalg.eig(TT)
                else:
                    w[:,t] = np.nan
                    v[:,:,t] = np.nan
                #w[:,t] = w[:,t]*min_element
                #Ocasionally, eigenvalues end up positive due to numerical mistake...if they are positive and super small (like less than 10**-20), we just convert them to 0
                #But if they are positive and appreciable, we print a warning
                #if np.max(w[:,t])>10**(-20):
                #    print('Warning: Found a positive eigenvalue of magnitude {} at temperature {}'.format(np.max(w[:,t]), self.Temperatures[t]))
            w[np.where(w>0)]=0
            self.eigenvalues.append(w)
            self.eigenvectors.append(v)

        
    def MFPT(self, temp, Trans_times, folded_cluster=0, starting_state = -1):
        """
        Compute MFPT to folded, fully synthesized state
        Input a list of times that tells you how long you spend at each length
        If you have L lengths, then this list of times shoudl only be L-1 long
        
        By default, you are assumd to start in the last cluster, which is generally the least folded one
        But for proteins such as MarR and CMK, you do not get reliable rates of transition from last cluster
        to penultimate one, since this involves folding of a simple antiparallel beta hairpin, which algorithm struggles with
        
        Thus, it may be worthwhile to set starting_state = -2 for these
        """

        self.delete_extraneous_clusters() 
        if ~hasattr(self, 'eigenvectors'):
            self.diagonalize() 
       
        #print('Length of rates list is {}'.format(len(self.rates)))
        for n, transmat in enumerate(self.rates):
            #print(n)
            temperatures=self.Temperatures[n]
            temp_index = np.where(temperatures==temp)[0][0]
            if n==0:   #At the initial length, you start in the unfolded state
                #print('baahhh')
                P0 = np.zeros(self.Nclusters[n])
                P0[starting_state]=1
            else: #propagate probability distribution from last timepoint of previous length to the corresponding clusters in the new length
                #print('Kawkawwwww')
                P0_old=P  
                P0=np.zeros(self.Nclusters[n])
                for i, p in enumerate(P0_old):
                    P0[self.connectivity[n-1][i]]+=p
            
            if n<len(self.rates)-1:  #Unless you are in the final length, we compute the time evolution normally
                #print('Broooo')
                tau = Trans_times[n]
                lambdas = self.eigenvalues[n][:, temp_index]
                lambdatau=np.array(np.dot(lambdas, tau), dtype=np.float32)  #convert to float 32 helps w numerics? 
                exp_lambdat = np.diag(np.exp(lambdatau))
                v = self.eigenvectors[n][:,:, temp_index]
                #Compute matrix exponential (time evolution operator)
                M = np.linalg.multi_dot((v, exp_lambdat, np.linalg.inv(v)))
                
                #Now time evolve the initial distribution
                P = np.dot(M,P0)
                #print(P)
            else: #For the final length, we set an absorbing boundary at the folded cluster so that we can compute a MFPT
                #print('MOO')
                MM=transmat[folded_cluster+1:,folded_cluster+1:,temp_index] #get rid of row/column corresponding to folded cluster, so that it becomes an absorbing state outside the system
                Nstates = np.shape(MM)[0]                
                MFPTs = np.dot( np.linalg.inv(MM.transpose()), -1*np.ones((Nstates, 1)))  #MFPT as a function of starting state
                MFPT = np.dot(P0[folded_cluster+1:], MFPTs)[0] #integrate over all states the product of the probability of starting in that state * MFPT given that state
                
        return MFPT + np.sum(Trans_times)



    def MFPT_and_yield(self, temp, Trans_times,folded_clusters=[0], starting_state = -1):
        """
        Compute MFPT to folded, fully synthesized state, as well as how much protein ends up in folded state
        after translation
        Input a list of times that tells you how long you spend at each length
        If you have L lengths, then this list of times should be L long,
        since the last element of this list tells you how long you go 
        between the time where you acquire folding properties of full protein 
        and end of translation (ex. for MarR, it would be between 110 and 144)
        
        By default, you are assumd to start in the last cluster, which is generally the least folded one
        But for proteins such as MarR and CMK, you do not get reliable rates of transition from last cluster
        to penultimate one, since this involves folding of a simple antiparallel beta hairpin, which algorithm struggles with
        
        Thus, it may be worthwhile to set starting_state = -2 for these
        
        folded_cluster is a list of clusters that count as folded...to compute yield, you compute
        sum of probaiblites of beign in one of these
        For MFPT, we set absorbing boundary at first one listed
        """
        folded_clusters = np.array(folded_clusters)        

        self.delete_extraneous_clusters() 
        if ~hasattr(self, 'eigenvectors'):
            self.diagonalize() 
       
        #print('Length of rates list is {}'.format(len(self.rates)))
        for n, transmat in enumerate(self.rates):
            #print(n)
            temperatures=self.Temperatures[n]
            temp_index = np.where(temperatures==temp)[0][0]
            if n==0:   #At the initial length, you start in the unfolded state
                #print('baahhh')
                P0 = np.zeros(self.Nclusters[n])
                P0[starting_state]=1
            else: #propagate probability distribution from last timepoint of previous length to the corresponding clusters in the new length
                #print('Kawkawwwww')
                P0_old=P  
                P0=np.zeros(self.Nclusters[n])
                for i, p in enumerate(P0_old):
                    P0[self.connectivity[n-1][i]]+=p
            
                #print('Broooo')
            tau = Trans_times[n]
            lambdas = self.eigenvalues[n][:, temp_index]
            lambdatau=np.array(np.dot(lambdas, tau), dtype=np.float32)  #convert to float 32 helps w numerics? 
            exp_lambdat = np.diag(np.exp(lambdatau))
            v = self.eigenvectors[n][:,:, temp_index]
            #Compute matrix exponential (time evolution operator)
            M = np.linalg.multi_dot((v, exp_lambdat, np.linalg.inv(v)))
            
            #Now time evolve the initial distribution
            P = np.dot(M,P0)
                #print(P)
            if n==len(self.rates)-1: #For the final length, we set an absorbing boundary at the folded cluster so that we can compute a MFPT
                #print('MOO')
                MM=transmat[folded_clusters[0]+1:,folded_clusters[0]+1:,temp_index] #get rid of row/column corresponding to folded cluster, so that it becomes an absorbing state outside the system
                Nstates = np.shape(MM)[0]                
                MFPTs = np.dot( np.linalg.inv(MM.transpose()), -1*np.ones((Nstates, 1)))  #MFPT as a function of starting state
                MFPT = np.dot(P0[folded_clusters[0]+1:], MFPTs)[0] #integrate over all states the product of the probability of starting in that state * MFPT given that state
                
        return MFPT + np.sum(Trans_times), np.sum(P[folded_clusters])
    
    def Yield(self, temp, Trans_times, folded_cluster=0, starting_state = -1):
        """
        Compute amount of folded protein at the end of synthesis
        Trans_times should include time to go from length i to i+1, as well
        as the last entry shoudl be the time from final length to termination of synthesis
        """

        self.delete_extraneous_clusters() 
        if ~hasattr(self, 'eigenvectors'):
            self.diagonalize() 
       
        #print('Length of rates list is {}'.format(len(self.rates)))
        for n, transmat in enumerate(self.rates):
            #print(n)
            temperatures=self.Temperatures[n]
            temp_index = np.where(temperatures==temp)[0][0]
            if n==0:   #At the initial length, you start in the unfolded state
                #print('baahhh')
                P0 = np.zeros(self.Nclusters[n])
                P0[starting_state]=1
            else: #propagate probability distribution from last timepoint of previous length to the corresponding clusters in the new length
                #print('Kawkawwwww')
                P0_old=P  
                P0=np.zeros(self.Nclusters[n])
                for i, p in enumerate(P0_old):
                    P0[self.connectivity[n-1][i]]+=p
                #print('Broooo')
            tau = Trans_times[n]
            lambdas = self.eigenvalues[n][:, temp_index]
            lambdatau=np.array(np.dot(lambdas, tau), dtype=np.float32)  #convert to float 32 helps w numerics? 
            exp_lambdat = np.diag(np.exp(lambdatau))
            v = self.eigenvectors[n][:,:, temp_index]
            #Compute matrix exponential (time evolution operator)
            M = np.linalg.multi_dot((v, exp_lambdat, np.linalg.inv(v)))
            
            #Now time evolve the  distribution
            P = np.dot(M,P0)

        return P[folded_cluster]
 
    
def compute_time_evolution(t, eigenvalues, eigenvectors, C, state=3):
    """
    Solves master equation and returns the time evolution of a desired state
    Default is 3, which is completely synthesized, folded state
    """
    if type(t)!=np.ndarray:
        t = np.array([t])
        
    lambda_t=np.dot(eigenvalues.reshape(len(eigenvalues),1), t.reshape(1, len(t))  )  #creates a matrix where horizontal dimension represents different time points in t. In each row of this matrix, the timepoints are all multiplied by a different eigenvalue
    exp_lambda_t=np.exp(lambda_t)   #exponential of the above matrix    

    coeff_times_exp=np.zeros(np.shape(exp_lambda_t))  #multiply each row of the above vector by the coefficient corresponding to it
    for n, c in enumerate(C):
        coeff_times_exp[n,:] = c*exp_lambda_t[n,:]

    #at a given column t, the row i of coeff_times_exp corresponds to  coefficient i multiplied by the exponential of the ith eigenvalue times time (t)

    #So now, to get the time evolution, we sipmly take the matrix product of our eigenvector matrix with coeff_times_exp

    time_evolution=np.dot(eigenvectors,coeff_times_exp)
    time_evolution=time_evolution[state,:] #time evolution of folded state
    return time_evolution

def plot_timecourse(X, temp, Trans_times, clusters_to_plot, starting_state = -1, final_time = None, colors= ['r','b', 'g', 'y', 'k'], linewidth = 8,colorchange = (None, None, None), scientific = True, ylim = (), ntimes = 200, ax = None,fontsize=30, labelsize = 25, timenorm = 1, ylabel = True, xlabel= True):
    """
    Compute amount of folded protein at the end of synthesis
    Trans_times should include time to go from length i to i+1, as well
    as the last entry shoudl be the time from final length to termination of synthesis
    

    
    clusters_to_plot is a list of lists
    the nth sublist tells you which cluster numbers you want to follow time evolution of 
    at the nth length
    
    But you may want to follow the time evolution of the summed probaiblity of being in
    either of two clusters, in which case some element of the sublist can itself be a list
    i
    
    ntimes is the number of timepoints at which you compute probabilities at each length--200 by default
    
    An important note involves color continuity...Note that the ith element of all
    of the sublists of clusters_to_plot will be plotted with the same color, 
    so make sure to order the list clusters_to_plot accordingly (i.e. probably
    you'll want to do it based on connectivity)
    

    
    You can also make times on the x axis dimensionless by dividing by the slowest folding time, given by timenorm
    By default, timenorm = 1, in which case you do not normalize 
    But if timenorm is set to something other than 1 (assumed to be slowest folding time), then everything is normalized by that slowest folding time
    
    By default, stops plotting once synthesis is complete. But perhaps you want to keep following the timecourses after
    synthesis completes. In that case, set final_time to something other than None. In that case, you will keep following time evolution
    for times up through final_time, assuming same folding kinetics as in the final lenght regime. 
    Note that final_time shold be inputted in normalized units, so if timenorm is set to something other than 1, that final time shoudl be divided by timenorm
    Note that final_time should have a value GREATER than the time value after synthesis ends.
    If it's less, the code won't do anything
    
    scientific indicates whether you want x axis in scientific notation
    
    A new fun feature: You can now indicate a COLORCHANGE! That is, for example in MarR, you may want unfolded cluster 1 (with only hairpin folded)
    to appear gold until length regime 1, at which point it will appear red. Then you can enter the keyword colorchange = (1, 1,'red') beyond which point,
    cluster 1 will appear red
    First element of tuple is length at whcih change should occur, second element is cluster that is to be affected by colorchange, third element is the new color
    
    """
    #TODO: change so that if x label has like 10*(11), make that exponent appear with the same fontsize as the other text!
    if ax == None: fig, ax = plt.subplots()
    X.delete_extraneous_clusters() 
    if ~hasattr(X, 'eigenvectors'):
        X.diagonalize() 
   
    #plt.figure()
    
    
    for n, transmat in enumerate(X.rates): #loop through lengths
        #print(n)
        if colorchange[0]!=None and colorchange[0]==n:
            colors[colorchange[1]]=colorchange[2]
            
        temperatures=X.Temperatures[n]
        temp_index = np.where(temperatures==temp)[0][0]
        if n==0:   #At the initial length, you start in the unfolded state
            #print('baahhh')
            P0 = np.zeros(X.Nclusters[n])
            P0[starting_state]=1
            
            if Trans_times[n]==0:
                taus = np.array([0])
            else:
                taus = np.arange(0, Trans_times[n],(Trans_times[n])/ntimes) #times at which you plot for this given length
            taus_cum = taus
        else: #propagate probability distribution from last timepoint of previous length to the corresponding clusters in the new length
            #print('Kawkawwwww')
            P0=np.zeros(X.Nclusters[n])
            for i, p in enumerate(P_final):
                P0[X.connectivity[n-1][i]]+=p 
            #taus = np.arange(0, Trans_times[n],(Trans_times[n])/ntimes)
            
            if Trans_times[n]==0:
                taus = np.array([0])
            else:
                taus = np.linspace(0, Trans_times[n],ntimes)
            taus_cum = taus+np.sum(Trans_times[0:n])

        
        lambdas = X.eigenvalues[n][:, temp_index]
        v = np.array(X.eigenvectors[n][:,:,temp_index], dtype=np.float32)
        C=np.linalg.solve(v, P0)
        
        for z, clust in enumerate(clusters_to_plot[n]):
            if type(clust)==list:
                Pi = np.zeros((len(clust), ntimes))
                for r, cc in enumerate(clust):
                    Pi[r, :] =  compute_time_evolution(taus, lambdas, v, C, state=cc)
                P_t = np.sum(Pi, axis = 0)
            else:
                P_t = compute_time_evolution(taus, lambdas, v, C, state=clust)
                
            ax.plot(taus_cum/timenorm, P_t, color=colors[z], linewidth = linewidth)

        #Figure out probability of each state at final time
        
        P_final = np.zeros(np.shape(P0))
        for k in range(len(P_final)):
            P_final[k] = compute_time_evolution(Trans_times[n], lambdas, v, C, state=k)
        
        if n!=0 and Trans_times[n]!=0:
            ax.axvline(x=taus_cum[0]/timenorm, color = 'k', linestyle = ':', linewidth = linewidth/2)  #draw a vertical line at the translation time
        
    
    #Keep plotting beyond synthesis, if desired
    if final_time!=None and final_time >np.max(taus_cum/timenorm):
        taus = np.linspace(0, final_time*timenorm - np.max(taus_cum),ntimes)
        taus_cum = taus + np.max(taus_cum)
        ax.axvline(x=taus_cum[0]/timenorm, color = 'k', linestyle = ':', linewidth = linewidth/2)
        P0 = P_final
        C=np.linalg.solve(v, P0)
        for z, clust in enumerate(clusters_to_plot[-1]): #plot the same clusters as at the final length regime
            if type(clust)==list:
                Pi = np.zeros((len(clust), ntimes))
                for r, cc in enumerate(clust):
                    Pi[r, :] =  compute_time_evolution(taus, lambdas, v, C, state=cc)
                P_t = np.sum(Pi, axis = 0)
            else:
                P_t = compute_time_evolution(taus, lambdas, v, C, state=clust)
                
            ax.plot(taus_cum/timenorm, P_t, color=colors[z], linewidth = linewidth)
   
    
    ax.tick_params(labelsize = labelsize)
    if scientific: ax.ticklabel_format(axis = 'x', style = 'scientific', scilimits = (0,0))
    if timenorm ==1 and xlabel:
        ax.set_xlabel('Time (MC sweeps)', fontsize=fontsize, labelpad = 20)
    elif timenorm!=1 and xlabel:
        ax.set_xlabel('Time / Slowest folding time', fontsize=fontsize, labelpad = 20)
    if ylabel: ax.set_ylabel('Probability', fontsize=fontsize, labelpad = 20)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    
    

def linemap( T_pause, T_posttrans, Total_transtimes_rare, Total_transtimes_opt,  rare_slowdown_factors, slowest_folding_rate, slowdowns_to_plot = [1,3,6,9],ax = None, legend = False, label = True,labelsize = 35, fontsize = 40, slowdown_labelsize = 35, ylim = None):
    """
    As a function the folding/synthesis rate, plots ratio of post-translational folding time to co-translatiaonl folding time for vairous possible rare codon slowdowns
    """
    if ax ==None: fig, ax = plt.subplots()
    #linestyles = ['-','--', '-.', ':']
    colors = ['rebeccapurple', 'darkorchid', 'darkviolet','purple', 'mediumorchid', 'violet', 'fuchsia', 'magenta'  ]
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':' ]
    for zz, slowdown in enumerate(slowdowns_to_plot):
        ind = np.where(rare_slowdown_factors==slowdown)[0][0]
        xxx = Total_transtimes_opt[:, ind]*slowest_folding_rate
        #label_value = np.round((Total_transtimes_rare[0,ind]/Total_transtimes_opt[0,ind] - 1)*100, 2)
        yy = T_posttrans/T_pause[:,ind]
        ax.plot(xxx, yy, linewidth = 4,label = 'Rare slowdown factor = {}'.format(slowdown-1), color = colors[zz], linestyle = linestyles[zz])
        #index = int(7*len(xxx)/8) #index for position at which text will be added
        index = np.argmax(yy)
        if label:
            if zz ==len(slowdowns_to_plot) - 1:
                ax.text( xxx[index] - 0.002, T_posttrans/T_pause[index, ind]+0.8, '{}% slowdown'.format( (slowdown-1)*100  ), fontsize = slowdown_labelsize , color = colors[zz]  )
            elif slowdown ==1:
                ax.text( xxx[index] - 0.002, T_posttrans/T_pause[index, ind]+0.8, 'No slowdown', fontsize = slowdown_labelsize , color = colors[zz]  )
        #elif zz ==1:
        #    ax.text( xxx[index] - 0.002, T_posttrans/T_pause[index, ind]+0.8, '{}% slowdown'.format( (slowdown-1)*100  ), fontsize = slowdown_labelsize , color = colors[zz]  )
            else:
                ax.text( xxx[index] + 0.002, T_posttrans/T_pause[index, ind]+0.0, '{}%'.format( (slowdown-1)*100  ), fontsize = slowdown_labelsize , color = colors[zz]  )
    ax.plot(xxx, [1 for i in range(len(xxx))], linestyle = ':', color = 'k')
    ax.set_xlabel('Slowest folding rate/Synthesis rate', fontsize = fontsize, labelpad = 20)
    ax.set_ylabel("τ$_{post-translational}$ / τ$_{co-translational} $", fontsize = fontsize+5, labelpad = 20)
    if ylim!=None:
        ax.set_ylim(ylim)
    ax.tick_params(labelsize = labelsize)
    if legend: ax.legend(fontsize = fontsize)
    plt.show()        

 #clusters_to_plot=[ [ [0,1]],  [4, 3, 0], [5, 4, [1,2], 0 ] ]

    