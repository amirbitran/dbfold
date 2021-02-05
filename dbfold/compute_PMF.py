#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:20:56 2020

@author: amirbitran


REFERENCES
[1] Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states.
J. Chem. Phys. 129:124105, 2008
http://dx.doi.org/10.1063/1.2978177
"""


import numpy as np
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
#import os
#import pickle
#import os.path
import joblib
import matplotlib.pyplot as plt
#import itertools
#import math
#import scipy
#from scipy import signal
#import copy as cp
import dbfold.analyze_structures as analyze_structures
import dbfold.load_data as load_data
import dbfold.utils as utils



def cluster_free_energies( clusters, x, free_energies, temperatures):
    """
    Compute log sum exp of free energies for configurations belong to each cluster, (obtained from analyzing unfolding simulations for fast kinetic exchnage) 
    and if desired, for each cluster, we also compute change in free energy between "trapped" and "non-trapped" configurations
    
    clusters is a list of clusters
    x is list of all unique topoological configurations, and free_energies is array of respective free energies, in same order (along second dimension)
    
    At each temperature, free energies are shifted such that lowest free energy is 0
    
    """
    print('Computing cluster free energies...')
    Nclusters=len(clusters)
    G=np.zeros((len(temperatures), Nclusters ))  #free energy of clusters
    #N0=min([len(a) for a in x]) #how many NATIVE substructures are present in the structure
    
    if len(np.shape(free_energies)) > 2:
        free_energies = free_energies[:,:,0] #For now, ignore that seocnd page, which is the uncertainties...may be useful to incorporate it later
    
    #Compute free energy of clusters by computing ln of partial trace over topological configurations in each cluster
    for c, clust in enumerate(clusters):
        Z=np.zeros(len(temperatures))  #partition function for each temperature
        for top in clust:
            if top in x:  #if not, then this is a topological configuration taht wasn't observed in replica exchanges...assume its PMF is vanishingly high at temperatures we care about
                index=x.index(top)
                Z+=np.exp(-free_energies[:, index])
        G[:, c]=-np.log(Z)
    
    for T in range(len(temperatures)):
        G[T,:] = G[T,:] - np.min(G[T,:])
    return G


def plot_1D_PMF(paths, temps_to_plot,  xlabel='Number of native contacts',legend = True, legend_loc='upper left', xlim=(), upper_cutoff=25, fontsize = 28,  labelsize = 28,title = True, errorbar = False):    
    """
    Plots PMF vs number of native contacts (or some other variable) in 1 dimension
    
    You can set errorbar to True, in which case the errors on the PMF are plotted,but then he variable
    free_energy better have a second page corresponding to those uncertainties! (i.e. shuold have computed it more recently once this feature was implemented)
    """
    plt.figure()
    colors = ['r', 'b', 'g', 'm', 'k', 'y' ,'c']
    for ppp, path in enumerate(paths):
        PMF_info = joblib.load(path)
        bin_centers = PMF_info['native values']
        free_energies = PMF_info['native free energies']
        temperatures = PMF_info['eq temps']
        if errorbar and len(np.shape(free_energies)) <3: print('WARNING!!!: uncertainties in free energies were not computed or not saved for path {}!'.format(path))
        if len(np.shape(free_energies))>2:
            uncertainties = free_energies[:,:,1]
            free_energies = free_energies[:,:,0]
        temperatures=np.unique(temperatures)
        indices_to_plot=[t for t in range(len(temperatures)) if temperatures[t] in temps_to_plot]
        for t in indices_to_plot:                
            plt.plot(bin_centers, free_energies[t], label=str('{}, T ={}'.format(path, temperatures[t])), color = colors[ppp] )
            if errorbar:
                indz = np.where(free_energies[t]<upper_cutoff)[0]
                plt.errorbar(bin_centers[indz], free_energies[t, indz], yerr = uncertainties[t, indz], fmt = 'none', ecolor = colors[ppp]) #fmt = 'none' to imply we don't want a line connecting points, but we still specify the ecolor to imply we want the errorbars there
        plt.ylim(0,upper_cutoff)
        if len(xlim)!=0: plt.xlim(xlim[0], xlim[1])
        plt.tick_params(labelsize=labelsize)
        #plt.ylim(min(free_energies[str(temp)]), min(free_energies[str(temp)])+25)
        
        if title: plt.title("Potential of mean force", fontsize=fontsize)
        plt.xlabel(xlabel, fontsize=fontsize, labelpad = 20)
        plt.ylabel('PMF ($k_{B}T$)', fontsize=fontsize, labelpad = 20)
    if legend: plt.legend(fontsize=28, loc=legend_loc)


def plot_natives_vs_temp(directories, proteins, PMF_infos = None,norm=None, ax = None, legend = True, d_cutoff = 6.5, min_seq_separation = 8, Native_structures = [],fontsize = 30, labelsize = 25, legend_fontsize = 25, labelpad = 20,colors = ['b', 'g', 'm', 'k'], linewidth = 8, xlabel = None, temp_norm = 1, filenames = []):
    """
    Keyword norm can have three modes:
        1. None: No normalization
        2. Fraction: Divide by number of contacts in native structure, to convert to fraction of native contacts..in this case, optional variable Native_strucutres must be entered with paths to respective native structures
        3. Monomer: Assumes we are plotting monomer and dimer melting curves: Divides dimer (second argument) native contacts by 2, to yield Native contacts per monomer,
            then normalizes everything by maximum value for dimer
    
    if temp_norm is some value other than 1, then you divide all temperatueres by that value (i.e. to normalize by melting temperature)
    
    If you wish, you may bypass entering directories (set this variable to an empty list []) if you instead provide your own PMF_info directory (or a list of PMF_info directories) via the optional argument
    """
    if ax == None: fig, ax = plt.subplots()
    if type(directories)==str: 
        directories = [directories]
    
    if len(directories)==0:
        directories = ['' for i in range(len(PMF_infos))]
    
    if type(proteins)==str: proteins = [proteins]
    mean_contacts = [[] for d in range(len(directories))]
    list_of_temps = []
    for d, directory in enumerate(directories):
        if len(directory)!=0:        
            if len(filenames)==0:
                filename = 'PMFs.dat'
            else:
                filename = filenames[d]
            PMF_info=joblib.load('{}/{}'.format(directory, filename))
        else:
            PMF_info = PMF_infos[d]
        free_energies = PMF_info['native free energies']
        natives = PMF_info['native values']
        temperatures = PMF_info['eq temps']
        #print(free_energies[20,8,0])
        if len(np.shape(free_energies))>2: free_energies = free_energies[:,:,0]
        if norm =='Fraction': 
            Ntotal = analyze_structures.count_contacts(Native_structures[d], d_cutoff, min_seq_separation)
            natives=natives/Ntotal
        elif norm == 'Monomer' and d == 1:
            natives = natives/2
        #temperatures=np.unique(temperatures)
        

        for t, temp in enumerate(temperatures):
            Z=np.sum(np.exp(-free_energies[t,:]))
            mean_contacts[d].append(1/Z*np.dot(np.exp(-free_energies[t,:]), natives))
        mean_contacts[d] = np.array(mean_contacts[d])
        list_of_temps.append(temperatures)
    
    if norm =='Monomer':
        for d, directory in enumerate(directories):
            mean_contacts[d] = mean_contacts[d]/np.max(mean_contacts[1])
    
    for d, directory in enumerate(directories):
        ax.plot(np.array(list_of_temps[d])/temp_norm, mean_contacts[d], label=proteins[d], linewidth = linewidth, color = colors[d])
    if norm =='Fraction':
        ax.set_ylabel('Fraction of native contacts', fontsize=fontsize, labelpad = labelpad)
    elif norm == 'Monomer':
        ax.set_ylabel('Fraction of native contacts per subunit', fontsize=fontsize, labelpad = labelpad)
    else:
        ax.set_ylabel('Mean native contacts', fontsize=fontsize, labelpad = labelpad)
    
    if xlabel == None:
        ax.set_xlabel('Temperature', fontsize=fontsize, labelpad = labelpad)
    else:
        ax.set_xlabel(xlabel, fontsize=fontsize, labelpad = 20)
        
    ax.tick_params(labelsize=labelsize)
    if legend: ax.legend(fontsize=legend_fontsize, loc='upper right')  




def plot_pfold_vs_length(protein, directories, fileroots, lengths, subs, temp, fontsize = 32, labelpad = 20,labelsize = 26,d_cutoff = 6.0, min_len = 0, linestyle = '-',markerstyle = 'o',label = None, min_seq_separation = 8,ylabel=True, xlabel = True, mode = 'substructures', norm = False, thresh = 20, markersize = 400, linewidth = 8,legend_fontsize = 26, include_label = True, color = 'C0', legend_loc = 'upper right',temp_norm = 1, legend = True, ax = None):
    """
    Plots the probabilty of being "folded" as a function of chain length. Can mean one of two things
    
    if mode =='natives':
        As a function of chain length, plots either:
            
        1. If norm = True, The probabilty that at least a fraction thresh of the native contacts
            that can form at that length are formed. thresh is 50% by default
        2. If norm = False, the probabilty that at least thresh native contacts are formed
        
        Note that thresh therefore has a different meaning depending on whether norm is True or False
    
    if mode =='substructures':
        As a function of chain length, plots probability that a certain set of desired substructures, subs, is folded
        For example, subs can be 'ade', in which case for each cain length, compute probability that at least a, d, and e are folded
        
        If you choose this mode, you have the option of entering a value for min_len...for all lengths prior to this, the substructure
        cannot form and so the probability will be plotted as 0. min_len shoudl be the length at which the first contact from substructure
        can form
        
        As for choosing what value to enter for subs: The way I've been doing it is, for each protein, I look at the full length protein and 
        I identify the cluster after the rate limiting step. Then, I look at that cluster and ask, what is the minimal set of substructures that 
        are present in all configurations assigned to that cluster? For instance for MARR, the cluster after the rate limiting step is
        [bcd, abc, bc] (at least as of this writing), and so the minimal set, and thus the value of subs, would be 'bc'
    
    Assumes last directory and fileroot correspond to native file
    
    """
    if ax==None: fig, ax = plt.subplots()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    
    #if type(subs)=='str': subs = [subs]
    
    Pfold=np.zeros(len(directories))
    for d, directory in enumerate(directories):
        if lengths[d]>=min_len:
            if mode =='substructures':
                [x, free_energies, temperatures]=joblib.load('{}/Substructure_PMF.dat'.format(directory))
                temperatures = np.unique(temperatures)
                #free_energies, x = Integrate_nonnative_contributions(free_energies, x, temperatures)        
                i=np.where(temperatures==temp)[0]                
                indices = []
            
                #FIgure out which substructures have been synthesized
                for t, top in enumerate(x):  
                    include = True
                    for b, bit in enumerate(top):
                        if bit!='1' and alphabet[b] in subs:
                            include = False
                    if include: indices.append(t)
    
            elif mode=='natives':
                [natives, free_energies, temperatures]=joblib.load('{}/Natives_PMF.dat'.format(directory))
                free_energies = free_energies[:,:,0]
                temperatures = np.unique(temperatures)
                i=np.where(temperatures==temp)[0]   
                indices = []
                
                if norm: 
                    norm = analyze_structures.count_contacts('{}/{}_0.100_Emin.pdb'.format(directories[d], fileroots[d]), d_cutoff = d_cutoff)
                    natives = natives/norm
                
                indices = np.where(natives>thresh)[0]
    
            
            indices = np.array(indices)
            Z=np.sum(np.exp(-free_energies[i,:]))
            
            if len(indices)>0:
                Pfold[d] = np.sum(np.exp(-free_energies[i, indices]))/Z
            else:
                Pfold[d]=0
        else:
            lengths[d]=np.nan
            Pfold[d]=np.nan

    
    
    if include_label:
        if temp_norm!=1 and label == None:
            label='T = {} $T_M$'.format( np.round( temp/temp_norm, 2))
        elif temp_norm ==1 and label == None:
            label = 'T = {}'.format(temp)
        ax.scatter(lengths, Pfold, s = markersize, label = label, color = color, marker = markerstyle)
    else:
        ax.scatter(lengths, Pfold, s = markersize, color = color, marker = markerstyle)
        
    ax.plot(lengths, Pfold, linewidth = linewidth, color = color, linestyle = linestyle)
    
    if min_len>0:
        xxxx=np.arange(0, min_len)  #solid line at 0 prior to min lenght
        ax.plot(xxxx, [0 for zzz in  range(len(xxxx))], color = color, linewidth = linewidth)
        
        
        #Line from 0 to first computed stability
        ax.plot([min_len - 1, np.nanmin(lengths)], [0, next(qq for qq in Pfold if not np.isnan(qq))], color = color, linewidth = linewidth )
    #for m, nsub in enumerate(nsubs):
    #    plt.scatter(lengths, Pfold[m,:], s=400)
    #    plt.plot(lengths, Pfold[m,:], linewidth=8, label = '$m$ = {}'.format(nsub))
    
    if ylabel:
        if norm and mode =='natives':
             ax.set_ylabel('Probability $\geq$ {}% native contacts'.format(np.round(100*thresh)), fontsize=fontsize, labelpad = labelpad )
        elif not norm and mode=='natives':
             ax.set_ylabel('Probability $\geq$ {} native contacts'.format(thresh), fontsize=fontsize, labelpad = labelpad )
        else:
             ax.set_ylabel('Probability folded', fontsize=fontsize, labelpad = labelpad )
    if xlabel: ax.set_xlabel('Length (AAs)', fontsize=fontsize, labelpad = labelpad)
    #plt.title('Probability that at least $m$ native substructures are folded', fontsize=32, y=1.02)
    ax.tick_params(labelsize=labelsize)
    if legend:
        leg = ax.legend(fontsize=legend_fontsize, loc = legend_loc)
        for item in leg.legendHandles:
            item.set_visible(False)



def plot_substructure_PMF(path, temps_to_plot, legend = True, legend_loc=None, upper_cutoff=25, integrate_out = [],states_to_plot = 'All' , linewidth = 1.3, alpha = 0.4,legend_fontsize = 40, ax = None, y_space=0.5, labelsize = 35,fontsize=30, label_fontsize=20, markersize = 120, temp_norm = 1):
    """
    path indicates where the PMFs that you want to plot are located...for instance 'ADK_umbrella_multistart/Substructure_PMF.dat'
    We do not plot any free eneriges with value higher than upper_cutoff
    
    y_space is how much vertical distance we require between labels
    
    Loops through every pair of points with consecutive N, (ex. N=1 and N=2), and if 
    they differ by only one substructure (ex. 01000 and 01100), draw a line between them,
    making sure lines are behind text to avoid collisions
    
    Can also apply the ad hoc function switch labels in case you want to, say, rename 110 as 101
    Then you would have optional argument switch=(1,2)
    
    You can also choose to plot only specific states, such as 'a, 'ab', 'null', etc
    BUt by default, you plot all states that have free energy less than upper cutoff 
    
    """
    PMF_info = joblib.load(path)
    x = PMF_info['tops']
    free_energies = PMF_info['top free energies'][:,:,0] #only take the first page, since second is uncertainties, which this function does not include
    temperatures = PMF_info['eq temps']
    
    

        
    temperatures=np.array(temperatures)
    
    x=np.array(x)
    
    for t,temp in enumerate(temps_to_plot):
        indices=np.where(temperatures==temp)[0][0]
        
        if type(indices)==np.int64:
            indices=np.where(temperatures==temp)[0]
            
        PMFs_to_plot=free_energies[np.array(indices), :]
        
        for f in range(len(x)):
            if np.min(PMFs_to_plot[:,f])>upper_cutoff:
                PMFs_to_plot[:,f]=np.nan
        
        
        if states_to_plot !='All':
            for s, state in enumerate(x):
                state_string = utils.barcodes_to_labels(state)
                
                if state_string not in states_to_plot:
                    PMFs_to_plot[:,s]=np.nan
        keep_indices=~np.isnan(PMFs_to_plot)
        #print(keep_indices)
        PMFs_to_plot=PMFs_to_plot[:,keep_indices[0,:]]
        x_to_plot=x[keep_indices[0,:]]
       
        N_substructures=[] #how many substructures are formed at each configuration
        
        for i, xx in enumerate(x_to_plot):
            if xx == '∅':
                N_substructures.append(0)
            else:
                N_substructures.append(len(xx))

        if ax ==None: fig, ax = plt.subplots()
        
        if temp_norm != 1:
            label='T={} $T_M$'.format(round(temp/temp_norm, 2))
        else:
            label = 'T = {}'.format(temp)
        ax.scatter(N_substructures, PMFs_to_plot,label = label, s = markersize)
        
        
        #We now annotate the plot, keeping track systematically of where we add each annotation, to avoid clashes
        N_substructures=np.array(N_substructures)
        for N in sorted(list(set(N_substructures))):   #loop throuhg all the N values that are present, in order
            annot_y_values=[] #heights for the existing annotations
            indices=np.where(N_substructures==N)[0] #Which configurations have the current value of N
            indices=indices[np.argsort(PMFs_to_plot[0,indices])]  #sort these configurations in order of increasing PMF
            
            if N<np.max(N_substructures):
                next_indices=np.where(N_substructures==N+1)[0]
            
            for i in indices:
                tentative_y=PMFs_to_plot[0,i]-0.05 #where we hope to add annotation: Right next to PMF y value
                if len(annot_y_values)>0 and tentative_y - np.max(annot_y_values)<y_space:
                        y=np.max(annot_y_values)+y_space
                else:
                    y=tentative_y
                annot_y_values.append(y)
                a=x_to_plot[i]  
                
                #if N==0: str_to_plot='$\emptyset$'
                ax.annotate(x_to_plot[i], xy=(N+0.1, y), fontsize=label_fontsize )

                
                #draw line between configurations if they differ by one substructure
                if N<np.max(N_substructures):
                    for j in next_indices:
                        b=x_to_plot[j]
                        if utils.connect_tops(a, b):
                            ax.plot([N, N+1], [PMFs_to_plot[0,i], PMFs_to_plot[0,j]], color='blue', linestyle=':', linewidth=linewidth, alpha=alpha)
                
            
        #plt.ylim((0, upper_cutoff)) 
        ax.tick_params(axis='both', labelsize=labelsize)
        ax.set_xlim(np.min(N_substructures)-0.5, np.max(N_substructures)+0.9)
        ax.set_xticks(np.arange(np.min(N_substructures), np.max(N_substructures)+1, 1))
        #plt.yticks()
        if legend:
            if legend_loc!=None:
                ax.legend(fontsize = legend_fontsize, loc = legend_loc, frameon = None)
            else:
                ax.legend( fontsize=legend_fontsize, frameon = None)
        ax.set_xlabel('Number of substructures formed', fontsize=fontsize, labelpad=15)
        ax.set_ylabel('Free energy($k_{B}T$)', fontsize=fontsize)

def plotTrajectoryTimecourses(paths, var, trajectories, legend_loc = None, title=None, ax = None, colors = ['b','g','r','c','m','y','k'], legend = True, fontsize = 30, labelsize = 30, sliding_windowlen = 1,legend_fontsize = 30,labels = [], time_offsets = None):
    """
    Loads data in All_data.dat file located in directory, and plots the desired variable for the given
    list of trajectories, for example ['0.575_0.', '0.575_10.']
    
    You input a path or a list of paths you want to plot, and it goes through them one by one and plots the same variable and 
    trajectories for each
    var is the variable you want to plot (ex. 'natives', 'energy') 
    
    By the way, if you plot multiple simulations here that have occured serially, you can set time_offset to a list other than None
    FOr instnace you ran 3 protien G simulations serially, each of which lasted 300 million MC steps, so you would enter time_offsetts = [0, 300000000, 600000000]
    
    Note: You can also do a sliding window average if you set sliding_windowlen to soethign other than 1
    """
    if type(paths)==str:
        paths = [paths]
    if ax == None:
        fig, ax = plt.subplots()
    data_to_plot = []
    times = []
    for p, path in enumerate(paths):
        data, temperatures, setpoints, log_files, tt, variables=load_data.load_log_data(path)
        i = variables.index(var)
        data_to_plot.append( data[:,:,i])        
        if time_offsets !=None:
            tt+=time_offsets[p]
        times.append(tt)
    
    data_to_plot = tuple(data_to_plot)
    data_to_plot = np.concatenate(data_to_plot, axis = 1)
    
    times = tuple(times)
    times = np.concatenate(times)
    
        
        #traj_indices = []
    for t, traj in enumerate(trajectories):
        #print(t)
        traj_index = [j for j, string in enumerate(log_files) if traj in string ][0]
        #traj_indices.append([j for j, string in enumerate(log_files) if traj in string ][0])

        roundindex = np.mod(t, 7)
        color = colors[roundindex]
        XXX = data_to_plot[traj_index, :]
        
        YYY = np.zeros(len(XXX)-sliding_windowlen)
        for TT in range(len(XXX)-sliding_windowlen):
            YYY[TT] = np.mean(XXX[TT:TT+sliding_windowlen])
        
        times = times[0:len(YYY)]

        if len(labels )==0:
            label = '{}, {}'.format(traj, var)
        else:
            label = labels[t]
        ax.plot(times,YYY , label =label, color = color )


    ax.set_xlabel('MC step', fontsize=fontsize)
    ax.ticklabel_format(axis = 'x', style = 'scientific', scilimits = (0,0))
    ax.tick_params(labelsize=labelsize)
    ax.set_ylabel(var, fontsize = fontsize)
    
    if legend:
        if legend_loc !=None:
            ax.legend(fontsize = legend_fontsize, loc = legend_loc)
        else:
            ax.legend(fontsize = legend_fontsize)
    if title!=None:
        ax.set_title(title, fontsize=fontsize)

def umbrella_PMF(x_kn, data_path, eq, k_bias, save_path,  temps_to_use='All', save=True, max_time = None):
    """
    CURRENT FUNCTION USED TO COMPUTE SUBSTRUCTURE AND 1D PMF 
    
    x_kn is some list of coordiante values along which you want to compute free energies
    We assume that x_kn is either an array where each row is a timecourse at a given condition, and the rows are in order of increasing reporter values (ex. temp or setpoint)
    OR
    that x_kn is a 1D array where a new trajectory starts every t values, and the trajectories are in order of increasing reporter values
    You can also enter None, in whcih case x_kn is simply the native contacts
    Alternatively you can enter x_kn as a string, for instnace 'rmsd', in whcih case the program extracts that variale from the data
    
    IMPORTANT: YOU NEED TO MAKE SURE THAT THE POINTS IN X_KN CORRESPOND TO SAME TIMEPOINTS AS THE POINTS IN THE DATA FILE...IF IT DOESN'T, THEN
    THE FUNCTION WILL TRY TO FIX IT BY SUBSAMPLING THE X_KN, BUT I DON'T TRUST THIS...
    
    data_path tells you where the native contacts data is located
    
    eq tells you how many steps you want to leave out initially while the simulations equilibrate
    
    k_bias is the spring constant
  
    
    save_path is where you want to save the results, for instance 'ADK_umbrella_multistart/Substructure_PMF.dat'
    
    As a model for this, see pymbar/examples/umbrella-sampling-pmf/umbrella-sampling.py
    May need to download this from github, not sure it's on home computer
    
    
    By the way, on 3/17/20, added a parameter max_time, which is the last MC step to be used in equliibrium calculations
    Typically, we use everything from eq till the end of the simulation (which is the case if max_time has its default value of None)
    But if you set some numerical value for max_time, then we'll use somethign else
    For instnace, if eq = 0 and max_time = 100000000, then we'll only use the first 100000000 MC timesteps to compute PMF
    BUt if eq is set to, say, 150000000 and max_time is kept at its default value of None, then we use everything from 150000000 and beyond
    to compute PMF
    """    
    print("# loading data...")
    
    log_file_data, temperatures, setpoints, log_files, times, variables=load_data.load_log_data(data_path)    
       
    energies_index=variables.index('energy') 
    energies=log_file_data[:,:,energies_index]
    #
    
    if 'natives' in variables:
        natives_index=variables.index('natives')
        natives=log_file_data[:, :, natives_index]
    else: #assume no umbrella biasing
        natives=np.zeros(np.shape(energies))
        k_bias=0
        
    
    if type(x_kn)==str:
        x_kn=log_file_data[:, :, variables.index(x_kn)]
    elif np.shape(x_kn)==(): 
        x_kn=natives
    
    
    setpoints=np.array(setpoints)
    temperatures=np.array(temperatures)    
    n_conditions=np.shape(natives)[0]

    x_kn=np.array(x_kn)
    
    del log_file_data
    
    if x_kn.ndim!=2:
        x_kn=np.reshape(x_kn,(n_conditions, int(len(x_kn)/n_conditions)))

    if temps_to_use!="All":
        indices_to_use=[t for t, temp in enumerate(temperatures) if temp in temps_to_use ]
        natives=natives[indices_to_use,:]
        energies=energies[indices_to_use, :]
        x_kn=x_kn[indices_to_use, :]
        n_conditions=len(indices_to_use)
        temperatures=temperatures[indices_to_use]
        setpoints=setpoints[indices_to_use]
    
    sample_frequency=int(np.shape(natives)[1]/np.shape(x_kn)[1])
    keep=np.arange(0,  np.shape(natives)[1], sample_frequency)
    natives=natives[:, keep]
    energies=energies[:, keep]
    
    times=np.array([times[t] for t in range(len(times)) if t in keep])
    
    eq_index=np.where(times==eq)[0][0]
    
    if max_time == None:
        natives=natives[:, eq_index:]
        energies=energies[:, eq_index:]
        x_kn=x_kn[:, eq_index:]
    else:
        max_index = np.where(times==max_time)[0][0]
        natives=natives[:, eq_index:max_index]
        energies=energies[:, eq_index:max_index]
        x_kn=x_kn[:, eq_index:max_index]
        
        print(times[eq_index:max_index])

   
    n_timepoints=np.shape(natives)[1]    
    

    print("# calculating potential...")
    
    
    N_k = np.array([n_timepoints for k in range(n_conditions)], np.int32)
    
    #u_kn=np.zeros()
    
    u_kln = np.zeros((n_conditions, n_conditions, n_timepoints))  
    #ukln tells you the reduced potential energy (energy/kbT + spring cost) that 
    #point n from condition k would experience if it were to occur in some (other) 
    #condition l
    for k in range(n_conditions):
        for n in range(n_timepoints):
            u_kln[k,:,n]=energies[k,n]/temperatures+k_bias * (natives[k,n] - setpoints)**2 
         #Have to add in bias by hand since the energies term from log files does not include that bias!
        
    print("# Computing normalizations...")
    mbar = pymbar.MBAR(u_kln, N_k)  #This initialization computes the log partition functions for all conditions (temperature/bias combinations)
    #dF = mbar.getFreeEnergyDifferences()[0][0,:]
    
    #In these previous steps, we compute the full trace (partition function) for all conditions (setpoint and temp combinations)...
    
    #We will now compute the free energies (partial trace over only snapshots assigned to a state) under a DIFFERENT condition (which was not represented
    #in the conditions whose normalizations we just calculated)--namely: the state in which you have no bias

    
    unique_temperatures=np.unique(temperatures)
    
    unique_x=np.unique(x_kn)

    
    print('Computing state free energies...')
    
    x_n=x_kn.flatten()
    nbins=len(unique_x)
    bin_n=np.array([np.where(unique_x==x)[0][0] for x in x_n])
    free_energies=np.zeros((len(unique_temperatures), len(unique_x)) )
    uncertainties = np.zeros((len(unique_temperatures), len(unique_x)) )
    
    for t, temp in enumerate(unique_temperatures):
        print("Computing free energy at T={}".format(temp))
        u_n=energies.flatten()/unique_temperatures[t]  #reduced potential energy at temperature we care about
        # we do NOT include bias in above formula because we want to compute the PMF specificlaly under the condition of no bias
        #f_i, df_i = mbar.computePMF(u_n, bin_n, nbins, uncertainties='from-normalization')
        f_i, df_i = mbar.computePMF(u_n, bin_n, nbins, uncertainties='from-lowest')
        #f_i, df_i = mbar.computePMF(u_n, bin_n, nbins, uncertainties='all-differences')
        
        
        
        f_i=f_i-np.min(f_i)  #set the lowest free energy to 0
        #f_i=f_i+np.log(np.sum(np.exp(-f_i)))  #normalize--this doesn't work well because you get overflow error
        free_energies[t,:]=f_i
        uncertainties[t,:] = df_i

    """
    As for the uncertainties: Since I was not computing these for many of my previous proteins, I do not want to make a new variable
    to avoid creating confusion with number of variables to be loaded by joblib
    Rather, what I will do from now on is append the uncertainties to a second page of the free_energies array
    Also, the uncertainties are scaled by sqrt(N/N_eff), where N is total nubmer of samples, and N_eff is effective number of uncorrelated samples
    """

    #First, compute statistical inefficiency using all data
    
    g = timeseries.statisticalInefficiencyMultiple(natives)
    NNN = len(x_n)
    N_eff = NNN/g
    uncertainties = uncertainties*np.sqrt(NNN/N_eff)
    
    free_energies = np.stack((free_energies, uncertainties), axis = 2)

    if save: joblib.dump([unique_x, free_energies, temperatures], save_path)
    return unique_x, unique_temperatures, free_energies


