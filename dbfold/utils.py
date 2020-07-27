#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 09:12:31 2020

@author: amirbitran
"""

import matplotlib.pyplot as plt
import numpy as np
import fnmatch



def barcodes_to_labels(tops):
    if type(tops)==str:
        tops = [tops]
        return_single=True
    else:
        return_single = False
        
    new_tops = []
    alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for top in tops:
        state_string = ''
        for zz, bit in enumerate(top):
            if bit =='1':
                state_string = '{}{}'.format(state_string, alphabet[zz])
        if len(state_string)==0:
            state_string = '∅'
        new_tops.append(state_string)
    
    if return_single: new_tops = new_tops[0]
    return new_tops
        

def connect_tops(a, b):
    """
    Tells you whether topological configs a and b differ by precisely one substructure
    """
    if a=='∅':
        a = ''
    elif b=='∅':
        b = ''
    
    if np.abs(len(a) - len(b))==1:
        if set(b).issubset(set(a)) or set(a).issubset(set(b)):
            return True
        else:
            return False
    else:
        return False
    
  
def double_exp_decay(t, lamb1, lamb2, a1, a2):
    return a1*np.exp(-lamb1*t)+ a2*np.exp(-lamb2*t)
    

def exp_decay(t, lamb):
    return np.exp(-lamb*t)


    

def filter_nans(sim_labels, mode='prev'):
    """
    Wherever you see a stream of Nan's, do one of two things:
    
    1. if mode='prev', replace everything in that stream with the last non-Nan value that came before
    2. If mode='next', replace everythign in that stream with the first non-Nan value that comes next
    """
    for n,l in enumerate(sim_labels):
        if np.isnan(l):
            if mode=='prev':
                if n>0:  #you have non-nan timesteps before this to work with
                    sim_labels[n]=sim_labels[n-1]
                else: #need to switch to "next" mode
                    remaining=sim_labels[n:]
                    sim_labels[n]=remaining[np.min(np.where(~np.isnan(remaining))[0])]
                    
                
                
                
            elif mode=='next': 
                remaining=sim_labels[n:]
                if len(np.where(~np.isnan(remaining))[0])==0: #In this case, everything from here to the end is Nan, so we go with what came before
                    sim_labels[n:]=sim_labels[n-1]
                else:
                    sim_labels[n]=remaining[np.min(np.where(~np.isnan(remaining))[0])]
    return sim_labels
            

def get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times



def get_trajectory(projections, PDB_files, traj_list):
    """
    Obtains PDB files that at some temperature or set of temperatures
    traj_list can either be a numpy array of temperature values (ex. np.arange),
    a single string corresponding to am temperature  (ex. '0.600'),
    a string with wildcards (ex. '0.6**'), or a list of strings, potentially with wildcards 
    (ex. ['0.6**', '0.7**'])
    
    Also works for kinetic simulations...just enter the trajectory number
    Ex. if you want files of the form adk_0.800_3._ _ _, just enter '3' for traj_list
    """
    if type(traj_list)==str: traj_list=[traj_list]
    if type(traj_list)==np.ndarray:
        traj_list=[num2str(t) for t in traj_list]
    #traj=[f for f in range(len(PDB_files)) for traj_num in traj_list if fnmatch.fnmatch( PDB_files[f], '*_{}*'.format(traj_num)) ]  #sequence of indices corresponding to current PDB file 
    traj=[f for f in range(len(PDB_files)) for traj_num in traj_list if fnmatch.fnmatch( PDB_files[f], '*{}*'.format(traj_num)) ]  #sequence of indices corresponding to current PDB file 
    traj_files=[PDB_files[f] for f in traj]
    traj_coords=np.array([projections[f] for f in traj])
    return traj_files, traj_coords


def histogram(labels, key):
    plt.figure()
    a=plt.hist(labels, bins=len(key)-1)


def lookup_files(label,labels, PDB_files, traj='All'):
    """
    returns all files that are assigned to a given label   
    You can use traj to specify that you only care about files in a given trajectory
    The defalt is 'All' (all trajectories are considered )
    """
    files=[f for n,f in enumerate(PDB_files) if labels[n]==label ]
    if traj!='All':
         files = [f for f in files if fnmatch.fnmatch( f, '*_{}*'.format(traj))]
    return files


def lookup_labels(file, PDB_files, labels, verbose = True):
    """
    returns the label for a given file or set of files
    file is a string that is contained within a single file, or a set of files
    (ex. '720.199000', or '0.850_15')
    """
    indices=[f for f in range(len(PDB_files)) if fnmatch.fnmatch( PDB_files[f], '*{}*'.format(file)) ]
    for index in indices:
        if verbose: print("{}: {} \n".format(PDB_files[index], labels[index]))
    return labels[indices[0]]       

def loopCluster(thresh, files, d_contacts, sort_orphans=False, min_clustersize=1, verbose=True):
    """
    Only keep clusters min_clustersize elements or more 
    Note, this function will (understandably) get very upset and loop forever if any node is completely disconnected from all other nodes
    """
    if d_contacts[1,0]!=d_contacts[0,1]:
        d_contacts=d_contacts+np.transpose(d_contacts)
    contacts_red=np.zeros(np.shape(d_contacts))
    for i, row in enumerate(d_contacts):
        for j, entry in enumerate(row):
            if entry<=thresh: contacts_red[i,j]=1
    
    for j,row in enumerate(contacts_red):
        if np.max(row)==0:
            print('WARNING! Node {} is not connected to any other node...this function is about to go BIZARK!!'.format(j))
    
    #We start with the 0th structure. the goal is to find all nodes in the network that can be reached starting the 0th node entirely by traversing edges whose weight is greater than thresh
    
    #First, we find all nodes that only require one such passage
    
    
    clusters=[]
    clustered_indices=[]
    unsorted_indices=np.array(range(len(files)))
    
    
    
    while len(unsorted_indices)>0:
        curr_cluster=np.unique(np.nonzero(contacts_red[unsorted_indices[0],:])[0])
        
        
        
        new_layer=curr_cluster
        curr_len=len(curr_cluster)        
        new_len=0
        
        
        while new_len != curr_len:
            curr_len=len(curr_cluster)
            for node in new_layer:
                new_layer=np.append(new_layer, np.nonzero(contacts_red[node,:])[0])
            new_layer=np.unique(new_layer)
            
            new_layer=np.array([x for x in new_layer if not np.any(curr_cluster==x)])
            
            curr_cluster=np.append(curr_cluster, new_layer)
            new_len=len(curr_cluster)
        
        
        curr_cluster=np.sort(curr_cluster)
        
        clusters.append(list(curr_cluster))
        clustered_indices.extend(curr_cluster)
        
        
        #clusters_zipped=[]
        unsorted_indices=np.array([x for x in range(len(files)) if x not in clustered_indices ])
    
    
    singlets=[x for x in clusters if len(x)==1]   #clusters with only one file--these are set aside
    clusters=[x for x in clusters if len(x)>=min_clustersize] #only keep the "fat" clusters with more than min_clustersize files
    
    clusters_tally=np.zeros((len(files), len(clusters)))   #note the assignments in a matrix: rows are files, columns are fat clusters

    for n,x in enumerate(clusters):
        for y in x:
            clusters_tally[int(y),n]=1
    
    clustered_indices=np.where(np.sum(clusters_tally, axis=1)==1)[0]
    #unclustered_indices=np.where(np.sum(clusters_tally, axis=1)==0)[0]
    
    
    if sort_orphans:
        spots_for_singlets=np.zeros((len(files), len(clusters)))
        
        for n, singlet in enumerate(singlets):
            row=d_contacts[int(singlet[0]),:]
            indices=np.where(row==np.min(row[clustered_indices]))[0]   #yields all indices in row such that row(index)=minimum distance to current file among those files that have been classified 
            best_index=[i for i in indices if i in clustered_indices][0]    # we want
            cluster_assignment=np.where(clusters_tally[best_index,:]==1)[0][0]
        
            spots_for_singlets[int(singlet[0]), cluster_assignment]=1
        
        clusters_tally=clusters_tally+spots_for_singlets
        
    
    clusters=[]
    clustered_files=[]
    for n, col in enumerate(clusters_tally.T):
        clusters.append([i for i in range(len(col)) if col[i]!=0])
        clustered_files.append([files[i] for i in clusters[n]])
        
        if verbose: print("Cluster {}: {} \n {} files total \n".format(n, clustered_files[n], sum(col)))
        
    #go through all pairs of files and figure out if two members of pair are in the same cluster or not.
    #If they are, add the distanec between these two files to the variable intracluster distances
    #Otherwise, add this distnce to intercluster_distances
    intracluster_distances=[]
    intercluster_distances=[]
    for i in range(len(d_contacts)):
        for j in range(i+1):
            if  np.size(np.where(clusters_tally[i,:])[0])!=0 and np.size(np.where(clusters_tally[j,:])[0])!=0:
                if np.where(clusters_tally[i])[0]==np.where(clusters_tally[j])[0]:  #they are in the same cluster
                    intracluster_distances.append(d_contacts[i,j])
                else:
                    intercluster_distances.append(d_contacts[i,j])
    
    mean_intracluster=np.mean(intracluster_distances)
    mean_intercluster=np.mean(intercluster_distances)
        
    if verbose: print('Mean distance within clusters: {}'.format(mean_intracluster))
    if verbose: print('Mean distance between clusters: {}'.format(mean_intercluster))
    

    return clusters, clustered_files, mean_intercluster, mean_intracluster


def num2str(num):
    string=str(num)
    while len(string)<5:
        string='{}0'.format(string)
    return string



def temp_to_str(temp):
    """
    converts temperature from format 0.1 to format 0.100 by adding three places after decimal
    """
    temp=str(temp)
    while len(temp)<5:
        temp='{}0'.format(temp)
    return temp
        
