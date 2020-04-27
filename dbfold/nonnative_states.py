#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:57:46 2020

@author: amirbitran
"""
import joblib
import numpy as np
import matplotlib.colors as cccc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import dbfold.analyze_structures as analyze_structures
import dbfold.utils as utils
import copy as cp




def cluster_nonnatives(nonnatives_path, native_file, d_cutoff=6.5, min_seq_separation = 8, filter_nonnatives = True, hist = False,filter_distance=2, thresh=100,cmap='Greys', frac_to_plot = 0.1):
    """
    Loads nonnative contact maps specified in nonnative_path, eliminates native contacts as well as any register shifts by using
    native file in native_file, computes hamming distance between every pair of nonnative snapshots, then does loop clustering with thresh
    Finally, plots all clusters that represent more than fract_to_plot of snapshots
    
    filter_distance tells you how much of a register shift you consider for native contacts
    
    if hist, make a histogram for # of nonnative contacts

    """
    native_contacts = analyze_structures.find_native_contacts(native_file, d_cutoff, min_seq_separation)
    [distance_maps, PDB_files, filescores]=joblib.load(nonnatives_path)
    
    if np.shape(distance_maps)[2]>len(PDB_files): #there is an extra page attached to end of the distance maps that tells you mean distances between residues
        #mean_distances = distance_maps[:,:,-1]
        distance_maps = distance_maps[:, :, 0:-1]
    
    N=np.shape(distance_maps)[2]    
    nn=np.shape(distance_maps)[0]
        
    if filter_nonnatives:
        Filter=cp.deepcopy(native_contacts)
        for d in range(-filter_distance, filter_distance+1):  #gets rid of register-shifted native contacts
            im1_to_add=np.roll(native_contacts, d, axis=1)
            if d<0:
                im1_to_add[:, d:]=0
            else:
                im1_to_add[:, 0:d]=0
            
            im2_to_add=np.roll(native_contacts, d, axis=0)
            if d<0:
                im2_to_add[d:,:]=0
            else:
                im2_to_add[0:d, :]=0
            Filter=Filter+im1_to_add + im2_to_add
            Filter[np.where(Filter)]=1

        nonnative_contacts=np.multiply(distance_maps, 1-np.repeat(Filter[:,:,np.newaxis], N, axis=2))
    else:
        nonnative_contacts = distance_maps
        
    N_nonnatives = np.sum(np.sum(nonnative_contacts, axis=0), axis=0)
    print('Mean number of nonnative contacts is {}'.format(  np.mean( N_nonnatives  )))
    if hist:
        plt.figure()
        plt.hist(N_nonnatives, bins = 50)
        
    distances=np.zeros((N,N))        
    #can speed up by paralelizing, doing "roll"
    for i in range(N):
        v1=nonnative_contacts[:,:,0:(N-i)]
        v2=nonnative_contacts[:,:,i:N]
        d=np.sum(np.sum((v1-v2)**2, axis=1), axis=0)
        distances[ np.arange( 0,(N-i), 1), np.arange(i,N,1)]=d

    clusters, unused, mean_intercluster, mean_intracluster=utils.loopCluster(thresh, np.arange(0,N,1), distances, verbose=False)

    percentages=[100*len(clust)/N for clust in clusters]
    clusters_to_plot=[c for c, clust in enumerate(clusters) if len(clust)/N>=frac_to_plot]
    
    M=len(clusters)


    mean_maps=np.zeros((nn,nn,M))
    for m, clust in enumerate(clusters):
        mean_maps[:,:,m] = np.mean(nonnative_contacts[:,:,clust], axis=2)
        if m in clusters_to_plot:
            plt.figure()
            im=plt.imshow(mean_maps[:,:,m], cmap=cmap, vmin = 0, vmax = 1)
            plt.title('Cluster {} ({}% of frames)'.format( m,  np.round(percentages[m],2)), fontsize=40, y = 1.08)
            cbar=plt.colorbar()
            cbar.ax.tick_params(labelsize=40) 
            plt.tick_params(labelsize=40)
        
    
    clustered_PDBs=[  [PDB_files[i] for i in clust ] for c, clust in enumerate(clusters) ]
    if M==1: mean_maps=mean_maps[:,:,0]
    
    return clustered_PDBs, [p for c,p in enumerate(percentages) ], mean_maps

def visualize_nonnatives(nonnatives_path, native_file, d_cutoff=6.5, min_seq_separation = 8, cmap='Greys', Return = False, cbar = True, filter_natives = True, filter_distance = 2, vmax = 1, alpha = 1,custom_filter = None, ax=None, labelsize = 40):
    """
    Reads a file of the form Distance_maps.dat and makes a contact map of nonnative contacts with shading according to frequency with whcih 
    that contact is observed
    
    d_cutoff is distance cutoff with which you identify NATIVE structures to subtract off from the nonnatives...sholud be
    the same as whatever was used to identify the nonnatives
    
    if filter_natives, then we ignore the native contacts, as well as a border around them given by filter_distance
    You also have the option to enter a Custom filter, which is a matrix of 1's at positions where you want to filter out the contact map...by default this is off and set to none
    Note that if a custom_filter is used, you still pad that filter with a border given by filter_distance
    
    If both filter_natives is set to and and you provide a custom filter, then the two filters are used in conjunction
    
    By the way, the variable vmax says what is the strongest value in the colorbar
    By default, it's 1, but you can also set it to None in which case it becomes the maximum value in the map
    """
    
    native_contacts = analyze_structures.find_native_contacts(native_file, d_cutoff, min_seq_separation)
    [distance_maps, PDB_files, filescores]=joblib.load(nonnatives_path)

    if np.shape(distance_maps)[2]>len(PDB_files): #there is an extra page attached to end of the distance maps that tells you mean distances between residues
        #mean_distances = distance_maps[:,:,-1]
        distance_maps = distance_maps[:, :, 0:-1]
    
    
    mean_nonnatives=np.mean(distance_maps, axis=2)
    #NN = np.shape(mean_nonnatives)[0]
    
    if filter_natives or np.shape(custom_filter)!=():
        
        if filter_natives and np.shape(custom_filter)==():
            Filter=cp.deepcopy(native_contacts)

        elif filter_natives and np.shape(custom_filter)!=():
            Filter = cp.deepcopy(native_contacts) + custom_filter
            zz = np.zeros(np.shape(Filter))
            zz[np.where(Filter>0)]=1
            Filter = zz
            
        else:
            Filter = custom_filter
            
        #plt.figure()
        #plt.imshow(Filter)

        for d in range(-filter_distance, filter_distance+1):  #gets rid of register-shifted native contacts
            
            im1_to_add=np.roll(Filter, d, axis=1)
            if d<0:
                im1_to_add[:, d:]=0
            else:
                im1_to_add[:, 0:d]=0
            
            im2_to_add=np.roll(Filter, d, axis=0)
            if d<0:
                im2_to_add[d:,:]=0
            else:
                im2_to_add[0:d, :]=0
            Filter=Filter+im1_to_add + im2_to_add
            Filter[np.where(Filter)]=1
        #plt.figure()
        #plt.imshow(Filter)
        mean_nonnatives = np.multiply(mean_nonnatives, 1 - Filter)

    if vmax == None:
        vmax = np.max(mean_nonnatives)
    normalize = cccc.Normalize(vmin = 0, vmax = vmax)
        
    if ax == None:
        fig, ax = plt.subplots()
    if cmap!=None:
        #im = ax.imshow(mean_nonnatives, cmap=cmap, norm = normalize, alpha = alpha, origin = 'lower')
        im = ax.imshow(mean_nonnatives + np.transpose(mean_nonnatives), cmap=cmap, norm = normalize, alpha = alpha, origin = 'upper') #changed to this on 1/10/19
    else:
        #im = ax.imshow(mean_nonnatives, norm = normalize, alpha = alpha, origin = 'lower')
        im = ax.imshow(mean_nonnatives + np.transpose(mean_nonnatives), norm = normalize, alpha = alpha, origin = 'upper') #changed to this on 1/10/19
    #im.set_clim((0, vmax))
    
    if cbar:
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize=labelsize) 
    ax.tick_params(labelsize=labelsize)
    
    ax.plot(np.arange(0, len(mean_nonnatives)), np.arange(0, len(mean_nonnatives)), color='gray', linestyle=':'  ) #added 1/10/19
        
    if Return: return im