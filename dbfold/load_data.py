#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:48:09 2020

@author: amirbitran
"""

import joblib
import natsort
import numpy as np



def combine_scores(path1, path2):
    """
    Loads unfolding score data in path1 and path2, and combines them into a single file by offsetting the trajectory nubmers in path2,
    then resaves as path1
    """
    Scores1,  PDB_files1, native_contacts, substructures= joblib.load(path1)
    Scores2,  PDB_files2, native_contacts, substructures= joblib.load(path2)
    
    
    traj_nums=[]  #x. file 'MARR/unfolding3/marr_1.000_399.99500000' will b given traj_num = 399 
    for f in PDB_files1:
        fields=f.split('_')
        traj_num=int(fields[-1].split('.')[0])
        traj_nums.append(traj_num)
        
    offset=np.max(traj_nums)+1  #to every traj number in PDB_files2, we add this value
    
    for i, F in enumerate(PDB_files2):
        fields=F.split('_')
        
        traj_num, time = fields[-1].split('.')
        traj_num=int(traj_num)+offset
        time=int(time)
        
        
        newstring=fields[0]
        for field in fields[1:-1]:
            newstring='{}_{}'.format(newstring, field)
        newstring='{}_{}.{}'.format(newstring, traj_num, time)
        PDB_files2[i]=newstring
    
    
    Scores=np.vstack((Scores1, Scores2))
    PDB_files=PDB_files1+PDB_files2
    #N_nonnative_substructures= np.hstack((N_nonnative_substructures1, N_nonnative_substructures2))
    
    joblib.dump([Scores, PDB_files, native_contacts, substructures], path1, compress=3 )
    


def get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(int(entries[-1]))
    return times


def load_log_data(path, reverse_order=False):
    """
    loads data and sorts files/data in order
    Use this function for data that was read ever since I started doing umbrella sampling

    reverse_order=True means you want to sort in order of decreasing setpoints
    """
    [data, variables, log_files, temperatures, setpoints, times]=joblib.load(path)
   
    
    
    #Figure out how to re-sort data to be in order of increasing temperatures, and increasing setpoints within each temperature
    order=[]
    unique_temperatures=np.array(list(set(temperatures))) 
    unique_temperatures=unique_temperatures[np.argsort(unique_temperatures)] #sort the temperatures first
    for temp in unique_temperatures:
        indices_at_temp=np.array([t for t in range(len(temperatures)) if temperatures[t]==temp ])
        
        if reverse_order:
            order_at_this_temp=indices_at_temp[np.argsort(-1*np.array(setpoints)[indices_at_temp])]
        else:
            order_at_this_temp=indices_at_temp[np.argsort(np.array(setpoints)[indices_at_temp])]
        for o in order_at_this_temp:
            order.append(o)  
        
    #order=np.argsort(reporters)
    data=data[order, :, :]
    log_files=[log_files[x] for x in order]
    temperatures=np.array(temperatures)[order]
    setpoints=np.array(setpoints)[order]
    return data, temperatures, setpoints, log_files, times, variables


def load_scores(score_path, thresh, convert_to_binary=True):
    """
    Thresh tells you how far residues in a substructure have to be relative to native distance, on average, for substructure to be considered formed
    """
    
    
    Scores, PDB_files, native_distances, substructures= joblib.load(score_path)
    
    #native contacts is actually DISTANCES
    
    a, Scores=sort_data(PDB_files, Scores)

    PDB_files=a

    
    mean_substructure_distances=[]
    for i in range(np.shape(substructures)[2]):
        x=np.multiply(substructures[:,:,i], native_distances)
        mean_substructure_distances.append(np.nanmean(x[np.where(x)]))
    mean_substructure_distances=np.array(mean_substructure_distances)
    if convert_to_binary:
        Scores=substructure_scores_to_barcodes(Scores/mean_substructure_distances, thresh)


    return Scores, PDB_files, substructures





def sort_data(PDB_files, data):
    """
    Sorts data by name according to PDB files
    """
    sort=natsort.natsorted(zip(PDB_files, data))
    PDB_files=[tuple[0] for tuple in sort]
    if type(data)==list:
        data=[tuple[1] for tuple in sort]
    else:
        data=np.vstack([tuple[1] for tuple in sort])
    return PDB_files, data





def substructure_scores_to_barcodes(data,  thresh ):
    """
    Input an array of substructure scores, data
    """
    barcodes=[]
    for i in range(len(data)):
        string=''
        for j in range(len(data[i,:])):
            if data[i,j]<=thresh:
                string='{}1'.format(string)
            else:
                string='{}0'.format(string)
        barcodes.append(string)
    return barcodes