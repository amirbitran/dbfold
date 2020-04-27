# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:59:20 2017

@author: amirbitran

Various functions that serve to compute the contacts matrix for a series of PDB snapshots

"""

import numpy as np
import glob
#from matplotlib import pyplot as plt
import sklearn
from sklearn import metrics
 

def read_PDB(file, atom):
    """
    extracts coordinates for some side chain atom in some PDB file
    For instance, atom will have value 'CA' if you care about the alpha carbons
    """
    openfile=open(file)
    
    resis=[]
    coords=[]
    
    
    for line in openfile.readlines():
        #print(line)
    
        line=line.rstrip('\n')
        entries=line.split()
        if len(entries)>1 and entries[2]==atom and entries[4] =='A' and entries[3]!='GLY':   #So long as the current residue is not a glycine, we append the coordinate for the carbon of interest              
            resis.append(entries[3])
            coords.append([float(entries[6]), float(entries[7]), float(entries[8])])        
        elif len(entries)>1 and entries[2]==atom and entries[4] !='A' and entries[3]!='GLY':
  			#first, we debug an error that sometimes happens
            if '-' in entries[5][1:-1]:  #occasionally, when the y coordinate has a negative sign and three digits or more (ex. -100), the tab between the x and y components dissappears and entry [6] mergies into entry [5] (ex. -50-100)
                x=entries[5]
                entries[5]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[6]=x[(x[1:-1].index('-')+1):]
            if '-' in entries[6][1:-1]:  #occasionally, when the z coordinate has a negative sign and three digits or more (ex. -100), the tab between the z and y components dissappears and entry [7] mergies into entry [6] (ex. -50-100)
                x=entries[6]
                entries[6]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[7]=x[(x[1:-1].index('-')+1):] 
            resis.append(entries[3])
            coords.append([float(entries[5]), float(entries[6]), float(entries[7])])
        elif len(entries)>1 and entries[2]=='CA' and entries[4] =='A' and entries[3]=='GLY':    #But if the current residue is a glycine, we can only append the alpha carbon since there is no side chain
            resis.append(entries[3])
            coords.append([float(entries[6]), float(entries[7]), float(entries[8])])        
        elif len(entries)>1 and entries[2]=='CA' and entries[4] !='A' and entries[3]=='GLY':   
  			#first, we debug an error that sometimes happens
            if '-' in entries[5][1:-1]:  #occasionally, when the y coordinate has a negative sign and three digits or more (ex. -100), the tab between the x and y components dissappears and entry [6] mergies into entry [5] (ex. -50-100)
                x=entries[5]
                entries[5]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[6]=x[(x[1:-1].index('-')+1):]
            if '-' in entries[6][1:-1]:  #occasionally, when the z coordinate has a negative sign and three digits or more (ex. -100), the tab between the z and y components dissappears and entry [7] mergies into entry [6] (ex. -50-100)
                x=entries[6]
                entries[6]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[7]=x[(x[1:-1].index('-')+1):]            
            resis.append(entries[3])
            coords.append([float(entries[5]), float(entries[6]), float(entries[7])])
            
    coords=np.array(coords)
    return coords, resis
    
    
    
def read_PDB_old(file):
    """
    extracts coordinates for alpha carbons
    """
    openfile=open(file)
    
    resis=[]
    coords=[]
    
    
    for line in openfile.readlines():
        #print(line)
    
        line=line.rstrip('\n')
        entries=line.split()
        if len(entries)>1 and entries[2]=='CA' and entries[4] =='A':
            
    
                
            resis.append(entries[3])
            coords.append([float(entries[6]), float(entries[7]), float(entries[8])])        
        elif len(entries)>1 and entries[2]=='CA' and entries[4] !='A':
            
      #first, we debug an error that sometimes happens
            if '-' in entries[5][1:-1]:  #occasionally, when the y coordinate has a negative sign and three digits or more (ex. -100), the tab between the x and y components dissappears and entry [6] mergies into entry [5] (ex. -50-100)
                x=entries[5]
                entries[5]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[6]=x[(x[1:-1].index('-')+1):]
            if '-' in entries[6][1:-1]:  #occasionally, when the z coordinate has a negative sign and three digits or more (ex. -100), the tab between the z and y components dissappears and entry [7] mergies into entry [6] (ex. -50-100)
                x=entries[6]
                entries[6]=x[0:x[1:-1].index('-')+1]   #we ignore the first element of enries 6 in case it is a - sign--we don't care about that one
                entries[7]=x[(x[1:-1].index('-')+1):]            
            
            resis.append(entries[3])
            coords.append([float(entries[5]), float(entries[6]), float(entries[7])])
    coords=np.array(coords)
    return coords, resis   
def compute_contacts_matrix_OLD(coords, thresh=7.5):
    """
    two residues are in contact if the distance between teh alpha carbons is less than thresh 
    """
    contacts=np.zeros((len(coords),len(coords)))
    for i in range(len(coords)):
        for j in range(i-1): #we do not consider residues that are separated by 1
            dist=np.sqrt(np.sum((coords[i,:]-coords[j,:])**2))
            if dist<thresh:
                contacts[i,j]=1
    return contacts


def coarsen_contacts_matrix(contacts, spacing):
    """
    Returns a  more coarse-grained contacts matrix by looking at blocks of the original contact matrix with shape spacing x spacing
    If there is at least one contact within that block, then the block is assigned a value of 1
    The coarse-grained contact matrix that is produced tells you whether there are contacts between each pair of blocks
    
    
    Arguments: 
    contacts is the origianl contacts matrix (np array)
    spacing is the size of the blocks 
    """
    base_inds=np.arange(0, np.shape(contacts)[0], spacing)
    
    coarse_contacts=np.zeros((len(base_inds), len(base_inds)))
    for n1, b1 in enumerate(base_inds):
        for n2, b2 in enumerate(base_inds[0:n1]): #create the block of the contact matrix thatthat we are interested in
            if n1==base_inds[-1]:
                block=contacts[b1:np.shape(contacts)[0], b2:b2+spacing]
            else:
                block=contacts[b1:b1+spacing, b2:b2+spacing]
            coarse_contacts[n1, n2]=np.max(block)
    return coarse_contacts

def compute_contacts_matrix(coords, mode='binary', thresh=7.5, min_seq_separation=3, spacing=1):
    """
    much faster computation
    min_seq_separation is minimum distnce the two residues must be apart in sequence for them to be counted
    
    You can specify either of two modes:
    
    1. 'binary': Returns 1 at positions where distance is less than or equal to thresh
    2. 'distances': Returns inter-residue distance wherever this distances is less than or equal to thresh
    
    """

    M=metrics.pairwise.pairwise_distances(coords)


    M=np.tril(M, -min_seq_separation)  #-min_seq_separation enures that we do not count residues that are closer than min_seq_separation

    
    if mode=='binary':
        contacts=np.zeros(np.shape(M))
        contacts[np.where((M<thresh) & (M!=0))]=1
    elif mode=='distances':
        contacts=np.zeros(np.shape(M))
        contacts[M>0]=M[M>0]
        
    if spacing>1: contacts=coarsen_contacts_matrix(contacts, spacing)
    return contacts





def find_native_contacts(native_file):
    """
    finds all native contacts from native PDB file
    """
    native_coords, resis=read_PDB(native_file, atom='CA')
    native_contacts=compute_contacts_matrix(native_coords)
    return native_coords, resis, native_contacts
    

def Plot_contacts_map(native_file):
    """
    finds all native contacts from PDB file
    """
    native_coords, resis=read_PDB(native_file)
    native_contacts=compute_contacts_matrix(native_coords)
    plt.figure()
    plt.imshow(native_contacts)


def analyze_trajectory(simulation_trajectory, native_contacts, trunc):
    """
    Goes through a list of files simulation_trajectory and obtains the total number of native, nonnative, and overall contacts at each snapshot
    Also obtains the pairwise matrix for these different contact types, averaged over all snapshots within the trajectory

    Finally, gives you the option of computing contacts vs time while ignoring some number of residues at the end--this number is specified by trunc    

    """
    
        
        
    native_contacts=native_contacts[0:(len(native_contacts)), 0:(len(native_contacts))]
    all_contacts=np.zeros((len(native_contacts), len(native_contacts), len(simulation_trajectory)))
    natives=np.zeros((len(native_contacts), len(native_contacts), len(simulation_trajectory)))
    nonnatives=np.zeros((len(native_contacts), len(native_contacts), len(simulation_trajectory)))
    n_contacts=np.zeros(len(simulation_trajectory))
    n_natives=np.zeros(len(simulation_trajectory))
    n_nonnatives=np.zeros(len(simulation_trajectory))
    
    #we will also compute the number of native contacts vs time by neglecting the last trunc residues
    n_contacts_trunc=np.zeros(len(simulation_trajectory))
    n_natives_trunc=np.zeros(len(simulation_trajectory))
    n_nonnatives_trunc=np.zeros(len(simulation_trajectory))
    
    for n, snap in enumerate(simulation_trajectory):
        coords, resis=read_PDB(snap)
        coords=coords[0:(len(coords))]   #truncate trunc resiudes
        matrix=compute_contacts_matrix(coords)
        n_contacts[n]=np.sum(matrix)
        n_contacts_trunc[n]=np.sum(matrix[0:(len(matrix)-trunc), 0:(len(matrix)-trunc)])
        matrix2=np.multiply(matrix, native_contacts)   #matrix of native pairwise contats
        matrix3=matrix-matrix2
        n_natives[n]=np.sum(matrix2)
        n_natives_trunc[n]=np.sum(matrix2[0:(len(matrix2)-trunc), 0:(len(matrix2)-trunc)])
        n_nonnatives[n]=np.sum(matrix3)
        n_nonnatives_trunc[n]=np.sum(matrix3[0:(len(matrix3)-trunc), 0:(len(matrix3)-trunc)])
        all_contacts[:,:,n]=matrix
        natives[:,:,n]=matrix2
        nonnatives[:,:,n]=matrix3
        
    
    
    
    mean_contacts=np.mean(all_contacts[:,:,50:-1], axis=2)  #start with time 50 to give system some time to equilibrate before computing average
    mean_natives=np.mean(natives[:,:,50:-1], axis=2)
    mean_nonnatives=np.mean(nonnatives[:,:,50:-1], axis=2)
    
    
    
    return n_contacts, n_contacts_trunc, n_natives,n_natives_trunc, n_nonnatives, n_nonnatives_trunc, mean_contacts, mean_natives, mean_nonnatives



def temp_to_str(temp):
    """
    converts temperature from format 0.1 to format 0.100 by adding three places after decimal
    """
    temp=str(temp)
    while len(temp)<5:
        temp='{}0'.format(temp)
    return temp
        



def temp_sweep(directory, fileroot, native_contacts, trunc, temps, times):
    """
    Goes throuhg a bunch of temperatures specified in temp
    Applies analyze_trajectory on each one
    The PDB files are saved at timesteps specified by times
    These times are assumed to be the same at all temperatures
    """

    #temps=np.linspace(min_temp, max_temp, ntemps)
    ntemps=len(temps)
    traj_length=len(times)  #all trajectories better be this same length
    n_contacts=np.zeros((ntemps, traj_length))
    n_natives=np.zeros((ntemps, traj_length))
    n_nonnatives=np.zeros((ntemps, traj_length))
    
    n_contacts_trunc=np.zeros((ntemps, traj_length))
    n_natives_trunc=np.zeros((ntemps, traj_length))
    n_nonnatives_trunc=np.zeros((ntemps, traj_length))
    
    
    trajectories=[]
    
    
    mean_contacts=np.zeros(((len(native_contacts)), (len(native_contacts)), ntemps))
    mean_natives=np.zeros(((len(native_contacts)), (len(native_contacts)), ntemps))
    mean_nonnatives=np.zeros(((len(native_contacts)), (len(native_contacts)), ntemps))
    
    
    
    for t, temp in enumerate(temps):
        print("Analyzing temperature {}".format(temp))
        simulation_trajectory=['{}/{}_{}.{}'.format(directory, fileroot, temp_to_str(temp), time) for time in times]        
        #simulation_trajectory=glob.glob('{}/{}_{}.*0'.format(directory, fileroot, temp_to_str(temp)))
        #if len(simulation_trajectory)!=traj_length: print('WARNING: Trajectores are not all the same length!!!')
        trajectories.append(simulation_trajectory)        
        n_contacts[t,:], n_contacts_trunc[t,:], n_natives[t,:], n_natives_trunc[t,:], n_nonnatives[t,:],n_nonnatives_trunc[t,:], mean_contacts[:,:,t], mean_natives[:,:,t], mean_nonnatives[:,:,t]=analyze_trajectory(simulation_trajectory, native_contacts, trunc)
        
    return n_contacts, n_contacts_trunc, n_natives, n_natives_trunc, n_nonnatives, n_nonnatives_trunc, mean_contacts, mean_natives, mean_nonnatives, trajectories





class Simulation_analysis:
    def __init__(self, native_file, directory, fileroot, min_temp, max_temp, n_temps, trunc ):
        self.native_file=native_file
        self.directory=directory
        self.protein=fileroot
        self.ignore_final_res=trunc
        self.temps=np.linspace(min_temp, max_temp, n_temps)
        self.times=np.linspace(0,199000000, 200, dtype='int')
        print('Hi!')
    def get_native(self):
        self.native_coords, self.resis, self.native_contacts = find_native_contacts(self.native_file)
    def analyze_trajectories(self):
        self.n_contacts, self.n_contacts_trunc, self.n_natives,self.n_natives_trunc, self.n_nonnatives, self.n_nonnatives_trunc, self.mean_contacts, self.mean_natives, self.mean_nonnatives, self.trajectories=temp_sweep(self.directory, self.protein, self.native_contacts, self.ignore_final_res, self.temps, self.times)
        
        






#mean_n_natives=np.mean(n_natives)
#mean_n_nonnatives=np.mean(n_nonnatives)
#plt.set_cmap('afmhot')
#plt.matshow(mean_contacts)






