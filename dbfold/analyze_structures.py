# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:59:20 2017

@author: amirbitran

Various functions that serve to compute the contacts matrix for a series of PDB snapshots

"""

import numpy as np
from matplotlib import pyplot as plt
import sklearn
from sklearn import metrics
from dbfold.utils import loopCluster
import joblib
import copy as cp
import matplotlib.colors as cccc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
  

def compute_contacts_matrix(coords, mode='binary', thresh=7.8, min_seq_separation=8):
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
    return contacts






def compute_RG(snapshot, atom='CA'):
    """
    Radius of gyration...
    """
    coords, resis = read_PDB(snapshot, atom)
    R_cm = np.mean(coords, axis=0)
    dR = coords - R_cm
    mag_R = np.sum(dR*dR, axis=1)
    RG = np.sqrt(np.mean(mag_R))
    return RG


def count_contacts(native_file, d_cutoff, min_seq_separation):
     coords, resis=read_PDB(native_file, 'CA')
     native_contacts=compute_contacts_matrix(coords, thresh=d_cutoff, min_seq_separation=min_seq_separation)
     return int(np.sum(native_contacts))





def create_substructure_PML(PDB_path, subs_to_plot, d_cutoff, min_clustersize, contact_sep_thresh,min_seq_separation=8, substructures = [], colours = []):
    """
    Identifies substructures, then creates a pymol .pml script that draws those substructures as colored contacts directly on the pymol
    
    Ex. Create_substructure_PML('MARR_umbrella3/marr_0.100_Emin.pdb', ['a','b','c','d','e','f'], 7.8, 7, 3)
    You can also pre-enter the substructures as an optional argument
    
    Otherwise, it will compute the substrucutres using PDB_path and save the file as PDB_path but with .pml instead of .pdb
    
    You can optinally input the sequence of colors you want to use to paint the substructures (using the fancy British spelling colours)
    Otherwise, it will color things automatically using the usual default sequence
    That optional argument, if used, needs to have len equal to thhat of subs_to_plot: one color per substructure to plot
    
    
    """
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    
    if len(substructures)==0:
        native_contacts, substructures = identify_native_substructures(PDB_path, d_cutoff, min_seq_separation, contact_sep_thresh, min_clustersize, plot=False)
    
    prefix = PDB_path.split('pdb')[0]
    PML_path = '{}pml'.format(prefix)

    Nsubs = np.shape(substructures)[2]
    
    file = open(PML_path, 'w')
    file.write('bg white \n color gray \n')
    
    
    if len(colours)==0:
        colors=cm.get_cmap('jet')
    
    
    counter = 0
    for s in range(Nsubs):        
        if alphabet[s] in subs_to_plot:  
            if len(colours)==0:
                curr_color=colors((s)/(Nsubs-1 ))
            else:
                curr_color = colours[counter]
            c_hex = cccc.to_hex(curr_color)
            c_hex = '0x{}'.format(c_hex.split('#')[1])
            sub = substructures[:,:,s]
            contacts = np.where(sub)   
            substr = 'sub{}'.format(alphabet[s])
            for z in range(len(contacts[0])):
                i = contacts[0][z]+1
                j = contacts[1][z]+1
                lines = "select aa, //resi {}/CA \n select bb, //resi {}/CA \n distance {}, aa, bb  \n hide labels, {}  \n set dash_color, {}, {} \n   ".format(i, j, substr, substr, c_hex, substr)
                file.write(lines)
            
            file.write('\n set dash_gap, 0.5 \n  set dash_radius, 0.2 \n')
            counter+=1
            
    
    file.close()
            

def find_native_contacts(native_file, thresh, min_seq_separation, mode = 'binary'):
    """
    finds all native contacts from native PDB file
    """
    native_coords, resis=read_PDB(native_file, atom='CA')
    native_contacts=compute_contacts_matrix(native_coords, thresh = thresh, min_seq_separation = min_seq_separation, mode = mode )
    return native_contacts
    


def identify_native_substructures(native_file, d_cutoff, min_seq_separation, contact_sep_thresh, min_clustersize,atom='CA', labelsize = 30, fontsize = 30, max_res = None, plot=True, ax = None, native_contacts=[], verbose=False):
    """
    Identify substructures within native file contact map
    Using the following strategy
    
    We produce a contact map which is a bunch of dots
    Contacts correspond to pairs of residues that are less than d_cutoff apart
    6 Angstroms is generally a good value, but may want a larger value for helical proteins where residues interact primarily
    via sidechains, and thus the alpha carbons are further apart
    
    We only count contacts if the residues are separated by min_seq_separation along the primary sequence
    We set min_seq_separation relatively high because we don't care to include intra-helix contacts within our contact map
    
    Ce can calculate the "Manhattan" distnace between every pair of dots on that contact map 
    and build a graph of contacts in which two contacts are connected if they are less than some
    threshold distance, contact_sep_thresh, apart in the contact map 
    
    Then, we find all connected components of this graph, each of which is a substructure
    But we only keep substructures whose number of residues is at least min_clustersize, to avoid sparse contacts here and there that we dont' care about
    
    Gives you option to input native contacts a priori, but by defualt you don't do this (value is set to None)
    
    You set Max_res to something other than None if you wish to plot only up to a certain residue number (ex. to depict what substructures can be formed when the first 100 AAs are synthesized)
    """
    if len(native_contacts)==0:
        coords, resis=read_PDB(native_file, atom)
        #we get a contact map with a min seq separation larger than usual to avoid helices
   
        native_distances=compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
        
        native_contacts=np.zeros(np.shape(native_distances))
        native_contacts[np.where((native_distances<d_cutoff) & (native_distances!=0))]=1

    
    positions=np.where(native_contacts==1) #which residues participate in contacts
    positions=np.transpose(positions)

    M=metrics.pairwise.pairwise_distances(positions, metric='manhattan')  #how far is each contact from each other contact?
    
    #To find connected components, I  use my loopCluster function by feeding in the positions ofr the contacts instead of the "files",
    #as well as above matrix M as d_contacts
    
    clusters, pairs_in_substructures, mean_intercluster, mean_intracluster=loopCluster(contact_sep_thresh, positions, M, sort_orphans=False, min_clustersize=min_clustersize, verbose=verbose)


    #pairs in substructures is a list of sublists, each of which correspodns to a given substructure
    #Within a given sublist, there are a bunch of pairs which tell you which pairs of residues belong to that substructure
    
    #The above is in a messy form, so we convert it into a form that allows for numpy indexing,
    #where we have a list of sublists, each sublist contains two arrays, the first of which gives the first indices for the interacting residues
    #pairs_in_substructures=[[np.array(C)[:,0], np.array(C)[:,1]] for C in pairs_in_substructures]
    pairs_in_substructures=[(np.array(C)[:,0], np.array(C)[:,1]) for C in pairs_in_substructures]
    
    
    
    
    nsubstructures=len(pairs_in_substructures)  #we now produce a set of matrices...the ith page tells you which contacts belong to the ith substructure
    substructures=np.zeros((np.shape(native_contacts)[0], np.shape(native_contacts)[1], nsubstructures))
    for n in range(nsubstructures):
        SS=np.zeros(np.shape(native_contacts))
        SS[pairs_in_substructures[n]]=1
        substructures[:,:,n]=SS
    if plot:
        visualize_substructures(native_contacts, substructures, max_res = max_res, ax = ax, labelsize = labelsize, fontsize = fontsize)
    #print(positions)
    return native_contacts, substructures

def PDB_contacts_matrix(PDB_file, thresh=7.8, min_seq_separation=8, plot = True,mode='binary'):
    """
    Input PDB file, plots contacts matrix
    
    """

    coords, resis   = read_PDB(PDB_file, 'CA')
    M=metrics.pairwise.pairwise_distances(coords)
    M=np.tril(M, -min_seq_separation)  #-min_seq_separation enures that we do not count residues that are closer than min_seq_separation
    if mode=='binary':
        contacts=np.zeros(np.shape(M))
        contacts[np.where((M<thresh) & (M!=0))]=1
    elif mode=='distances':
        contacts=np.zeros(np.shape(M))
        contacts[M>0]=M[M>0]

    if plot:
        plt.figure()
        plt.imshow(contacts)
        plt.title(PDB_file)
    
    return contacts







def read_PDB(file, atom):
    """
    extracts coordinates for some side chain atom in some PDB file
    For instance, atom will have value 'CA' if you care about the alpha carbons
    
    TODO: Fix this so it can deal with chain labels
    Right now if the PDB has a chain label in the fifth column, this will give nonsense results
    """
    openfile=open(file)
    
    resis=[]
    coords=[]
    
    
    for line in openfile.readlines():
        #print(line)
    
        line=line.rstrip('\n')
        entries=line.split()
        if entries[0]=='ATOM':
            if entries[2]==atom and entries[4] =='A' and entries[3]!='GLY':   #So long as the current residue is not a glycine, we append the coordinate for the carbon of interest              
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
        


def score_snapshot(snapshot, substructures, atom='CA', min_seq_separation=8 ):
    """
    Assigns a set of scores for a snapshot
    the ith score tells you what is the average distnace between pairs of residues residues that participate in the ith substructure, in this snapshto
    If the score is close to the characteristic contact distnace, then the substructure should be mostly formed
    """
    coords, resis=read_PDB(snapshot, atom)
    distances=compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
    length=np.shape(distances)[0]
    len_substructures=np.shape(substructures)[0]
    if length>len_substructures: #We are applying substructures from a smaller protein to analyze a larger protein, so we only keep the part of the larger protein that is encompassed by these substructures
    	distances=distances[0:len_substructures, 0:len_substructures]
    nsubstructures=np.shape(substructures)[2]
    scores=np.zeros((nsubstructures))
    for s in range(nsubstructures): 
        sub=substructures[:,:,s]
        participation=np.multiply(distances, sub)#gives the overlap between native substrucutres and this snapshot's contacts
        scores[s]=np.mean(participation[np.nonzero(participation)])
    return scores, distances


def visualize_nonnatives(nonnatives_path, native_file, d_cutoff=6.5, cmap='Greys', Return = False, cbar = True, filter_natives = True, filter_distance = 2, vmax = 1, alpha = 1,custom_filter = None, ax=None, labelsize = 40):
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
    
    native_contacts, substructures = identify_native_substructures(native_file, d_cutoff=d_cutoff, plot=False)
    [distance_maps, PDB_files, filescores]=joblib.load(nonnatives_path)

    if np.shape(distance_maps)[2]>len(PDB_files): #there is an extra page attached to end of the distance maps that tells you mean distances between residues
        mean_distances = distance_maps[:,:,-1]
        distance_maps = distance_maps[:, :, 0:-1]
    
    
    mean_nonnatives=np.mean(distance_maps, axis=2)
    NN = np.shape(mean_nonnatives)[0]
    
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


    #if filter_natives: mean_nonnatives=np.multiply(mean_nonnatives, 1 - native_contacts)
    
    #Commented all this out September 3 2019
    #if cmap != 'Greys':
    #    for i in range(NN):
    #        for j in range(NN):
    #            if mean_nonnatives[i,j]==0:
    #                mean_nonnatives[i,j] = np.nan #makes points without any contact probability show up as white rather than peach red
    
    
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





def visualize_substructures( native_contacts, substructures, max_res = None, ax = None, labelsize = 30, fontsize = 30):
    """
    Visualizes substructures as follows
    Everything that is a native contact but not part of any substructure will have value -1 on shown image
    (Update 10/1/18, actually will only show contacts that are part of substructures)
    Meanwhile, everything that is part of substructure i (i ranges from 0 to N_substructures-1) will have value i
    Finally, all non-contacts will just be Nans and appear white
    
    Edited this on 2/4/19 so that substructures are labeled by letter rather than number
    Also reinstated the feature that contacts unassigned to substructures are visualized
    
    On 2/10/2020, Changed a bit how the script work
    Made it a bit simpler
    Also made it so unassigned contacts are now shown in gray
    """
    alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'            
    substructure_image=np.zeros(np.shape(native_contacts))    
    native_contacts=native_contacts+np.transpose(native_contacts)    
    unassigned_contacts = cp.deepcopy(native_contacts)
    for s in range(np.shape(substructures)[2]):      
        substructure_image+=(s+1)*(substructures[:,:,s]+substructures[:,:,s].transpose())
        #substructure_image+=(s+1)*substructures[:,:,s]        
        unassigned_contacts-=substructures[:,:,s]+substructures[:,:,s].transpose()

    substructure_image[substructure_image==0] = np.nan #Set background to white
    #im[im<0]=np.nan #10/1    
    #im[np.diag_indices(len(native_contacts))]=0    
    colors=cm.get_cmap('jet')    
    if ax ==None: fig, ax  = plt.subplots()
    #ax.imshow(im, cmap='jet')
    ax.imshow(substructure_image, cmap='jet')    
    ax.tick_params(labelsize=labelsize)    
    for s in range(np.shape(substructures)[2]):  
        #Let's annotate
        #y_pos=np.where(substructures[:,:,s])[0][0]-3
        y_pos=np.where(substructures[:,:,s])[0][0]+4 #2/4, decreased this from +6
        x_pos=np.where(substructures[:,:,s])[1][0]+5  #2/4, increased this from +5
        #curr_color=colors((s+1)/(np.max(substructure_image) ))
        curr_color=colors((s)/(np.nanmax(substructure_image)-1 ))
        #print(np.max(substructure_image)-1)
        ax.annotate('{}'.format(alphabet[s]), (x_pos, y_pos), fontsize=fontsize, color=curr_color)
        ax.annotate('{}'.format(alphabet[s]), (y_pos-5, x_pos-8), fontsize=fontsize, color=curr_color)
    
    nsubstructures=np.shape(substructures)[2]
    nbins=nsubstructures+1 #number of colors we are showing...add 1 to accoutn for unassigned contacts

    unassigned_contacts[unassigned_contacts==0] = np.nan
    ax.imshow(unassigned_contacts, cmap = 'gray', alpha = 0.5)    
    ax.plot(np.arange(0, len(native_contacts)), np.arange(0, len(native_contacts)), color='gray', linestyle=':'  )
    if max_res !=None:
        ax.set_xlim(( max_res, 0))
        ax.set_ylim((0, max_res))






