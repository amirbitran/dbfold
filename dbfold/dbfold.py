#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:29:18 2020

@author: amirbitran
"""
import dbfold.load_data as load_data
import dbfold.analyze_structures as analyze_structures
import dbfold.compute_PMF as compute_PMF
import dbfold.folding_rates as folding_rates
import os.path
import joblib
import numpy as np
import dbfold.utils as utils
class Protein():
    def __init__(self, protein_name, native_structure):
        """
        native_structure should be the path to the equilibrated native structure for the protein
        
        eq_datapath should be the path to a .dat file that contains All_Data
        eq_scorepath should be the path to a .dat file that contains the substructure scores
        
        The data in eq_datapath and eq_scorepath must correspond to the same timespoints!
        """
        print('hi')
        self.native_structure = native_structure
        self.name = protein_name
        #self.eq_datapath = eq_datapath
        #self.eq_scorepath = eq_scorepath
    
#    def add_nonnative_path(self, name, path):
#        """
#        Add path to nonnative contact map
#        
#        name would be what you refer to these nonnatives as...it may make sense to 
#        let name refer to the cluster or topological configuration for which we want
#        to keep track of those maps
#        
#        Path refers to the path to the .dat file where the nonnative maps are contained
#        """
#        if not hasattr(self, nonnatives_dir):
#            self.nonnatives_dir = {}
#        self.nonnatives_dir.update({name: path})
#    




    def count_contacts(self, d_cutoff=6, min_seq_separation=8 ):
        """
        TODO: These parameters shoudl be made global for the protein! The d_cutoff and the min_seq_separation
        """
        return analyze_structures.count_contacts(self.native_structure, d_cutoff, min_seq_separation) 
    
    def form_clusters(self, T_A = 100):
        print("A value of T_A = {} will be used to form clusters".format( T_A))
        
        combined_trajs, clusters = folding_rates.cluster(T_A, self.temp_unfolding_info['unique labels'], 
                                                         self.temp_unfolding_info['kinetic distances'], 
                                                         self.temp_unfolding_info['state trajs'], 
                                                         self.temp_unfolding_info['key'])
        
        clusters_dic = {}
        for c, cluster in enumerate(clusters):
            clusters_dic.__setitem__(c, cluster)
        
        print("We have obtained the following clusters: \n {} \n Please make sure that unfolding rates fit well to Arrhenius equation before proceeding ".format(clusters_dic))
        self.temp_unfolding_info.update(  {"clusters": clusters,
                               "combined_trajs": combined_trajs,
                               "T_A": T_A,
                               "clusters dic": clusters_dic})
        

    def generate_subs(self, d_cutoff=6, min_seq_separation=8, contact_sep_thresh=7, min_clustersize=6, savepath = None):
        native_contacts, substructures = analyze_structures.identify_native_substructures(self.native_structure, d_cutoff, min_seq_separation, contact_sep_thresh, min_clustersize)
        native_distances = analyze_structures.find_native_contacts(self.native_structure, d_cutoff, min_seq_separation, mode = 'distances') #computes pairwise distances between pairs of residues in native structure
        self.substructures = substructures
        Nsubs = np.shape(substructures)[2]
        alphabet ='abcdefghijklmnopqrstuvwxyz' 
        folded_state = ''
        for n in range(Nsubs):
            folded_state = '{}{}'.format(folded_state, alphabet[n])
        self.folded_state = folded_state
        if savepath != None:
            joblib.dump([substructures, native_distances, d_cutoff, min_seq_separation], '{}'.format(savepath))
        
    def obtain_folding_rates(self, Arrhenius_temps = [], unfolding_transitions_of_interest = [], N_trials=1000, min_trans = 5, overwrite = False ):
        """
        Computes folding rate and error distribution (due to sampling noise in unfolding simulation) from detailed balance
        If overwrite is True, then we recompute the folding rates from scratch even if we have already computed them before, and no previous info remains
        Otherwise, we load info that we had computed previously
        
        Note that you can also use the attribute update_folding_rates to update only certain folding transitions
        while keeping others intact, for instance if you run new unfolding simulations
        """
        if os.path.exists('{}/Folding_info.dat'.format(self.eq_dir)) and not overwrite:
            print('Folding info already exists! Loading now...')
            folding_info = joblib.load('{}/folding_info.dat'.format(self.eq_dir))
            print("Folding info successfully loaded")
        else:
            data = self.temp_unfolding_info['combined_trajs']
            clusters = self.temp_unfolding_info['clusters']
            PDB_files = self.temp_unfolding_info['PDB files']
            protein_name = self.name
            
            tops = self.PMF_info['tops']
            top_free_energies = self.PMF_info['top free energies']
            eq_temperatures = self.PMF_info['eq temps']
            
            
            folding_info = folding_rates.all_folding_rates(Arrhenius_temps, unfolding_transitions_of_interest,
                                                               data, clusters, PDB_files,  protein_name, tops,
                                                               top_free_energies, eq_temperatures,  N_trials,  min_trans )
            folding_info.update({'clusters': clusters, 
                             'PDB files': PDB_files,
                             'substructures': self.substructures,
                             'combined unfolding trajs': data,
                             'clusters_dic': self.temp_unfolding_info['clusters dic']})
            joblib.dump(folding_info, '{}/Folding_info.dat'.format(self.eq_dir))
        
        self.folding_info = folding_info


    def obtain_PMFs(self, eq_step=0,  k_bias = 0.02, overwrite = False, eq_dir = '', log_data_filename = 'Equilibrium_log_data.dat', score_filename = 'Equilibrium_scores.dat'):
        """
        Computes potentials of mean force as a function of native contacts and also topological configuration
        
        eq_step is the first timestep at which we wish to compute PMF. Should be a time at which simulation already appears converged
        
        
        Set overwrite = True if we want to recompute PMFs even if they have been computed before...this may be useful, for instnace, if
        you want to recompute with new parameters such as a new eq_step value
        """
        if not hasattr(self,"eq_dir"):
            if len(eq_dir)==0:
                print('No equilibrium directory has been set! Need to set equilibrium directory with self.set_eq_path, or provide the directory using the optional argument eq_dir!')
            else:
                self.set_eq_path(eq_dir, log_data_filename = log_data_filename, score_filename = score_filename )
                
            
        
        if os.path.exists('{}/PMFs.dat'.format(self.eq_dir)) and not overwrite:
            print('PMFs already exist! Loading now...')
            PMF_info = joblib.load('{}/PMFs.dat'.format(self.eq_dir))
            print("PMFs successfully loaded")
            
            if not hasattr(self, "substructures"):
                self.substructures=PMF_info["substructures"]
                print("Substructures set for protein")
            if not hasattr(self, "f"):
                self.f = PMF_info["f"]
                print("Value of f = {} set for protein".format(self.f))
        else:
            print('Need to compute PMFs using MBAR!')
            print('Computing PMFs vs natives..')
            natives, temperatures, native_free_energies = compute_PMF.umbrella_PMF(None, self.log_datapath, eq_step, k_bias,  None, save=False)

            print('Computing PMFs vs topological configuration..')
            scores, PDB_files, substructures = load_data.load_scores(self.eq_scorepath, self.f)
            tops, temperatures, top_free_energies = compute_PMF.umbrella_PMF(scores, self.log_datapath, eq_step, k_bias, None, save=False)
            tops = utils.barcodes_to_labels(tops)
        
            print("Finished computing PMFs!")
        
            PMF_info = {"eq_step": eq_step, 
               "f": self.f, 
               "k_bias": k_bias, 
               "tops": tops, 
               "top free energies":top_free_energies, 
               "eq temps": temperatures, 
               "native values": natives,
               "native free energies": native_free_energies,
               "PDB files": PDB_files,
               "substructures": substructures}
        
            
            joblib.dump(PMF_info, '{}/PMFs.dat'.format(self.eq_dir))
        
        self.PMF_info = PMF_info

    def plot_melting_curve(self):
        compute_PMF.plot_natives_vs_temp([], self.name, PMF_infos = [self.PMF_info] )

        
    def set_f(self,f):
        """
        Sets value of f, the maximum value of the 
        average distance between residues assigned to a substructure within a snapshot
        dived by that same average distance in the native file such that a substructure is considered formed
        
        We set it here and then always use this same value for consistency
        """
        self.f=f
    
    def runHMM(self, s = 0, m = 0.1, starting_state = 'folded', unfolding_dir = '',  score_filename = 'Unfolding_scores.dat'):
        """
        See folding_rates.runHMM for description of all parameters
        setting starting_state = 'folded' will cause the starting_state variable to resort to the topological config in which all substructures are formed
        """
        if not hasattr(self, 'unfolding_path'):
            if len(unfolding_dir)==0:
                print('No unfolding data path has been set! Need to set unfolding path with self.set_unfolidng_path, or provide the directory using the optional argument unfolding_dir!')
            else:
                self.set_unfolding_path(unfolding_dir, score_filename = score_filename)
            
        if starting_state=='folded': starting_state = self.folded_state
        print("The following parameters will be used to run HMM: \n f = {} \n s = {} \n m = {} \n starting_state = {} ".format(self.f,s, m,starting_state))

        unique_labels, distance, state_trajs, key, PDB_files = folding_rates.runHMM(self.unfolding_scorepath, self.f,s, m,starting_state, plot = True)
        
        print('From this resulting distance map, you may choose to cluster by choosing an adjancy threshold T_A that results in a significant separation of exchange timescales')
        
        self.temp_unfolding_info = {"unique labels": unique_labels,
                                    "kinetic distances": distance,
                               "state trajs": state_trajs,
                               "key": key,
                               "s": s,
                               "m": m,
                               "starting state": starting_state,
                               "PDB files": PDB_files
                }
        
    def set_eq_path(self, eq_dir, log_data_filename = 'Equilibrium_log_data.dat', score_filename = 'Equilibrium_scores.dat' ):
        self.eq_dir = eq_dir
        self.log_datapath = '{}/{}'.format(eq_dir, log_data_filename)
        self.eq_scorepath = '{}/{}'.format(eq_dir, score_filename)
        print('The equilibrium log datapath has been set as {} and the scorepath as {}'.format(self.log_datapath, self.eq_scorepath))
    
    def set_unfolding_path(self, unfolding_dir, score_filename = 'Unfolding_scores.dat' ):
        self.unfolding_dir = unfolding_dir
        self.unfolding_scorepath = '{}/{}'.format(unfolding_dir, score_filename)
        print('The unfolding path has been set as {}'.format(self.unfolding_scorepath))
        
    def update_folding_rates(self, Arrhenius_temps, unfolding_transitions_of_interest,clusters_to_replace, N_trials=1000, min_trans = 5 ):
        """
        IN PROGRESS!
        Does a new folding rate calculation using the currently stored in tedmp_unfolding_info, assuming some folding rates had already been computed
        In general, the cluster numbering, and the identity of the clusters, may be different in this new unfolding info
        Thus, we must replace the old info with this new one
        This is accomplished with clusters_to_replace
        For concreteness
        Suppose that we have some initial unfolding trajectory:
            
        Cluster 0     Cluster 1           Cluster 2          Cluster 3         Cluster 4
        [111111] -> [111110,101110] -> [011110, 001110] -> [010110 ,000110] -> [000000]
        
        Now, suppose we've run a new trajectory starting from [111110] states obsreved
        in equilibrium simulation, and the clustering is as follows in that simulation:
            
        Cluster 0'    Cluster 1'         Cluster 2'  Cluster 3'
        [111110] -> [011110, 001110] -> [010110] -> [000000]
        
        
        With new rates, and moreoever for some reason we stopped seeing the fast exchange in some
        of the clusters. 
        Suppose we onluy care to replace the 1->2 and 2->3 transitions in the original clusters
        with the new 0->1 and 1->2 transitions
        Then we would set unfolding_transitions_of_interest = [(0,1),(1,2)]
        clusters_to_replace = {0:1, 1:2, 2:3} 
        That is, clusters_to_replace[i] would tell you the cluster, in the old numbering, to which
        cluster i in the new numbering corresponds
        Then the old cluster 1 would remain cluster 1 in number, but the states it contains will now be those of the
        new cluster 0'
        Likewise for 2 and 3, whoes states would be replaced by those in 1' and 2', respectively
        
        We then update the corresponding unfolding parameters
        We then compute the new free energies of ALL clusters
        
        Then from the new free energies of all clusters as well as unfolding parameters for
        all transitions, we re-compute unfolding rates and folding rates for all transitions
        
        Doing this globally accounts for the fact that cluster definitions have chnaged, and thus
        folding transitions for clusters that we may not explicitly be trying to update (for instnace,
        in this case the 1->0 transition, in the old numbering) would change, by virtue of the fact
        that the definition of cluster 1 has changed...previously it was [111110,101110],
        now it is [111110], which may have a different free energy than [111110,101110]

        """

        data_update = self.temp_unfolding_info['combined_trajs']
        clusters_update = self.temp_unfolding_info['clusters']
        PDB_files_update = self.temp_unfolding_info['PDB files']
        protein_name = self.name
        
        tops = self.PMF_info['tops']
        top_free_energies = self.PMF_info['top free energies']
        eq_temperatures = self.PMF_info['eq temps']
        
        
        activation_energies_update, prefactors_update, Bootstrapped_EAs_update, Bootstrapped_intercepts_update= folding_rates.bootstrap(Arrhenius_temps, unfolding_transitions_of_interest, 
                                                                                                                                        data_update, PDB_files_update, protein_name, 
                                                                                                                                        legend_loc='upper right', Ntrials=N_trials, 
                                                                                                                                        min_trans= min_trans)    
        
        
        
        EAs = self.folding_info['activation energies']
        prefactors = self.folding_info['prefactors']
        Bootstrapped_EAs = self.folding_info['Bootstrapped EAs']
        Bootstrapped_intercepts = self.folding_info['Bootstrapped intercepts']
        clusters = self.folding_info['clusters']
        

        
        for nn, trans in enumerate(unfolding_transitions_of_interest):
            i = trans[0]
            j = trans[1]
            
            ii = clusters_to_replace[i]
            jj = clusters_to_replace[j]
            
            EAs[ii,jj]=activation_energies_update[i,j]
            prefactors[ii,jj] = prefactors_update[i,j]
            Bootstrapped_EAs[ii,jj,0:N_trials]=Bootstrapped_EAs_update[i,j,:]
            Bootstrapped_EAs[ii,jj,N_trials:]=np.nan
            Bootstrapped_intercepts[ii,jj, 0:N_trials] = Bootstrapped_intercepts_update[i,j,:]
            Bootstrapped_intercepts[ii,jj,N_trials:]=np.nan
            
            
            clusters[ii] = clusters_update[i]
            if nn==len(trans):
                clusters[jj]=clusters_update[j]
            
        
        G=compute_PMF.cluster_free_energies( clusters, tops, top_free_energies, eq_temperatures)
        
        fold_rates, unfolding_rates, unused=folding_rates.infer_folding_rates(clusters,EAs, prefactors, G, eq_temperatures)
        
        
        self.folding_info['activation energies'] = EAs
        self.folding_info['prefactors'] = prefactors
        self.folding_info['Bootstrapped EAs'] = Bootstrapped_EAs
        self.folding_info['Bootstrapped intercepts'] = Bootstrapped_intercepts
        self.folding_info['clusters'] = clusters
        
        self.folding_info['cluster free energies'] = G
        self.folding_info['folding rates'] = fold_rates
        self.folding_info['unfolding rates'] = unfolding_rates
        
        self.folding_info['PDB files'] = self.folding_info['PDB files'] + PDB_files_update
        self.folding_info['combined unfolding trajs'] = np.concatenate((self.folding_info['combined unfolding trajs'], data_update))
        
        clusters_dic = {}
        for c, cluster in enumerate(self.folding_info['clusters']):
            clusters_dic.__setitem__(c, cluster)
        self.folding_info['clusters_dic'] = clusters_dic
        

        joblib.dump(self.folding_info, '{}/Folding_info.dat'.format(self.eq_dir))
        
        
        
        
        
        
        


        
            
            #joblib.dump(folding_info, '{}/Folding_info.dat'.format(self.eq_dir))
        
        #self.folding_info = folding_info
        



        #TODO:Check that times for both type of data match...this is a bit nontrivial to do since the 
        #substructure scores has PDB and times for all trajectories, so need to extract individual trajectories
        
#        scoretimes = load_data.get_times(PDB_files)
#        if scoretimes != times:
#            print("Warning! the times in eq_datapath and eq_scorepath do not match! \n The times in eq_datapath are {} \n while the times in eq_scorepath are {}".format(times, scoretimes))
#        else:
#            print("the times in eq_datapath and eq_scorepath do match!")
#        
#        return
        
        #TODO: A more robust system of adding new data where you check if times already exist, 
        #maybe come up with a way of combining serial simulations, etc...