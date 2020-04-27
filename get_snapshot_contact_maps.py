import numpy as np
import joblib
import argparse
import dbfold
import dbfold.load_data
import dbfold.analyze_structures
import dbfold.utils
import glob


 
 
"""
Loops through all PDB files that you care about* and saves contacct maps
provided that the distances are at most d_cutoff
Note that this includes both natives and nonnative contacts!

*By this, we mean that the PDB files should be at some temperature of interest and potentially topological configuration 


An example instance would be:
python get_snapshot_contact_maps.py --directory=../1igd/MultiUmbrella17 --native_file=../1igd/files/1igd_0.200_20_Emin.pdb --root=1igd --temperatures='0.800, 0.825' --configs_of_interest=',b' --min_step=0 --max_step=100000000 --output_filename='testmap_null.dat'  --min_seq_separation=3 --score_filename=Equilibrium_scores.dat
This would read snapshots assigned to the null or the b state at T=0.800 or T=0.825
"""  

parser = argparse.ArgumentParser(description='Hi')
parser.add_argument("--directory", help = 'This is the path to the directory containing the PDB files we want to analyze...' )
parser.add_argument("--native_file", help = 'Path to the native file for the protein whose snapshots you want to analyze')
parser.add_argument("--root", help = 'The root for PDB files...all PDB files will be named {root}_{temperature}_{setpoint}.{MC timestep} for equilibrium simulations, and {root}_{temperature}_{traj number}.{MC timestep} for unfolding simulations')
parser.add_argument("--temperatures", default = "*.***", type = str, help = "Temperatures at which you want to run analysis, as a comma-separated string. For instnace you can type --temperatures='0.800, 0.900' or --temperatures = '0.8**' or --temperatures='*.***', the latter being the default ")
parser.add_argument("--configs_of_interest", default = "All", type = str, help = "We want to analyze snapshots assigned to these topological configurations. This should be a comma-separated string with no spaces. For instnace, you can type --configs_of_interest='bc' or --configs_of_interest=',bc' (where the former in this case is the null configuration) or --configs_of_interest='All' (the default) in which case all configs are analyzed  ")
parser.add_argument("--min_step", default = '0', help = 'minimum MC step to analyze. Defaults to 0.')
parser.add_argument("--max_step", default = 'inf', help = 'maximum MC step to analyze. Defaults to infinity.')
parser.add_argument("--score_filename", default = "Substructure_scores.dat", help = "The script expects the native substructures scores to be stored in this file within directory, which by default is assumed to be called Substructure_scores.dat  ")
parser.add_argument("--output_filename", default = "Nonnative_maps.dat", help = "A file with this name, which contains nonnative contact maps for all snapshots, will be saved in directory. By default, this will be called Nonnative_maps.dat")
parser.add_argument("--d_cutoff",  default = '6.5', help = 'distance cutoff (in angstroms) used to compute contacts. Defaults to 6.5.')
parser.add_argument("--min_seq_separation",  default = '8', help = 'min separation in sequence between two residues for them to be counted as being in contact. Defaults to 8')
parser.add_argument("--f", default=1.7, help='Maximum ratio of average distance between residues assigned to a substructure in a snapshot to that same distance in the native file such that the substructure is conisdered formed in the snapshot. Defaults to 1.7')

args = parser.parse_args()
native_file=args.native_file
directory=args.directory
fileroot = args.root
scores_of_interest= [item for item in args.configs_of_interest.split(',')]
if '' in scores_of_interest:
	zzzzz=scores_of_interest.index('')
	scores_of_interest[zzzzz]='∅'
min_step=float(args.min_step)
max_step=float(args.max_step)
score_file=args.score_filename
distance_map_name=args.output_filename
d_cutoff = float(args.d_cutoff)
min_seq_separation=int(args.min_seq_separation)
d_thresh=float(args.f)


temps = [item for item in args.temperatures.split(',')]
scores, score_files, substructures= dbfold.load_data.load_scores('{}/{}'.format(directory, score_file), d_thresh, convert_to_binary=True)
score_files = [file.split('/')[-1] for file in score_files]

PDB_files=[]

for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))
times = dbfold.utils.get_times(PDB_files)
PDB_files = [file for f, file in enumerate(PDB_files) if times[f]>=min_step and times[f]<max_step ]
PDB_files = [file.split('/')[-1] for file in PDB_files]



files_of_interest=[]

filescores=[]  #score for each file that was read

if scores_of_interest!='All':  
	#print(score_files)
	for k, f in enumerate(PDB_files):
		ind=score_files.index(f)
		curr_score = dbfold.utils.barcodes_to_labels(scores[ind])
		if curr_score in scores_of_interest:
			files_of_interest.append(f)
			filescores.append(curr_score)
else:
	files_of_interest=PDB_files
        


print('{} files total'.format(len(files_of_interest)))

distance_maps=np.zeros((np.shape(substructures)[0], np.shape(substructures)[1], len(files_of_interest)))
  
  

      
for n,snapshot in enumerate(files_of_interest):
    if np.mod(n,1000)==0:
        print(n)    
    coords, resis=dbfold.analyze_structures.read_PDB('{}/{}'.format(directory, snapshot), 'CA')
    distance_maps[:,:,n]=dbfold.analyze_structures.compute_contacts_matrix(coords, mode='distances',thresh=d_cutoff, min_seq_separation=min_seq_separation)


contacts = np.zeros(np.shape(distance_maps)) 
contacts[np.where((distance_maps<d_cutoff) & (distance_maps!=0))]=1  #page i gives the contacts matrix for the ith snapshot

#we staple a page at the end of this contacts array that gives the average distnace between residue i and j, looking only that those snapshots where i and j form a bona fide contact
page = np.multiply(distance_maps, contacts)
page[np.where(page==0)]=np.nan
page = np.nanmean(page, axis=2)


contacts= np.dstack((contacts, page ))
            
joblib.dump([contacts, files_of_interest, filescores ], open('{}/{}'.format(directory, distance_map_name), 'wb'), compress=6)
    



