# MCPU
//start	

MCPU - Monte-Carlo protein simulation program. Creates and runs a grid of parallel 
Monte-Carlo simulations based on PDB files for a single-chain protein with no hydrogen 
atoms. The simulation grid has two axes--namely simulation temperature T and umbrella
setpoint s. This setpoint biases the protein to explore configurations with a native 
contacts value close to s by adding a term to the potential function given by

	E_umbrella = 1/2*k_bias*(N(x) - s)^2 (Eq. 1)

Where N(x) is the number of native contacts in the current protein configuration x,
and k_bias is the biasing spring constant is desired (e.g. in unfolding simulations). 
For each temperature, we run a number of simulations given by nodes_per_temp with 
umbrella setpoints s spaced out by a value contacts_step, with the highest value given
by s= number_of_contacts_max. For instance, if number_of_contacts_max = 100,
nodes_per_temp = 11, and contacts_step = 10, then the setpoints at each temperature will
be given by s=100, s= 90, s= 80, s = 70 ... x = 10, s = 0.All parameters are set by the 
cfg file (see step 3 in section Steps in implementation below) 


There is also the option of implementing replica exchange, where a core can exchange with 
its neighbors in the grid along either the temperature and set point directions.

NOTE: At the moment, this code works for PDB files with up to 8000 atoms and 
1000 residues. To increase this, change the values for MAX_ATPMS and MAXSEQUENCE, 
respectively, in define.h and recompile

					Directories
											
1. mcpu_prep - the directory containing the code to create input files

2. sim - the directory with files prepared for simulations of any specific protein. This is 
also where output is stored

3. src_mpi_umbrella - the directory with source code. Contains configuration files
src_mpi_umbrella/cfg 

4. config_files - the directory with parameters




					Steps in implementation


1. Create necessary input files: 
	<PDB_ID>.triple
	<PDB_ID>.sctorsion
	<PDB_ID>.sec_str

To create the first two files, run save_triple.c (in the mcpu_prep directory): 
	./save_triple <PDB_ID>
with triple.energy, sct.energy, and <PDB_ID>.fasta in the directory. This may take a few minutes to run.
NOTE: It is important that the FASTA file have 80 characters per line

Create <PDB_ID>.sec_str manually. File contains secondary structure assignment for each protein residue (see publication [1]).
first line: use input secondary structure? (9/0 = yes/no)
second line: secondary structure type (H/E/C = helix/sheet/coil)
For most applications, the first line is entirely 0's (no input secondary structure) and the second line is set to entirely C's 

Place input files, along with the pdb file, in the directory sim/DHFR/files/
(currently contains sample input files for DHFR)


2. Edit configuration options in cfg file. The most relevant options (without changing the 
potential) are:

					NATIVE PROTEIN DATA				
	NATIVE_DIRECTORY -- Contains a set of PDB files that will be used to initialize 
	simulation for each respective core. The PDB files in this directory should be named 
	0.pdb, 1.pdb, 2.pdb, etc., and the ith core will initialized with PDB file i.pdb. 
	If you wish to initialize all cores with the same input file, set this option to None, 
	and simply edit NATIVE_FILE and STRUCTURE_FILE below 
	
	NATIVE_FILE and STRUCTURE_FILE -- input PDB file for simulations (folded structure, 
	single chain, no hydrogens). Even if the NATIVE_DIRECTORY option above is set to 
	something other than None, these options should still be specified to allow for RMSD 
	computation
	
	PDB_OUT_FILE -- path to simulation output. Should be formatted as {
	path to output directory}/{protein name}, such that all simulation output will be 
	saved in {path to output directory} and the output files will all incorporate 
	{protein name} in their name.
									
					MONTE-CARLO PARAMETERS
	MC_STEPS -- length of the simulation in MC steps
	
	MC_PDB_PRINT_STEPS -- frequency of outputting coordinates to a pdb file
	
	MC_PRINT_STEPS -- frequency of outputting info to log file
	IMPORTANT!!!: The dbfold analysis package assumes that MC_PRINT_STEPS and 
	MC_PDB_PRINT_STEPS are set to the same value
	

					Replica Exchange Parameter
	MC_REPLICA_STEPS -- frequency of replica exchange. To turn off exchange, set to a 
	value greater than MC_STEPS.


					SIMULATION PARAMETERS
	MC_TEMP_MIN -- Lowest simulation temperature in grid
	
	TEMP_STEP -- Spacing between successive temperatures in grid
	
	NODES_PER_TEMP -- How many cores are assigned to each temperature. Note that the TOTAL
	number of temperatures will be given by n/nodes_per_temp, where N is the number of
	cores requested when initializing the simulation (see step 5 below)
	
	USE_CLUSTER -- Gives the frequency at which knowledge-based moves are attempted, given 
	that the function LoopBackboneMove has been called to make a move (see [4] for 
	description of knowledge-based moves). The function LoopBackboneMove is only called 
	with probability 0.33, so the overall chance of attempting a knowledge move 
	is 0.33*USE_CLUSTER. At MC steps above MAX_CLUSTERSTEP (see below), this value is 
	automatically set to 0 and knowledge-based moves are no longer used.
		NOTE: Knowledge-based moves violate detailed balance and thus, steps that 
		incorporate them should not be used to compute thermodynamic properties. 
		But these moves may nonetheless be useful to incorporate at the beginning of 
		a simulation to help the simulation find energy minima at intermediate numbers of 
		native contacts (as in [5])
		NOTE: At the moment, if the secondary structure file has any helices 
		H characters), then the value in the cfg file will be ignored and USE_CLUSTER will
		be set to 0.5 by default. This can be changed if desired by modifying the function 
		LoopBackboneMove (in Move.h) and recompiling.
		
	MAX_CLUSTERSTEP -- Largest step at which knowledge-based moves (see [4]) are to be 
	used. For all MC steps beyond this, these moves will be turned off. As above, this 
	requires the secondary structure file to have no H characters.
	



					 Umbrella parameters

	UMBRELLA--indicates whether or not umbrella sampling is to be used. 1 if so, 0 if not. 
	All subsequent parameters in this section are moot if set to 0. If set to 1, then all
	output files will be named as {PDB_OUT_FILE}_{temperature}_{setpoint}.{timestep}.
	If set to 0, then output files will be named as
	{PDB_OUT_FILE}_{temperature}_{node number}.{timestep}.
	
	
	K_BIAS -- Spring constant for umbrella biasing. See equation (Eq. 1) above. This value
	can be set to 0 as another way of turning off umbrella biasing, but in this case,
	output files will stil be named as if umbrella biasing were on
	
	NUMBER_OF_CONTACTS_MAX	-- Highest set point to be used in umbrella biasing
	
	CONTACTS_STEP -- Separation between set points. These parameters set up a simulation 
	grid with nodes_per_temp cores at each temperature, whose set points range from 
	NUMBER_OF_CONTACTS_MAX and descending in increments of CONTACTS_STEP
	
	MIN_SEQ_SEP -- Minimum separation in sequence between residues for pairs of residues 
	in the native PDB file to define a contact. Setting to values above 4 ensures that 
	short range contacts in alpha helices are not included in contacts definition.
	
	CONTACT_CALPHA_CUTOFF --Maximum distance, in angstroms, for two alpha carbons to be 
	considered in contact
	
					




	 				PARAMETER FILES
	-- direct these to the correct input files produced in step 1 . 
		- TRIPLET_ENERGY_FILE is <PDB_ID>.triple (see step 1)
		- SIDECHAIN_TORSION_FILE is <PDB_ID>.sctorsion
		- SECONDARY_STRUCTURE_FILE is <PDB_ID>.sec_str
	


Note: In publications [4] and [5], we run three kinds of simulations, each with different
parameters:

	I. Energy minimization simulations run at a temperature of MC_TEMP_MIN=0.1 with 
	NUMBER_OF_CONTACTS_MAX set to the total number of contacts in the starting structure,
	rounded to the nearest 10, and CONTACT_STEP = 10. NODES_PER_TEMP is is set
	such that the lowest setpoint is 0 (for instnace, 3 nodes if 
	NUMBER_OF_CONTACTS_MAX is 20). The energy minimization simulation is typically run for 
	20 million MC steps, as given by MC_STEPS. Upon completion, files of the form 
	{PDB_OUT_FILE}_{temperature}_{node number}_Emin.pdb will be produced. 
	The Emin file with the lowest energy (as indicated by the log files) is used as the
	equilibrated native structure in subsequent simulations.
	Alternatively, the script Simulated_annealing.sh can be used to run a simulated
	annealing protocol for minimization.
	
	
	II. Equlibrium simluations are run starting from the equilibrated native structure in 
	the above simulations with NUMBER_OF_CONTACTS_MAX set to the total number of contacts 
	in that structure, rounded to the nearest 10, and CONTACT_STEP = 10.
	NODES_PER_TEMP is is set such that the lowest setpoint is 0
	MC_TEMP_MIN is set to 0.4 and TEMP_STEP is set to 0.025, and the total number of 
	nodes requested is 25*NODES_PER_TEMP, such that the highest temperature is  
	T=1.000.
	In some cases (i.e. proteins in  publication [5]) USE_CLUSTER is set to 1 and
	MAX_CLUSTERSTEP set set to values around 200 million such that the initial portion
	of the simulation uses knowledge-based moves.
	These simulations are typically run for about 1 billion MC steps total, as set 
	by the parameter MC_STEPS. Depending on the size of the protein, this may require
	a few weeks of computation time. If it is not possible to run such long
	simulations, then shorter simulations may be run serially, and the output of 
	one simulation may be used as input for the next one using the script
	sim/Starting_files_from_previous_simulation.sh, which produces files of the form
	0.pdb, 1.pdb etc. to initialize subsequent simulations. In this case, the parameter
	NATIVE_DIRECTORY should point to the directory containing these files.
	
	III. Unfolding simulations are run with umbrella biasing turned off (UMBRELLA set to
	0) and likewise without replica exchange (MC_REPLICA_STEPS should be set to a value 
	greater than MC_STEPS). These simulations are run for 100 million MC steps with
	~50-120 nodes per temp at ~5 temperatures that are slightly above the melting 
	temperature (which can be guessed or determined by analyzing the equilibrium
	simulations above).  Typically, the lowest temperature used for unfolding simulation
	ranges from T=0.85 to T=0.95 depending on the protein, and TEMP_STEP is set to 0.025.
	These may either be initialized from the equilibrated native structure, or 
	from a set of snapshots drawn from replica simulations, in which case 
	NATIVE_DIRECTORY should point to the location of starting files with names
	0.pdb, 1.pdb, etc. 
	 
	


In all of the above simulations except III, UMBRELLA is set to 1, K_BIAS is set to 0.02,
and MC_REPLICA_STEPS is set to 500000

3. Additional parameters related to the energy function can be modified in define.h if 
desired, but this will require the program to be recompiled.
The following are weights for different energy terms (see publications [1], [3]): 

	POTNTL_WEIGHT -- contact potential

	HBOND_WEIGHT -- hydrogen bonding

	TOR_WEIGHT -- torsional energy for amino acid triplets
	
	SCT_WEIGHT -- side chain torsional energy
	
	ARO_WEIGHT -- relative orientations of aromatic residues


4. Compile (if necessary) and run
If it is necessary to re-compile the code, one can do so from src_mpi_umbrella directory
 by simply typing ./compile (which runs the commands  
 gcc -c rng.c and mpicc -O3 -o fold_potential_mpi backbone.c -lm rng.o)

To run:
mpiexec -n <# of procs> ./fold_potential_mpi <cfg filename>	



5. Data analysis
MCPU produces log files, which contain energy, number of native contacts, total number
of contacts (not typically analyzed), and RMSD as a function of MC step. Likewise,
PDB files are produced. Data from both filetypes can be analyzed using dbfold


Publications:
[1] J.S. Yang et al., Structure 15, 53 (2007)
[2] J. Tian et al., PLOS Comp. Bio., in press
[3] J. Xu, L. Huang, E. I. Shakhnovich, Proteins 79, 1704 (2011)
[4] W. Chen, J.S. Yang, E. I. Shakhnovich, Proteins 66, 682 (2007)
[5] A. Bitran, W. M. Jacobs, X. Zhai, E. I Shakhnovich, PNAS (2020)
[6]A. Bitran, W. M. Jacobs, E.I. Shakhnovich, in press


//end