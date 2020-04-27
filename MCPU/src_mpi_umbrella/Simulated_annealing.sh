#!bin/bash
#Runs MCPU iteratively at gradually decreasing temperatures starting from some input file of the from protein_0.100_Emin.pdb

#first input should be protein name (ex. FABG)
#second input should be protein root used to name output files, typically in lower case (ex. fabg)


#Assumes each temperature runs for 2 million MC steps

protein=$1
protein_root=$2

module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
temps=(0.400 0.350 0.300 0.250 0.200 0.150 0.100)

cfg_file="cfg_${protein}_sim_annealing"
echo "The cfg file is ${cfg_file}"
root='/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim'  #master directory where all PDB files are stored

cp ${root}/${protein}/files/${protein_root}_0.100_Emin.pdb ${root}/${protein}/Simulated_annealing/0.pdb

for temp in ${temps[*]}
do
	echo "Modifying cfg file for temperature ${temp}"
	awk -v var="$temp" '{OFS="\t\t\t"}{if ($1=="MC_TEMP_MIN"){$2=var}} {print}' $cfg_file>temporary && mv temporary $cfg_file  #Edit structure file to change temperature
	
	echo "Running simulation"
	mpiexec -n 1 ./fold_potential_mpi ./$cfg_file > out.txt 32> err.txt 
	
	echo "Renaming file ${protein_root}_${temp}.2000000 as 0.pdb"
	cp ${root}/${protein}/Simulated_annealing/${protein_root}_${temp}.2000000 ${root}/${protein}/Simulated_annealing/0.pdb
done
