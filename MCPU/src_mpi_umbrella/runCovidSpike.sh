#! /bin/bash

#SBATCH -n 450
#SBATCH -N 8
#SBATCH -J replicaCovidSpike_constraint
#SBATCH -o replicaCovidSpike.out 
#SBATCH -e replicaCovidSpike.err
#SBATCH -p shakhnovich 
#SBATCH --mem=85000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/8.2.0-fasrc01 openmpi/4.0.1-fasrc01
mpiexec -n 450 ./fold_potential_mpi ./cfg_CovidSpike_constraint > out.txt 32> err.txt 
cp -r /n/holyscratch01/shakhnovich_lab/amirbitran/sim/CovidSpike/MultiUmbrellaConstraint2 /n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/CovidSpike
