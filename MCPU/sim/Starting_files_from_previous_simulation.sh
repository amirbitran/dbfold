#!bin/bash

#This script takes the last files of one simulation whose output path is specified in the first input variable, 
#and saves them as input files for another simulation in a directory specified in second input variable
#The MC timestep constituting the last files of the simulation is given by the third argument

#MCPU will always recognize files of the form 0.pdb, 1.pdb...as such, this reads myrank for each output file by
#looking at the log file, and then creates a new file myrank.pdb



#Sample instance would be . ./Starting_files_from_previous_simulation.sh 1igd/MultiUmbrella 1igd/MultiUmbrella 99500000

echo "WARNING: THIS IS GONNA GET CONFUSED BY ANY PROTEIN WHOSE DIRECTORY HAS A CAPITAL T IN IT, FOR INSTANCE CTS1...FIX THE CODE TO RECOGNIZE _T_ RATHER THAN JUST T"
input_dir=$1
output_dir=$2
final_time=$3

array=( $( ls $input_dir/*.log ) )

for logfile in "${array[@]}"
do
	IFS='T' read -r protein everything_else <<< "$logfile" 
	IFS='_' read -r blank temp setpoint_log <<<"$everything_else"
	IFS='.' read -r setpoint log <<< "$setpoint_log"  
	line=$(grep 'myrank' "$logfile")
	IFS=' ' read -r blah colon myrank <<<"$line"
	IFS=',' read -r myrank blah <<<"$myrank"
	IFS='/' read -r blah blah protein <<<"$protein"
	protein=${protein%?}  
	cp $input_dir/${protein}_${temp}_$setpoint.$final_time $output_dir/$myrank.pdb
done 
