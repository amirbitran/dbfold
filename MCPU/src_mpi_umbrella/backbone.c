#include "backbone.h"

int main(int argc, char *argv[]) {

  int i;
  




  if(argc != 2) {
    printf("ERROR! Usage: ./fold_potential config_file\nargc : %d\n", argc);
    exit(1);
  }

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  mpi_world_comm = MPI_COMM_WORLD;
  MPI_Comm_size(mpi_world_comm, &nprocs);
  MPI_Comm_rank(mpi_world_comm, &myrank);


  /* initialize data */
  SetProgramOptions(argc, argv);
  /*MC_TEMP = MC_TEMP_MIN + myrank*TEMP_STEP; */ // Temperature for this processor

  /* Note: fopen(std_file, "w") has been put into SetProgramOptions */

  seed = time(NULL); 
  seed += (int) (MC_TEMP*1000) + getpid() + myrank;
  // srand48(seed);
  set_threefry_array((unsigned long int) seed);  
  //fprintf(STATUS, "nprocs = %d, myrank = %d\n", nprocs, myrank);
  fprintf(STATUS, "---GENERAL---\n");
  fprintf(STATUS, "  seed:\t\t%ld\n\n", seed);
  fflush(STATUS);


  weight_potential = POTNTL_WEIGHT;
  weight_hbond = HBOND_WEIGHT;

  STEP_SIZE = 2*deg2rad;

  InitializeMatrix();  /*AB: Initializes a matrix whose function I'm not totally sure of at the moment*/
  InitializeProtein();
  for (i=0;i<nresidues;i++)
    total_ntorsions += native_residue[i].ntorsions;

  if (CLUSTER_MOVE && USE_CLUSTER){
    fprintf(STATUS, "WARNING! KNOWLEDGE MOVES ARE IMPLEMENTED! THESE DO NOT SATISFY DETAILED BALANCE!  \n ");
  }
  if (YANG_MOVE) {
  	fprintf(STATUS, "WARNING! YANG (LOCAL) MOVES ARE CURRENTLY ON! ACCEPTANCE CRITERION NEEDS TO BE CORRECTED, OTHERWISE SIMULATION WILL NOT SATISFY DETAILED BALANCE \n ");
  }
  
  fprintf(STATUS, "Note that if any secondary structure character is H, there's a chance that knowledge based moves (break detailed balance) are being implemented! \n");
  
  fprintf(STATUS, "Starting configuration: %s\n", native_file);
  fprintf(STATUS, "Target configuration:   %s\n", structure_file);
  fprintf(STATUS, "Template Information:   %s\n", template_file);
  fprintf(STATUS, "Potential data:         %s\n", potential_file);
  fprintf(STATUS, "Atom typing scheme:     %s\n\n", atom_type_file);
  fprintf(STATUS, "Read potential:         %d\n", READ_POTENTIAL);
  fprintf(STATUS, "Skip local contact:     %d\n", SKIP_LOCAL_CONTACT_RANGE);
  fprintf(STATUS, "Skip bb contact:        %d\n\n", SKIP_BB_CONTACT_RANGE);
  fprintf(STATUS, "Alpha:               %7.2f\n", ALPHA);
  fprintf(STATUS, "Beta:                %7.2f\n", beta);
  fprintf(STATUS, "Lambda:              %7.2f\n\n", LAMBDA);
  fprintf(STATUS, "STEP_SIZE:           %7.2f\n\n", STEP_SIZE*rad2deg);
  fprintf(STATUS, "YANG_MOVE:           %7.2f\n", YANG_MOVE);
  fprintf(STATUS, "YANG_SCALE:          %7.2f\n\n", YANG_SCALE);
  fprintf(STATUS, "CLUSTER_MOVE:        %7.2f\n", CLUSTER_MOVE);
  fprintf(STATUS, "USE_CLUSTER:         %7.2f\n", USE_CLUSTER);
  fprintf(STATUS, "NOCLUSTERS:          %4d\n", NOCLUSTERS);
  fprintf(STATUS, "CLUSTER_NOISE:       %7.2f\n\n", CLUSTER_NOISE);
  fprintf(STATUS, "beta_favor:          %7.2f\n", beta_favor);
  fprintf(STATUS, "HB_CUTOFF:           %7.2f\n", HB_CUTOFF);
  fprintf(STATUS, "HB_INNER:            %7.2f\n", HB_INNER);
  fprintf(STATUS, "HB_PENALTY:          %7.2f\n", HB_PENALTY);
  fprintf(STATUS, "SEQ_DEP_HB:          %7d\n", SEQ_DEP_HB);
  fprintf(STATUS, "RDTHREE_CON:         %7.2f\n", RDTHREE_CON);
  fprintf(STATUS, "AROMATIC_DISTANCE:   %7.2f\n\n", AROMATIC_DISTANCE);
  fprintf(STATUS, "USE_ROTAMERS         %4d\n", USE_ROTAMERS);
  fprintf(STATUS, "USE_ROT_PROB         %4d\n\n", USE_ROT_PROB);
  fprintf(STATUS, "SIDECHAIN_MOVES      %4d\n\n", SIDECHAIN_MOVES);
  fprintf(STATUS, "replica_directory:       %s\n", path_dir);
  fprintf(STATUS, "MC_REPLICA_STEPS:        %7d\n", MC_REPLICA_STEPS);
  fprintf(STATUS, "exchanges in each step:  %7d\n\n", MAX_EXCHANGE);
  fprintf(STATUS, "MC_STEPS:            %10ld\n", MC_STEPS);
  fprintf(STATUS, "MC_ANNEAL_STEPS:     %10ld\n", MC_ANNEAL_STEPS);
  fprintf(STATUS, "MC_PRINT_STEPS:      %10d\n", MC_PRINT_STEPS);
  fprintf(STATUS, "MC_PDB_PRINT_STEPS:  %10d\n\n", MC_PDB_PRINT_STEPS);
  fprintf(STATUS, "Total side torsion no.: %4d\n\n", total_ntorsions);
  Fold();

  MPI_Finalize();
  return 0;
}


