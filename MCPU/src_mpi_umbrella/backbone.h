#include "define.h"
#include "misc_util.h"
Float rot_mat[3][3];            /* this is used in vector.h
                                 
                                 */
#include "vector.h"
#include "atom.h"
#include "amino.h"
#include "mpi.h"
#include "getrms.h"             /* getrms.h was moved up from includes at the bottom because it provides struct backbone */

/* global dummy variables */
FILE *DATA, *DATA_OUT, *STATUS;

/* static data structures */
Float *radii;
long int **hard_core;
struct cutoff {
  long int a;
  long int b;
} **contact_distance;
struct cell ***the_matrix;
struct amino *amino_acids; 
long int seed;
struct atom_type {
  char res_name[4];
  char atom_name[4];
  int type_num;
} atom_type_list[200];
struct angles {Float chis[100][4];};
float  prob_ang[20][100];
float  deviation_ang[20][100][4];
int    no_chi_list[20];
struct angles *rotamer_angles;
int MAX_TYPES, natom_type_list, bb_O_type, bb_N_type, bb_OXT_type;


/* contact matrices */
short nclashes, ncontacts;
short delta_nclashes, delta_contacts;
struct contact_data {
  unsigned char contacts:1;
  unsigned char delta_contacts:1;
  unsigned char clashes:1;
  unsigned char delta_clashes:1;
  unsigned char disulfide:1;
  unsigned char check_contacts:1;
  unsigned char check_clashes:1;
  unsigned char check_hbond:1;
  long int hbond;
  long int delta_hbond;
  int closehb; // if 0, r<2.5
  int delta_closehb;
} **data, **struct_data;


/* hydrogen bonding */
struct hydrogen_bond {

  float max_D2;
  float min_D2;
  float D_int;

  struct atom *donor;
  struct atom *donor2;
  struct atom *donor3;
  struct atom *donor4;
  struct atom *donor5;
  struct atom *donor6;
  struct atom *donor7;
  struct atom *acceptor;
  struct atom *acceptor2;
  struct atom *acceptor3;
  struct atom *acceptor4;
  struct atom *acceptor5;
  struct atom *acceptor6;
  struct atom *acceptor7;
  struct atom *acceptor8;

};
short **type_contacts;
short *****sct_E;
short *hbond_E;
struct hydrogen_bond **hbonds;
int n_hbond_int;
long int NO_HBOND;
struct pair *hbond_pair;
int total_hbond_pairs;
Float prev_E_hbond, E_hbond, dE_hbond;
float d_memory;

int   SEQ_DEP_HB;
float seq_hb[3][20][20];
int   helix_sheet;

char  secstr[MAXSEQUENCE];
char  path_dir[500];		/* 300>500 */

/* contact matrices for debugging */
#if DEBUG
unsigned char **debug_contacts;
unsigned char **debug_dcontacts;
short debug_ncontacts;
unsigned char **debug_clashes;
short debug_nclashes;
#endif


/* protein data */
/* These will be pointers to an array of struct atoms, representing the protein atoms */
struct atom *native, *prev_native, *orig_native, *native_Emin, *native_RMSDmin;
struct residue *native_residue;
/* VZ:
 * The following two used to be local vars in fold.h, and struct_f2 also local to move.h
 * struct_f1 holds reference structure backbone.
 * I've moved them to here because I'd like the functions in move.h to access them
 * Yes even more global namespace pollution unfortunately.
 * A struct backbone holds five struct vector's named N, CA, C, CB, O 
 * MAXSEQUENCE is macro for 350
 */
struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE];
int natoms, nresidues;    /* The number of atoms and residues in the structure */
Float Rg;   /* radius of gyration */
Float sc_rms, native_rms, new_rms, rms_Emin, rms_RMSDmin;
int backbone_accepted;
int first_atom_res;
int *helix;

/* potential globals */
Float **potential;
Float mu, hydrogen_bond;

/* globals for contacts.h */
long int distance;
long int X_int, Y_int, Z_int;
short M, N, O, P;
short *A, *B;
struct int_vector temp_xyz_int;
struct atom *temp_atom, *temp_prev_atom;
struct cell *temp_cell, *temp_cell3;
struct cell **temp_cell_array;
unsigned char *is_rotated;

/* globals for move.h */
long int move_cycle=0;
struct triplet {short a; short b; short c;};
struct triplet *residue_triplets;
short total_triplets;
struct pair {short a; short b;};
struct pair *ab, *cd;
short total_pairs, total_pairs2;
short TOTAL_SINGLE_LOOP_MOVES, TOTAL_DOUBLE_LOOP_MOVES, TOTAL_TRIPLE_LOOP_MOVES;
float jacobi_before, jacobi_after;
int   soln_no_before, soln_no_after;
int   total_ntorsions;

/* alignment data structures */
struct align_cutoff {
  Float a;
  Float b;
} **align_con_dist;
Float **align_hard_core;
char alignment_file[500];	/* 100 to 500 vzhao */
int  nalign, struct_natoms, struct_nresidues;
short struct_ncontacts, struct_nclashes; 
struct atom *struct_native; 
struct residue *struct_residue;
int *seq_to_struct, *struct_to_seq;
int *map_to_seq, *map_to_struct;
struct segment{
  int a;
  int b;
} str_segment[100], seq_segment[100];
int nseg;
char amino_letters[20] = "ARNDCQEGHILKMFPSTWYV";
/* ----------------------- */



/*Bias*/
double k_bias;
short * orig_contactstring;
short ** substructures; //AB
int number_of_contacts_setpoint;
int umbrella; /*Is umbrella simulation turned on?*/
int number_of_contacts_max;  /*Number of contacts setpoint for first node*/
int contacts_step;  /*By how much does stepoint decrease for each node? */ 
int natives;  /*Current number of native contacts for a given node*/
int min_seq_sep;  /*How far apart in sequence do two residues have to be to qualify as a contact? Only relevant for native contact calculation*/
int new_natives;
int diff_number_of_contacts_new;
int diff_number_of_contacts_current;
int n_substructures;
int *substructure_sizes; //AB
double contact_calpha_cutoff = 7.; //AB

/* program option variables */
int PRINT_PDB;
long int MC_STEPS, MC_INIT_STEPS, MC_ANNEAL_STEPS;
int SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, MC_REPLICA_STEPS, MAX_EXCHANGE;
Float MC_TEMP;
Float MC_TEMP_MIN;
Float TEMP_STEP;
int NODES_PER_TEMP;
Float STEP_SIZE;
Float ALPHA, LAMBDA, beta;
Float NON_NATIVE_REPULSION, NON_SPECIFIC_ENERGY, NATIVE_ATTRACTION;
Float weight_clash, weight_potential, weight_rms;
Float weight_hbond;
Float LATTICE_SIZE;
Float USE_GLOBAL_BB_MOVES;
Float YANG_MOVE;
Float YANG_SCALE;
int  USE_SIDECHAINS, USE_ROTAMERS, USE_ROT_PROB, NO_NEW_CLASHES;
int  SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE;
Float SIDECHAIN_NOISE;
unsigned char MATRIX_SIZE, HALF_MATRIX_SIZE;
short  READ_POTENTIAL, READ_DENSITY, USE_GO_POTENTIAL;
char PROTEIN_NAME[100];

Float rmsd_constraint;

Float dE_tor, prev_E_tor, E_tor;
Float dE_sct, prev_E_sct, E_sct;
Float dE_aro, prev_E_aro, E_aro;
Float dE_Rg, prev_E_Rg, E_Rg;
double a_PCA[MAXSEQUENCE], a_bCA[MAXSEQUENCE], phim[MAXSEQUENCE], psim[MAXSEQUENCE];
double cur_phi[MAXSEQUENCE], cur_psi[MAXSEQUENCE];
short torsion_E[MAXSEQUENCE][6][6][6][6];
short aromatic_E[9];
int   Naromatic;
int   Res_aromatic[MAXSEQUENCE];
float cluster_phi[20][NOCLUSTERS];
float cluster_psi[20][NOCLUSTERS];

/* parameter file names */
char std_file[500], std_prefix[500]; /* vzhao: added these here from backbone.c */
/* vzhao: I increased the length limits to 500 from 100/150 for these: */
char native_file[500], structure_file[500], native_directory[1000], triplet_file[500], sctorsion_file[500], sec_str_file[500], template_file[500], pdb_out_file[500], amino_data_file[500], rotamer_data_file[500], potential_file[500], atom_type_file[500], helicity_data[500], hydrogen_bonding_data[500];
char substructure_path[500]; //AB
/* file names added 12DEC02 IAH */
char min_etot_file[500], min_drms_file[500];
char aromatic_file[500];

/* MC variables */
long int mcstep, mcstep_Emin, mcstep_RMSDmin;
int sidechain_step;
int *cur_rotamers, old_rotamer;
Float low_E, prev_E, E, native_E, Emin, Emin_pot, Emin_hbond, Emin_tor, Emin_sct, Emin_aro;
float E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro;
float E_RMSDmin;
Float dE, dE_pot, prev_E_pot, E_pot;
int naccepted, n_sidechain_accepted;
int nrejected, nothers, nomove;
int sidemovedone;
struct monte_carlo {
  int backbone_move;
  int sel_res_num;
  int sel_rotamer;
  int is_phi;
  int torsion;
  Float delta_angle[4];
  Float delta_phi_angle[3];
  Float delta_psi_angle[3];
  int loop_size;
  int sel_triplet;
  int selected[3];
} mc;

struct monte_carlo_flags {
  unsigned char clashed:1;
  unsigned char init:1;
} mc_flags;


/* rotation data structures */
Float rot_mat_00, rot_mat_10, rot_mat_20, rot_mat_01, rot_mat_11, rot_mat_21, rot_mat_02, rot_mat_12, rot_mat_22;
short **rotate_natoms;       /* number of atoms affected by the rotation */
short ***rotate_atom;   /* atoms affected by the rotation */
short **loop_rotate_natoms;       /* number of atoms affected by the rotation */
short ***loop_rotate_atoms;   /* atoms affected by the rotation */
short **loop_int_rotate_natoms;       /* number of atoms affected by the rotation */
short ***loop_int_rotate_atoms;   /* atoms affected by the rotation */
char ***loop_not_rotated;
char ***not_rotated;
short **rotate_sidechain_natoms;
short ***rotate_sidechain_atom;
char ***sidechain_not_rotated;
short ***sidechain_torsion;
short *all_rotated_atoms; /* pointer to atoms rotated by the current step */
short all_rotated_natoms;
short yang_rotated_natoms;
short *yang_rotated_atoms; 
char *yang_not_rotated;

MPI_Comm mpi_world_comm;
MPI_Status mpi_status;
int nprocs, myrank, sel_num, ierr;
int *accepted_replica, *rejected_replica, *replica_index;
//int **accepted_replica, **rejected_replica, *replica_index; //AB changed to 2D array
float *Enode, *Tnode;
int  *Nnode, *Cnode;
float delta_E, delta_T, delta_N, delta_all;
float  *buf_in, *buf_out;

/*Things defined by AB related to knowledge-based moves*/
float CLUSTER_MOVE =0.33;  //originally 0.3333 
float USE_CLUSTER = 0.0;       // 0.0 to turn knowledge based moves off, was 0.1 originally (knowledge moves on)
int MAX_CLUSTERSTEP = 0;  //Added by AB...all MC steps after this will always have USE_CLUSTER set to 0


/* includes */
#include "protein_util.h"
#include "lattice_util.h"
#include "tripep_closure.h"
#include "jac_local.h"
#include "contacts.h"
#include "loop.h"
#include "pdb_util.h"
#include "rotate.h"
#include "energy.h"
#include "debug.h"
#include "init.h"
#include "align.h"
#include "hbonds.h"
#include "update.h"
#include "move.h"
#include "mu_potential.h"
#include "fold.h"
