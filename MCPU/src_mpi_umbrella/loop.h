#define MAX_RES 500
#define MSAR    15 //Max Sidechain Atoms per Residue

#include "in_out.h"

int res_atomno[MAX_RES];
int is_template[MAX_RES];

void integloop(float step_size, int *n_soln);
void get_orgco(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3], double rorg_s[5][MSAR][3], int ns[5], int n0);
void get_orgat(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3], double rorg_s[5][MSAR][3], int ns[5], int r, int i);
void get_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int n0);
void get_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int r, int i);
void put_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int n0);
void put_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int r, int i);
void yangloop(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], double b_len[6], double b_ang[7], int n0, float dih_ch, char res_name[5][4], int *n_soln);
void yang_rotate(int n0, int is_phi);
void get_template();
void check_bb();
void check_phipsi();
void initialize_torsion();
void initialize_sct();
void initialize_aromatic();
int  will2edo(int );
void read_cluster();
void initialize_secstr();

//========================================================================================================
void integloop(float step_size, int *n_soln) //Does a local move...based on the Dill paper
 {
  double r_n[5][3], r_a[5][3], r_c[5][3], r_o[5][3], r_s[5][MSAR][3];
  double rorg_n[5][3], rorg_a[5][3], rorg_c[5][3], rorg_o[5][3], rorg_s[5][MSAR][3];
  char out_pdb[100];
  char res_name[5][4];
  int  write_out_pdb = 0;

  double b_len[6], b_ang[7];
  double t_len[6], t_ang[7];
  double o_len[6], o_ang[7];
  double small_value = 0.001;
  int ns[5];
  int sel_res, n0;
  static int step;
//  double temp1[3], temp2[3], temp3[3], ang_nac, ang_side, len_on, len_ns;
//  int s2;
  static int success = 0;
//--------------------------------------------------------------------------------------------------------
  *n_soln = 0; //Number of solutions to 16th degree polynomial
  step++;
	  
  mc.loop_size = 1;
  all_rotated_natoms = 0;
  total_pairs=total_pairs2 = 0;
  total_hbond_pairs = 0;

  res_atomno[nresidues] = natoms;
    
  sel_res = (int) (threefryrand()*nresidues); //Select a residue to initialize local move
//  sel_res = 20;
  mc.is_phi = (int) (threefryrand()*2);  //With 50% probability, rotate about phi angle, otherwise psi angle
//  mc.is_phi = 1;
//fprintf(STATUS, "%f\n", deg2rad*2.);
  step_size = YANG_SCALE*2.*deg2rad*GaussianNum(); //Magnitude  of driver rotation, in radians
//  step_size = threefryrand()*1.0*pi - 0.5*pi;
//  step_size = 0.;

  if (mc.is_phi)
// residues are to be formed along (-) direction.
   {
    n0 = sel_res-3;
    if (sel_res < 4) //Don't make a move, since you don't have enough residues before to work with 
     {
//      fprintf(STATUS, "Too left! n0: %d\n", n0);
      return;
     }
    if ((native_residue[sel_res].amino_num==14) || (native_residue[sel_res-1].amino_num==14) ||
        (native_residue[sel_res-2].amino_num==14) || (native_residue[sel_res-3].amino_num==14))
     {
//      fprintf(STATUS, "Can't rotate phi angle of proline! n0: %d\n", n0);
      return;
     }
    if ((is_template[sel_res]==1) || (is_template[sel_res-1]==1) ||
        (is_template[sel_res-2]==1) || (is_template[sel_res-3]==1))
      return;
   }
  else
// residues are to be formed along (+) direction.
   {
    n0 = sel_res+1;
    if (sel_res > nresidues-5)
     {
//      fprintf(STATUS, "Too right! n0: %d\n", n0);
      return;
     }
    if ((native_residue[sel_res+1].amino_num==14) || (native_residue[sel_res+2].amino_num==14) ||
        (native_residue[sel_res+3].amino_num==14))
     {
//      fprintf(STATUS, "Can't rotate phi angle of proline! n0: %d\n", n0);
      return;
     }
    if ((is_template[sel_res]==1) || (is_template[sel_res+1]==1) ||
        (is_template[sel_res+2]==1) || (is_template[sel_res+3]==1))
      return;
   }
//  fprintf(STATUS, "Local move at %d\n", n0);
  strcpy(res_name[0], native_residue[n0-1].res);
  strcpy(res_name[1], native_residue[n0].res);
  strcpy(res_name[2], native_residue[n0+1].res);
  strcpy(res_name[3], native_residue[n0+2].res);
  strcpy(res_name[4], native_residue[n0+3].res);

  get_orgco(rorg_n, rorg_a, rorg_c, rorg_o, rorg_s, ns, n0);
  get_coord(r_n, r_a, r_c, r_o, r_s, ns, n0);

  c_bnd_len(rorg_a[1], rorg_c[1], &o_len[0]);
  c_bnd_len(rorg_c[1], rorg_n[2], &o_len[1]);
  c_bnd_len(rorg_n[2], rorg_a[2], &o_len[2]);
  c_bnd_len(rorg_a[2], rorg_c[2], &o_len[3]);
  c_bnd_len(rorg_c[2], rorg_n[3], &o_len[4]);
  c_bnd_len(rorg_n[3], rorg_a[3], &o_len[5]);
  c_bnd_ang(rorg_n[1], rorg_a[1], rorg_c[1], &o_ang[0]);
  c_bnd_ang(rorg_a[1], rorg_c[1], rorg_n[2], &o_ang[1]);
  c_bnd_ang(rorg_c[1], rorg_n[2], rorg_a[2], &o_ang[2]);
  c_bnd_ang(rorg_n[2], rorg_a[2], rorg_c[2], &o_ang[3]);
  c_bnd_ang(rorg_a[2], rorg_c[2], rorg_n[3], &o_ang[4]);
  c_bnd_ang(rorg_c[2], rorg_n[3], rorg_a[3], &o_ang[5]);
  c_bnd_ang(rorg_n[3], rorg_a[3], rorg_c[3], &o_ang[6]);
  
  c_bnd_len(r_a[1], r_c[1], &b_len[0]);
  c_bnd_len(r_c[1], r_n[2], &b_len[1]);
  c_bnd_len(r_n[2], r_a[2], &b_len[2]);
  c_bnd_len(r_a[2], r_c[2], &b_len[3]);
  c_bnd_len(r_c[2], r_n[3], &b_len[4]);
  c_bnd_len(r_n[3], r_a[3], &b_len[5]);
  c_bnd_ang(r_n[1], r_a[1], r_c[1], &b_ang[0]);
  c_bnd_ang(r_a[1], r_c[1], r_n[2], &b_ang[1]);
  c_bnd_ang(r_c[1], r_n[2], r_a[2], &b_ang[2]);
  c_bnd_ang(r_n[2], r_a[2], r_c[2], &b_ang[3]);
  c_bnd_ang(r_a[2], r_c[2], r_n[3], &b_ang[4]);
  c_bnd_ang(r_c[2], r_n[3], r_a[3], &b_ang[5]);
  c_bnd_ang(r_n[3], r_a[3], r_c[3], &b_ang[6]);
  
//Now we call yangloop, which actually does the rotation  
  
  yangloop(r_n, r_a, r_c, r_o, r_s, ns, b_len, b_ang, n0, step_size, res_name, n_soln);
  c_bnd_len(r_a[1], r_c[1], &t_len[0]);
  c_bnd_len(r_c[1], r_n[2], &t_len[1]);
  c_bnd_len(r_n[2], r_a[2], &t_len[2]);
  c_bnd_len(r_a[2], r_c[2], &t_len[3]);
  c_bnd_len(r_c[2], r_n[3], &t_len[4]);
  c_bnd_len(r_n[3], r_a[3], &t_len[5]);
  c_bnd_ang(r_n[1], r_a[1], r_c[1], &t_ang[0]);
  c_bnd_ang(r_a[1], r_c[1], r_n[2], &t_ang[1]);
  c_bnd_ang(r_c[1], r_n[2], r_a[2], &t_ang[2]);
  c_bnd_ang(r_n[2], r_a[2], r_c[2], &t_ang[3]);
  c_bnd_ang(r_a[2], r_c[2], r_n[3], &t_ang[4]);
  c_bnd_ang(r_c[2], r_n[3], r_a[3], &t_ang[5]);
  c_bnd_ang(r_n[3], r_a[3], r_c[3], &t_ang[6]);
  if(*n_soln > 0)
   {
    if(fabs(t_len[0]-o_len[0])>small_value)
     {
//      fprintf(STATUS, "distance 0: %6.3lf %6.3lf\n", t_len[0], o_len[0]);
      return;
     }
    if(fabs(t_len[1]-o_len[1])>small_value)
     {
//      fprintf(STATUS, "distance 1: %6.3lf %6.3lf\n", t_len[1], o_len[1]);
      return;
     }
    if(fabs(t_len[2]-o_len[2])>small_value)
     {
//      fprintf(STATUS, "distance 2: %6.3lf %6.3lf\n", t_len[2], o_len[2]);
      return;
     }
    if(fabs(t_len[3]-o_len[3])>small_value)
     {
//      fprintf(STATUS, "distance 3: %6.3lf %6.3lf\n", t_len[3], o_len[3]);
      return;
     }
    if(fabs(t_len[4]-o_len[4])>small_value)
     {
//      fprintf(STATUS, "distance 4: %6.3lf %6.3lf\n", t_len[4], o_len[4]);
      return;
     }
    if(fabs(t_len[5]-o_len[5])>small_value)
     {
//      fprintf(STATUS, "distance 5: %6.3lf %6.3lf\n", t_len[5], o_len[5]);
      return;
     }

    if(fabs(t_ang[0]-o_ang[0])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 0:  %6.3lf %6.3lf\n", t_ang[0]*rad2deg, o_ang[0]*rad2deg);
      return;
     }
    if(fabs(t_ang[1]-o_ang[1])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 1:  %6.3lf %6.3lf\n", t_ang[1]*rad2deg, o_ang[1]*rad2deg);
      return;
     }
    if(fabs(t_ang[2]-o_ang[2])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 2:  %6.3lf %6.3lf\n", t_ang[2]*rad2deg, o_ang[2]*rad2deg);
      return;
     }
    if(fabs(t_ang[3]-o_ang[3])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 3:  %6.3lf %6.3lf\n", t_ang[3]*rad2deg, o_ang[3]*rad2deg);
      return;
     }
    if(fabs(t_ang[4]-o_ang[4])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 4:  %6.3lf %6.3lf\n", t_ang[4]*rad2deg, o_ang[4]*rad2deg);
      return;
     }
    if(fabs(t_ang[5]-o_ang[5])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 5:  %6.3lf %6.3lf\n", t_ang[5]*rad2deg, o_ang[5]*rad2deg);
      return;
     }
    if(fabs(t_ang[6]-o_ang[6])>small_value*deg2rad*100.)
     {
//      fprintf(STATUS, "angle 6:  %6.3lf %6.3lf\n", t_ang[6]*rad2deg, o_ang[6]*rad2deg);
      return;
     }
//    fprintf(STATUS, "passed\n");
   }
  
  if(*n_soln > 0)
   {
    put_coord(r_n, r_a, r_c, r_o, r_s, ns, n0);
    yang_rotate(n0, mc.is_phi);
    all_rotated_natoms = yang_rotated_natoms;
    all_rotated_atoms = yang_rotated_atoms;
    UpdateLattice(yang_rotated_natoms, yang_rotated_atoms);
    NewDeltaContacts(yang_rotated_natoms, yang_rotated_atoms, yang_not_rotated);
    success++;
   }

/*  s2 = 51;
  temp1[0] = (double)native[res_atomno[s2]+0].xyz.x;
  temp1[1] = (double)native[res_atomno[s2]+0].xyz.y;
  temp1[2] = (double)native[res_atomno[s2]+0].xyz.z;
  temp2[0] = (double)native[res_atomno[s2]+1].xyz.x;
  temp2[1] = (double)native[res_atomno[s2]+1].xyz.y;
  temp2[2] = (double)native[res_atomno[s2]+1].xyz.z;
  temp3[0] = (double)native[res_atomno[s2]+2].xyz.x;
  temp3[1] = (double)native[res_atomno[s2]+2].xyz.y;
  temp3[2] = (double)native[res_atomno[s2]+2].xyz.z;
  c_bnd_ang(temp1, temp2, temp3, &ang_nac);
  temp1[0] = (double)native[res_atomno[s2]+2].xyz.x;
  temp1[1] = (double)native[res_atomno[s2]+2].xyz.y;
  temp1[2] = (double)native[res_atomno[s2]+2].xyz.z;
  temp2[0] = (double)native[res_atomno[s2]+4].xyz.x;
  temp2[1] = (double)native[res_atomno[s2]+4].xyz.y;
  temp2[2] = (double)native[res_atomno[s2]+4].xyz.z;
  temp3[0] = (double)native[res_atomno[s2]+5].xyz.x;
  temp3[1] = (double)native[res_atomno[s2]+5].xyz.y;
  temp3[2] = (double)native[res_atomno[s2]+5].xyz.z;
  c_bnd_ang(temp1, temp2, temp3, &ang_side);
  temp1[0] = (double)native[res_atomno[s2+1]+0].xyz.x;
  temp1[1] = (double)native[res_atomno[s2+1]+0].xyz.y;
  temp1[2] = (double)native[res_atomno[s2+1]+0].xyz.z;
  temp2[0] = (double)native[res_atomno[s2]+3].xyz.x;
  temp2[1] = (double)native[res_atomno[s2]+3].xyz.y;
  temp2[2] = (double)native[res_atomno[s2]+3].xyz.z;
  c_bnd_len(temp1, temp2, &len_on);
  temp1[0] = (double)native[res_atomno[s2]+0].xyz.x;
  temp1[1] = (double)native[res_atomno[s2]+0].xyz.y;
  temp1[2] = (double)native[res_atomno[s2]+0].xyz.z;
  temp2[0] = (double)native[res_atomno[s2]+4].xyz.x;
  temp2[1] = (double)native[res_atomno[s2]+4].xyz.y;
  temp2[2] = (double)native[res_atomno[s2]+4].xyz.z;
  c_bnd_len(temp1, temp2, &len_ns);*/
//  if (step%10000 == 0)
//    fprintf(STATUS, "step: %9d res. no: %3d res. name: %s %9.3f\n", step, s2, native_residue[s2].res, ang_nac*rad2deg);

//  if (step%10000 == 0)
//    fprintf(STATUS, "step:%8d  success:%8d  n0:%4d is_phi: %d  driver ang.:%9.5f  n_soln: %d  ang_nac:%9.5f  ang_side:%9.5f  len_on:%9.5f  len_ns:%9.5f\n", step, success, n0, mc.is_phi, step_size*rad2deg, *n_soln, ang_nac*rad2deg, ang_side*rad2deg, len_on, len_ns);

  if(write_out_pdb)
   {
    sprintf(out_pdb, "data/135l_loop_min.pdb");
    fprintf(STATUS, "\nRecording the solution at min. rmsd in %s\n", out_pdb);
    write_pdb_backbone(out_pdb, res_name, r_n, r_a, r_c, r_o, r_s, n0, n0+4);
   }
  //fprintf(STATUS,  "Yang move attempted \n");
  return;
 }

//========================================================================================================
void get_orgco(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3], double rorg_s[5][MSAR][3], int ns[5], int n0)
 {
  int cur_res, cur_atom;
  int i;
//--------------------------------------------------------------------------------------------------------
 
  for (i=0;i<5;i++)
    ns[i] = 0;
  
  for (cur_res=0; cur_res<5; cur_res++)
    for (cur_atom=res_atomno[n0-1]; cur_atom<res_atomno[n0+4]; cur_atom++)
      if(native[cur_atom].res_num == n0-1+cur_res)
        get_orgat(rorg_n, rorg_a, rorg_c, rorg_o, rorg_s, ns, cur_res, cur_atom);
   
  return;
 }

//========================================================================================================
void get_orgat(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3], double rorg_s[5][MSAR][3], int ns[5], int r, int a)
 {
  if(strcmp(native[a].atomname, "N")==0)
   {
    rorg_n[r][0] = orig_native[a].xyz.x;
    rorg_n[r][1] = orig_native[a].xyz.y;
    rorg_n[r][2] = orig_native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "CA")==0)
   {
    rorg_a[r][0] = orig_native[a].xyz.x;
    rorg_a[r][1] = orig_native[a].xyz.y;
    rorg_a[r][2] = orig_native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "C")==0)
   {
    rorg_c[r][0] = orig_native[a].xyz.x;
    rorg_c[r][1] = orig_native[a].xyz.y;
    rorg_c[r][2] = orig_native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    rorg_o[r][0] = orig_native[a].xyz.x;
    rorg_o[r][1] = orig_native[a].xyz.y;
    rorg_o[r][2] = orig_native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    rorg_o[r][0] = orig_native[a].xyz.x;
    rorg_o[r][1] = orig_native[a].xyz.y;
    rorg_o[r][2] = orig_native[a].xyz.z;
   }
  else
   {
    rorg_s[r][ns[r]][0] = orig_native[a].xyz.x;
    rorg_s[r][ns[r]][1] = orig_native[a].xyz.y;
    rorg_s[r][ns[r]][2] = orig_native[a].xyz.z;
    ns[r]++;
   }
  return;
 }   
//========================================================================================================
void get_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int n0)
 {
  int cur_res, cur_atom;
  int i;
//--------------------------------------------------------------------------------------------------------
 
  for (i=0;i<5;i++)
    ns[i] = 0;
  
  for (cur_res=0; cur_res<5; cur_res++)
    for (cur_atom=res_atomno[n0-1]; cur_atom<res_atomno[n0+4]; cur_atom++)
      if(native[cur_atom].res_num == n0-1+cur_res)
        get_atom(r_n, r_a, r_c, r_o, r_s, ns, cur_res, cur_atom);
   
  return;
 }

//========================================================================================================
void get_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int r, int a)
 {
  if(strcmp(native[a].atomname, "N")==0)
   {
    r_n[r][0] = native[a].xyz.x;
    r_n[r][1] = native[a].xyz.y;
    r_n[r][2] = native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "CA")==0)
   {
    r_a[r][0] = native[a].xyz.x;
    r_a[r][1] = native[a].xyz.y;
    r_a[r][2] = native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "C")==0)
   {
    r_c[r][0] = native[a].xyz.x;
    r_c[r][1] = native[a].xyz.y;
    r_c[r][2] = native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    r_o[r][0] = native[a].xyz.x;
    r_o[r][1] = native[a].xyz.y;
    r_o[r][2] = native[a].xyz.z;
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    r_o[r][0] = native[a].xyz.x;
    r_o[r][1] = native[a].xyz.y;
    r_o[r][2] = native[a].xyz.z;
   }
  else
   {
    r_s[r][ns[r]][0] = native[a].xyz.x;
    r_s[r][ns[r]][1] = native[a].xyz.y;
    r_s[r][ns[r]][2] = native[a].xyz.z;
    ns[r]++;
   }
  return;
 }   
//========================================================================================================
void put_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int n0)
 {
  int cur_res, cur_atom;
  int i;
//--------------------------------------------------------------------------------------------------------
 
  for (i=0;i<5;i++)
    ns[i] = 0;
  
  for (cur_res=0; cur_res<5; cur_res++)
    for (cur_atom=res_atomno[n0-1]; cur_atom<res_atomno[n0+4]; cur_atom++)
      if(native[cur_atom].res_num == n0-1+cur_res)
        put_atom(r_n, r_a, r_c, r_o, r_s, ns, cur_res, cur_atom);
   
  return;
 }

//========================================================================================================
void put_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], int r, int a)
 {
  if(strcmp(native[a].atomname, "N")==0)
   {
    native[a].xyz.x = r_n[r][0];
    native[a].xyz.y = r_n[r][1];
    native[a].xyz.z = r_n[r][2];
   }
  else if(strcmp(native[a].atomname, "CA")==0)
   {
    native[a].xyz.x = r_a[r][0];
    native[a].xyz.y = r_a[r][1];
    native[a].xyz.z = r_a[r][2];
   }
  else if(strcmp(native[a].atomname, "C")==0)
   {
    native[a].xyz.x = r_c[r][0];
    native[a].xyz.y = r_c[r][1];
    native[a].xyz.z = r_c[r][2];
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    native[a].xyz.x = r_o[r][0];
    native[a].xyz.y = r_o[r][1];
    native[a].xyz.z = r_o[r][2];
   }
  else if(strcmp(native[a].atomname, "O")==0)
   {
    native[a].xyz.x = r_o[r][0];
    native[a].xyz.y = r_o[r][1];
    native[a].xyz.z = r_o[r][2]; 
   }
  else
   {
    native[a].xyz.x = r_s[r][ns[r]][0];
    native[a].xyz.y = r_s[r][ns[r]][1];
    native[a].xyz.z = r_s[r][ns[r]][2];
    ns[r]++;
   }
  return;
 }   
//========================================================================================================
void yang_rotate(int n0, int is_phi)
 {
  int i, cur_atom;
  int another = 0;

  for (i=0;i<natoms;i++)
    yang_not_rotated[i] = 1;

  is_rotated[native_residue[n0].C]=1;
  yang_not_rotated[native_residue[n0].C]=0;
  yang_rotated_atoms[another++]=native_residue[n0].C;
  is_rotated[native_residue[n0].O]=2;
  yang_not_rotated[native_residue[n0].O]=0;
  yang_rotated_atoms[another++]=native_residue[n0].O;
  is_rotated[native_residue[n0+1].N]=2;
  yang_not_rotated[native_residue[n0+1].N]=0;
  yang_rotated_atoms[another++]=native_residue[n0+1].N;
  is_rotated[native_residue[n0+1].CA]=3;
  yang_not_rotated[native_residue[n0+1].CA]=0;
  yang_rotated_atoms[another++]=native_residue[n0+1].CA;
  is_rotated[native_residue[n0+1].C]=4;
  yang_not_rotated[native_residue[n0+1].C]=0;
  yang_rotated_atoms[another++]=native_residue[n0+1].C;
  is_rotated[native_residue[n0+1].O]=4;
  yang_not_rotated[native_residue[n0+1].O]=0;
  yang_rotated_atoms[another++]=native_residue[n0+1].O;
  is_rotated[native_residue[n0+2].N]=5;
  yang_not_rotated[native_residue[n0+2].N]=0;
  yang_rotated_atoms[another++]=native_residue[n0+2].N;

  if(is_phi){
    is_rotated[native_residue[n0+2].CA]=6;
    yang_not_rotated[native_residue[n0+2].CA]=0;
    yang_rotated_atoms[another++]=native_residue[n0+2].CA;
    is_rotated[native_residue[n0+2].C]=6;
    yang_not_rotated[native_residue[n0+2].C]=0;
    yang_rotated_atoms[another++]=native_residue[n0+2].C;
    is_rotated[native_residue[n0+2].O]=6;
    yang_not_rotated[native_residue[n0+2].O]=0;
    yang_rotated_atoms[another++]=native_residue[n0+2].O;
  } else {
    is_rotated[native_residue[n0-1].O]=6;
    yang_not_rotated[native_residue[n0-1].O]=0;
    yang_rotated_atoms[another++]=native_residue[n0-1].O;
    is_rotated[native_residue[n0].N]=6;
    yang_not_rotated[native_residue[n0].N]=0;
    yang_rotated_atoms[another++]=native_residue[n0].N;
    is_rotated[native_residue[n0].CA]=6;
    yang_not_rotated[native_residue[n0].CA]=0;
    yang_rotated_atoms[another++]=native_residue[n0].CA;
  }

  for (cur_atom=res_atomno[n0];cur_atom<res_atomno[n0+3];cur_atom++){
    if(native[cur_atom].is_sidechain){
      is_rotated[cur_atom]=2*(native[cur_atom].res_num-n0)+1;
      yang_not_rotated[cur_atom]=0;
      yang_rotated_atoms[another++]=cur_atom;
    }
  }

  yang_rotated_natoms = another;
//  fprintf(STATUS, "yang: %5d\n", yang_rotated_natoms);
    
  return;
 }
//========================================================================================================
void yangloop(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int ns[5], double b_len[6], double b_ang[7], int n0, float dih_ch, char res_name[5][4], int *n_soln)
 {
  int i, j, k, z;
  double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3], r_soln_o[max_soln][3][3], r_soln_s[3][MSAR][3];
  double rmsd, sum, dr[3];
  double r0_n[3][3], r0_a[3][3], r0_c[3][3], r0_o[3][3], r0_s[3][MSAR][3];
  double r0drms_n[3][3], r0drms_a[3][3], r0drms_c[3][3];
  double rot1_n[3][3], rot1_a[3][3], rot1_c[3][3], rot1_o[3][3], rot1_s[3][MSAR][3];
  double rot2_n[3][3], rot2_a[3][3], rot2_c[3][3], rot2_o[3][3], rot2_s[3][MSAR][3];
  double rot3_a[3][3], rot3_c[3][3], rot3_o[3][3], rot3_s[3][MSAR][3];
  double rot4_c[3][3], rot4_s[3][MSAR][3];
  double r1_n[5][3], r1_a[5][3], r1_c[5][3], r1_o[5][3];
  
  int  calc_rmsd=1;
  double t_ang[2];
//  double len_f[6], ang_f[7];

  int index = 1000;
  double min_value = 1000.;
  int n_soln_before = 0;
  double r_n_before[3][3], r_ca_before[3][3], r_c_before[3][3];
  double r_n_after[3][3], r_ca_after[3][3], r_c_after[3][3];
  double before_jacobi, after_jacobi;
//  double oang_i, oang_f;
//  double olen_i[3], olen_f[3];
//  double slen_i[3], slen_f[3];
//!-----------------------------------------------------------------------
  *n_soln = 0;
/*  
  c_bnd_len(r_c[1], r_o[1], &olen_i[0]);
  c_bnd_len(r_c[2], r_o[2], &olen_i[1]);
  c_bnd_len(r_c[3], r_o[3], &olen_i[2]);
  c_bnd_len(r_a[1], r_s[1][0], &slen_i[0]);
  c_bnd_len(r_a[2], r_s[2][0], &slen_i[1]);
  c_bnd_len(r_a[3], r_s[3][0], &slen_i[2]);
*/
  t_ang[0] = pi; //The omega dihedral angle, always has value 180 degrees
  t_ang[1] = pi;

	//First, we solve the loop closure problem for the current set of dihedrals, pre rotation
  initialize_loop_closure(b_len, b_ang, t_ang);  //Feed this function values that are to be fixed (by basic chemistry), namely bond angles, bond lengths, and the omega dihedral angle, always 180 degrees
  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, &n_soln_before); //This is the function that solves the loop closture problem..Need to feed it coordinates for nitrogen/C_alpha atoms in first and third residue...see tripep_closure.h
  soln_no_before = n_soln_before;  //Although n_soln_before is initialized to zero, its value  gets updated in the above solve_3p3_poly
  if(soln_no_before==0)
    soln_no_before=1;
  for (i=0;i<3;i++)
    for (k=0;k<3;k++) //Keep track of the atom positions that we just computed, prior to rotation
     {
      r_n_before[i][k] = r_n[i+1][k];
      r_ca_before[i][k] = r_a[i+1][k];
      r_c_before[i][k] = r_c[i+1][k];
     }
  loop_Jacobian(r_n_before, r_ca_before, r_c_before, &before_jacobi); //Compute Jacobi pre-rotation
  jacobi_before = before_jacobi;

  for (i=0;i<3;i++)  //Not totally sure what r0drms is..
    for (j=0;j<3;j++)
     {
      r0drms_n[i][j]=r_n[i+1][j];
      r0drms_a[i][j]=r_a[i+1][j];
      r0drms_c[i][j]=r_c[i+1][j];
     }
  //
// driver angle
// We do a driver rotation
// Note that the second argument of the function driver_rot is coordiate of
//the target atom AFTER rotation, and this one is computed in the driver_rot function
//See tripep_closure.h for what all these arguments mean

/*I THINK the general idea is that the "driver" function applies a rotation either to FINAL
phi angle of the triplet, or to the FIRST psi angle
If the former, then we change positions of final C and Ca atoms
If the later, we change position of first N and CA (although I had somehow thoguht these were fixed..)
*/
  if(mc.is_phi) //We rotate the phi dihedral, which goes from N to CA
   {
    driver_rot(r_c[3], r1_c[3], r_a[4], r_n[4], dih_ch); //We plug in the current value of the carbon, and get the updated value (second argument) as a result of a dihedral rotation of value dih_ch...honeslty not totally clear how this works from Dill paper..
    driver_rot(r_a[3], r1_a[3], r_a[4], r_n[4], dih_ch);
    for (i=0;i<3;i++)
     {
      r_c[3][i] = r1_c[3][i];
      r_a[3][i] = r1_a[3][i];
     }
   }
  else //psi
   {
    driver_rot(r_o[0], r1_o[0], r_a[0], r_c[0], dih_ch);
    driver_rot(r_n[1], r1_n[1], r_a[0], r_c[0], dih_ch);
    driver_rot(r_a[1], r1_a[1], r_a[0], r_c[0], dih_ch);
    for (i=0;i<3;i++)
     {
      r_o[0][i] = r1_o[0][i];
      r_n[1][i] = r1_n[1][i];
      r_a[1][i] = r1_a[1][i];
     }
   }

  for (i=0;i<3;i++)  //not totally sure what these for loops are doing
    for (j=0;j<3;j++)
     {
      r0_n[i][j]=r_n[i+1][j];
      r0_a[i][j]=r_a[i+1][j];
      r0_c[i][j]=r_c[i+1][j];
      r0_o[i][j]=r_o[i+1][j];
      for (z=0;z<ns[i+1];z++)
        r0_s[i][z][j]=r_s[i+1][z][j];
     }

  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, n_soln); //Find new loop closure solution post rotation

     if (calc_rmsd)  //This clause is ridiculous since calc_rmsd is initialized to 1 before so will always be true...but it's an important clause since in it, you choose the solution to accept etc...this can be cleaned up but I didn't want to mess with code too much...
      {
       for (k=0;k<*n_soln;k++)
        {
         for (i=0;i<3;i++)
	   for (j=0;j<3;j++)
	    {
	     r_n[i+1][j] = r_soln_n[k][i][j];
	     r_a[i+1][j] = r_soln_a[k][i][j];
	     r_c[i+1][j] = r_soln_c[k][i][j];
	    } 
	   sum = 0.0e0;
	   
         for (i=0;i<3;i++)
	 {
	   for (j=0;j<3;j++)
	     dr[j] = r_soln_n[k][i][j] - r0drms_n[i][j];
	   sum += dot_product(dr, dr);
	   for (j=0;j<3;j++)
	     dr[j] = r_soln_a[k][i][j] - r0drms_a[i][j];
	   sum += dot_product(dr, dr);
	   for (j=0;j<3;j++)
	     dr[j] = r_soln_c[k][i][j] - r0drms_c[i][j];
	   sum += dot_product(dr, dr);
	  }
         rmsd = sqrt(sum/9.0e0);

	 if(rmsd<min_value)
	  {
	   min_value = rmsd;
	   index = k;
	  }
        }
//       fprintf(STATUS, "\nsoln no.: %d   min. rmsd:%7.3f\n", index+1, min_value);
       if(*n_soln>0)
        {
         index = (int) (threefryrand()*(*n_soln));
         if(index>=(*n_soln))
          {
           fprintf(STATUS, "ERROR!!!: index %d is greater than equal to no. of solutions %d.\n", index, *n_soln);
           exit(1);
          }
        }
        //AB decided to move the Jacobi calculation to later, at which point there is a loop that naturlaly allows r_n_after etc to be defined

//       if(*n_soln>0)
       if((*n_soln > 0)&&(*n_soln <= deg_pol)) //This is a huge clause, and I did not look at what it does, althoguh I should
        {
	 if (mc.is_phi)
// This is the case where driver angle is res. 4. (0~4)
	  {
// rotation with n0-a0
           for (z=0;z<ns[1];z++)
             get_rot(r0_s[0][z], r_soln_s[0][z], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);
           get_rot(r0_o[0], rot1_o[0], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);
           get_rot(r0_n[1], rot1_n[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);
           get_rot(r0_a[1], rot1_a[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);
           for (z=0;z<ns[2];z++)
             get_rot(r0_s[1][z], rot1_s[1][z], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);
           get_rot(r0_c[1], rot1_c[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], r0_c[0], r_soln_c[index][0]);

// rotation with a0-c0	 
	   get_rot(rot1_o[0], r_soln_o[index][0], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot1_n[1], r_soln_n[index][1]);
	   get_rot(rot1_a[1], rot2_a[1], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot1_n[1], r_soln_n[index][1]);
           for (z=0;z<ns[2];z++)
	   get_rot(rot1_s[1][z], rot2_s[1][z], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot1_n[1], r_soln_n[index][1]);
	   get_rot(rot1_c[1], rot2_c[1], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot1_n[1], r_soln_n[index][1]);

// rotation with c0-n1
           for (z=0;z<ns[2];z++)
	     get_rot(rot2_s[1][z], rot3_s[1][z], r_soln_a[index][0], r_soln_c[index][0], r_soln_n[index][1], rot2_a[1], r_soln_a[index][1]);
	   get_rot(rot2_c[1], rot3_c[1], r_soln_a[index][0], r_soln_c[index][0], r_soln_n[index][1], rot2_a[1], r_soln_a[index][1]);

// rotation with n1-a1
           for (z=0;z<ns[2];z++)
	     get_rot(rot3_s[1][z], r_soln_s[1][z], r_soln_c[index][0], r_soln_n[index][1], r_soln_a[index][1], rot3_c[1], r_soln_c[index][1]);
	 
// rotation with a[+]-n[+]
           driver_rot(r0_o[2], r_soln_o[index][2], r_a[4], r_n[4], dih_ch);
           for (z=0;z<ns[3];z++)
             driver_rot(r0_s[2][z], rot1_s[2][z], r_a[4], r_n[4], dih_ch);
           driver_rot(r0_n[2], rot1_n[2], r_a[4], r_n[4], dih_ch);
           driver_rot(r0_c[1], rot1_c[1], r_a[4], r_n[4], dih_ch);
           driver_rot(r0_o[1], rot1_o[1], r_a[4], r_n[4], dih_ch);
           driver_rot(r0_a[1], rot1_a[1], r_a[4], r_n[4], dih_ch);

// rotation with c2-a2	 
           for (z=0;z<ns[3];z++)
             get_rot(rot1_s[2][z], r_soln_s[2][z], r_n[4], r_soln_c[index][2], r_soln_a[index][2], rot1_n[2], r_soln_n[index][2]);
           get_rot(rot1_c[1], rot2_c[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], rot1_n[2], r_soln_n[index][2]);
           get_rot(rot1_o[1], rot2_o[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], rot1_n[2], r_soln_n[index][2]);
           get_rot(rot1_a[1], rot2_a[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], rot1_n[2], r_soln_n[index][2]);

// rotation with a2-n2
           get_rot(rot2_o[1], rot3_o[1], r_soln_c[index][2], r_soln_a[index][2], r_soln_n[index][2], rot2_c[1], r_soln_c[index][1]);
           get_rot(rot2_a[1], rot3_a[1], r_soln_c[index][2], r_soln_a[index][2], r_soln_n[index][2], rot2_c[1], r_soln_c[index][1]);
	   
// rotation with n2-c1
	   get_rot(rot3_o[1], r_soln_o[index][1], r_soln_a[index][2], r_soln_n[index][2], r_soln_c[index][1], rot3_a[1], r_soln_a[index][1]);
	  }
	 else
// This is the case where driver angle is res. 0. (0~4)
	  {
// rotation with n[-]-a[-]
           for (z=0;z<ns[1];z++)
             driver_rot(r0_s[0][z], rot1_s[0][z], r_a[0], r_c[0], dih_ch);
           driver_rot(r0_c[0], rot1_c[0], r_a[0], r_c[0], dih_ch);
           driver_rot(r0_o[0], rot1_o[0], r_a[0], r_c[0], dih_ch);
           driver_rot(r0_n[1], rot1_n[1], r_a[0], r_c[0], dih_ch);
           driver_rot(r0_a[1], rot1_a[1], r_a[0], r_c[0], dih_ch);
           for (z=0;z<ns[2];z++)
             driver_rot(r0_s[1][z], rot1_s[1][z], r_a[0], r_c[0], dih_ch);
           driver_rot(r0_c[1], rot1_c[1], r_a[0], r_c[0], dih_ch);

// rotation with n0-a0
           for (z=0;z<ns[1];z++)
             get_rot(rot1_s[0][z], r_soln_s[0][z], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);
           get_rot(rot1_o[0], rot2_o[0], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);
           get_rot(rot1_n[1], rot2_n[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);
           get_rot(rot1_a[1], rot2_a[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);
           for (z=0;z<ns[2];z++)
             get_rot(rot1_s[1][z], rot2_s[1][z], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);
           get_rot(rot1_c[1], rot2_c[1], r_c[0], r_soln_n[index][0], r_soln_a[index][0], rot1_c[0], r_soln_c[index][0]);

// rotation with a0-c0	 
	   get_rot(rot2_o[0], r_soln_o[index][0], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot2_n[1], r_soln_n[index][1]);
	   get_rot(rot2_a[1], rot3_a[1], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot2_n[1], r_soln_n[index][1]);
           for (z=0;z<ns[2];z++)
	     get_rot(rot2_s[1][z], rot3_s[1][z], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot2_n[1], r_soln_n[index][1]);
	   get_rot(rot2_c[1], rot3_c[1], r_soln_n[index][0], r_soln_a[index][0], r_soln_c[index][0], rot2_n[1], r_soln_n[index][1]);

// rotation with c0-n1
           for (z=0;z<ns[2];z++)
	     get_rot(rot3_s[1][z], rot4_s[1][z], r_soln_a[index][0], r_soln_c[index][0], r_soln_n[index][1], rot3_a[1], r_soln_a[index][1]);
	   get_rot(rot3_c[1], rot4_c[1], r_soln_a[index][0], r_soln_c[index][0], r_soln_n[index][1], rot3_a[1], r_soln_a[index][1]);

// rotation with n1-a1
           for (z=0;z<ns[2];z++)
	     get_rot(rot4_s[1][z], r_soln_s[1][z], r_soln_c[index][0], r_soln_n[index][1], r_soln_a[index][1], rot4_c[1], r_soln_c[index][1]);
	 
           for (i=0;i<3;i++)
             r_soln_o[index][2][i] = r0_o[2][i];

// rotation with c2-a2	 
           for (z=0;z<ns[3];z++)
             get_rot(r0_s[2][z], r_soln_s[2][z], r_n[4], r_soln_c[index][2], r_soln_a[index][2], r0_n[2], r_soln_n[index][2]);
           get_rot(r0_c[1], rot1_c[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], r0_n[2], r_soln_n[index][2]);
           get_rot(r0_o[1], rot1_o[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], r0_n[2], r_soln_n[index][2]);
           get_rot(r0_a[1], rot1_a[1], r_n[4], r_soln_c[index][2], r_soln_a[index][2], r0_n[2], r_soln_n[index][2]);

// rotation with a2-n2
           get_rot(rot1_o[1], rot2_o[1], r_soln_c[index][2], r_soln_a[index][2], r_soln_n[index][2], rot1_c[1], r_soln_c[index][1]);
           get_rot(rot1_a[1], rot2_a[1], r_soln_c[index][2], r_soln_a[index][2], r_soln_n[index][2], rot1_c[1], r_soln_c[index][1]);
	   
// rotation with n2-c1
	   get_rot(rot2_o[1], r_soln_o[index][1], r_soln_a[index][2], r_soln_n[index][2], r_soln_c[index][1], rot2_a[1], r_soln_a[index][1]);
	  }

         for (i=0;i<3;i++)
	   for (j=0;j<3;j++)
	    {
	    //The first three are added by AB on July 15 2019
	     r_n_after[i][j] = r_soln_n[index][i][j];
	     r_ca_after[i][j] = r_soln_a[index][i][j];
	     r_c_after[i][j] = r_soln_c[index][i][j];
	     //end AB
	     
	    
	     r_n[i+1][j] = r_soln_n[index][i][j];
	     r_a[i+1][j] = r_soln_a[index][i][j];
	     r_c[i+1][j] = r_soln_c[index][i][j];
	     r_o[i+1][j] = r_soln_o[index][i][j];
             for (z=0;z<ns[i+1];z++)
	       r_s[i+1][z][j] = r_soln_s[i][z][j];
	    } 
	 
        }
      }


       loop_Jacobian(r_n_after, r_ca_after, r_c_after, &after_jacobi);
       jacobi_after = after_jacobi;


  return;
 }

//========================================================================================================

void get_template()
 {
  FILE *the_file;

  char line[250], temp[10];
  int initial_template[MAX_RES], mod_template[MAX_RES], len_template[MAX_RES];
  int i, j;
  int prev_tmpl = 0;
  int len_tmpl = 0;

  for (i=0;i<nresidues;i++)
   {
    initial_template[i] = 0;
    len_template[i] = 0;
   }

  if((the_file = fopen(template_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", template_file);
    exit(1);
   }
  while(fgets(line, 100, the_file) != NULL)
   {
    if (strncmp(line, "ATOM", 4) == 0)
     {
      strncpy(temp, &(line[22]), 4);
      strcpy(&temp[4], "\0");
      initial_template[atoi(temp)-1] = 1;
     }
   }
  fclose(the_file);

  for (i=0;i<nresidues;i++)
   {
    if ((prev_tmpl==0)&&(initial_template[i]==1))
      len_tmpl = 1;
    else if ((prev_tmpl==1)&&(initial_template[i]==1))
      len_tmpl++;
    else if ((prev_tmpl==1)&&(initial_template[i]==0))
     {
      for (j=i-len_tmpl;j<i;j++)
	len_template[j] = len_tmpl;
      len_tmpl = 0;
     }
    if ((i==nresidues-1)&&(initial_template[i]==1))
      for (j=i-len_tmpl+1;j<=i;j++)
	len_template[j] = len_tmpl;
    prev_tmpl = initial_template[i];
   }
  
  for (i=0;i<nresidues;i++)
    if (len_template[i]<MIN_LEN_TEMPLATE)
      mod_template[i] = 0;
    else
      mod_template[i] = initial_template[i];

  prev_tmpl = 0;
  for (i=0;i<nresidues;i++)
   {
    if ((prev_tmpl==0)&&(mod_template[i]==1)&&(i!=0))
      is_template[i] = 0;
    else if ((prev_tmpl==1)&&(mod_template[i]==0))
      is_template[i-1] = 0;
    else
      is_template[i] = mod_template[i];
    prev_tmpl = mod_template[i];
   }
  for (i=0;i<nresidues;i++){
    fprintf(STATUS, "%5d%5d%5d%5d%5d\n", i, initial_template[i], len_template[i], mod_template[i], is_template[i]);
    fflush(STATUS);
  }
  return;
 }

//========================================================================================================
void check_bb()
 {
  struct vector tmp1, tmp2, tmp3, tmp4, plane1, plane2, bisect1, bisect2;
  int i;
   
  for (i=0;i<nresidues-2;i++)
   {
    MakeVector(native[native_residue[i].CA].xyz, native[native_residue[i].N].xyz, &tmp1);
    MakeVector(native[native_residue[i].CA].xyz, native[native_residue[i].O].xyz, &tmp2);
    MakeVector(native[native_residue[i+2].CA].xyz, native[native_residue[i+2].N].xyz, &tmp3);
    MakeVector(native[native_residue[i+2].CA].xyz, native[native_residue[i+2].O].xyz, &tmp4);
    CrossProduct(tmp1, tmp2, &plane1);
    CrossProduct(tmp3, tmp4, &plane2);
    a_PCA[i] = Angle(plane1, plane2);
    a_PCA[i] *= rad2deg;
    Normalize(&tmp1);
    Normalize(&tmp2);
    Normalize(&tmp3);
    Normalize(&tmp4);
    bisect(tmp1, tmp2, &bisect1);
    bisect(tmp3, tmp4, &bisect2);
    a_bCA[i] = Angle(bisect1, bisect2);
    a_bCA[i] *= rad2deg;

    MakeVector(native[native_residue[i].C].xyz, native[native_residue[i+1].N].xyz, &tmp1);
    MakeVector(native[native_residue[i+1].N].xyz, native[native_residue[i+1].CA].xyz, &tmp2);
    MakeVector(native[native_residue[i+1].CA].xyz, native[native_residue[i+1].C].xyz, &tmp3);
    phim[i] = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    phim[i] *= rad2deg;
    MakeVector(native[native_residue[i+1].N].xyz, native[native_residue[i+1].CA].xyz, &tmp1);
    MakeVector(native[native_residue[i+1].CA].xyz, native[native_residue[i+1].C].xyz, &tmp2);
    MakeVector(native[native_residue[i+1].C].xyz, native[native_residue[i+2].N].xyz, &tmp3);
    psim[i] = struct_calc_dih_ang(tmp1, tmp2, tmp3);
    psim[i] *= rad2deg;
    phim[i] += 180.;
    psim[i] += 180.;

//    fprintf(STATUS, "%4d %4d %9.3f %9.3f %9.3f %9.3f %9.3f\n", i+1, i+2,  dih_CA[i], ang_CA[i], ang_CA[i+1], len12_CA[i], len03_CA[i]);
   }
  
  return;
 }

//========================================================================================================
void initialize_torsion()
 {
  int r, cur_res, value;
  int x, y, z, w;

  FILE *ftor;
  char line[1000];
      
  for (r=0;r<nresidues;r++)
    for (x=0;x<6;x++)
      for (y=0;y<6;y++)
        for (z=0;z<6;z++)
          for (w=0;w<6;w++)
            torsion_E[r][x][y][z][w] = 1000;

  fprintf(STATUS, "Opening the file: %s\n", triplet_file);
  if((ftor=fopen(triplet_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", triplet_file);
    exit(1);
   }
  while(fgets(line, 200, ftor)!=NULL)
   {
    sscanf(line, "%d %*d%*d%*d %d %d %d %d %d %*s%*s", &cur_res, &x, &y, &z, &w, &value); 
    torsion_E[cur_res][x][y][z][w] = value;
   }
  fclose(ftor);
   
  return;
 }

//========================================================================================================
void initialize_sct()
 {
  int r, cur_res, value;
  int x, y, z, w;

  FILE *ftor;
  char line[1000];

  for (r=0;r<nresidues;r++)
    for (x=0;x<12;x++)
      for (y=0;y<12;y++)
        for (z=0;z<12;z++)
          for (w=0;w<12;w++)
            sct_E[r][x][y][z][w] = 1000;

  fprintf(STATUS, "Opening the file: %s\n", sctorsion_file);
  if((ftor=fopen(sctorsion_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", sctorsion_file);
    exit(1);
   }
  while(fgets(line, 200, ftor)!=NULL)
   {
    sscanf(line, "%d %*d%*d%*d %d %d %d %d %d %*s%*s", &cur_res, &x, &y, &z, &w, &value);
    sct_E[cur_res][x][y][z][w] = value;
//    fprintf(STATUS, "%10d\n", sct_E[cur_res][x][y][z][w]);
   }
  fclose(ftor);

  return;
 }

//========================================================================================================
void initialize_aromatic()
 {
  int value;
  int x, r;

  FILE *ftor;
  char line[1000];

  fprintf(STATUS, "Opening the file: %s\n", aromatic_file);
  if((ftor=fopen(aromatic_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", aromatic_file);
    exit(1);
   }
  while(fgets(line, 200, ftor)!=NULL)
   {
    sscanf(line, "%d %d %*s%*s", &x, &value); 
    aromatic_E[x] = value;
//    fprintf(STATUS, "%2d %5d\n", x, aromatic_E[x]);
   }
  fclose(ftor);

  for (r=0;r<nresidues;r++)
    if((strcmp(native_residue[r].res, "PHE")==0) || (strcmp(native_residue[r].res, "TRP")==0) ||
       (strcmp(native_residue[r].res, "RING")==0))
     {
      Res_aromatic[Naromatic] = r;
//      fprintf(STATUS, "%5d %5d\n", Naromatic, Res_aromatic[Naromatic]);
      Naromatic++;
     }
  return;
 }

//========================================================================================================
int will2edo(int edo) {
  
  if (edo == 1)
    return 7;
  else if (edo == 2)
    return 0;
  else if (edo == 3)
    return 19;
  else if (edo == 4)
    return 10;
  else if (edo == 5)
    return 9;
  else if (edo == 6)
    return 15;
  else if (edo == 7)
    return 16;
  else if (edo == 8)
    return 3;
  else if (edo == 9)
    return 6;
  else if (edo == 10)
    return 2;
  else if (edo == 11)
    return 5;
  else if (edo == 12)
    return 11;
  else if (edo == 13)
    return 1;
  else if (edo == 14)
    return 4;
  else if (edo == 15)
    return 12;
  else if (edo == 16)
    return 8;
  else if (edo == 17)
    return 13;
  else if (edo == 18)
    return 18;
  else if (edo == 19)
    return 17;
  else if (edo == 20)
    return 14;
  else
    return 0;

};

//========================================================================================================
void read_cluster()
 {
  FILE *fcluster;
  int cur_res;
//  int i, j;
  int temp_value=0;
  float x, y;
  char cluster_file[100], line[1000];
  
  sprintf(cluster_file, "../config_files/center%d_state", NOCLUSTERS); 
  fprintf(STATUS, "Opening the file: %s\n", cluster_file);
  if((fcluster=fopen(cluster_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", cluster_file);
    exit(1);
   }
  while(fgets(line, 200, fcluster)!=NULL)
   {
    temp_value = temp_value%NOCLUSTERS;
    sscanf(line, "%d %f %f", &cur_res, &x, &y); 
    cluster_phi[will2edo(cur_res)][temp_value] = x;
    cluster_psi[will2edo(cur_res)][temp_value] = y;
    temp_value++;
   }
  fclose(fcluster);
//  for (i=0;i<20;i++)
//    for (j=0;j<NOCLUSTERS;j++)
//      fprintf(STATUS, "%d %d %8.5f %8.5f\n", i, j, cluster_phi[i][j], cluster_psi[i][j]);
  return;
 }

//========================================================================================================
void check_phipsi()
 {
  double res_CA[MAXSEQUENCE][3];
  double res_N[MAXSEQUENCE][3];
  double res_C[MAXSEQUENCE][3];
  int i;
   
  for (i=0;i<nresidues;i++)
   {
    res_N[i][0] = native[native_residue[i].N].xyz.x; 
    res_N[i][1] = native[native_residue[i].N].xyz.y; 
    res_N[i][2] = native[native_residue[i].N].xyz.z; 
    res_CA[i][0] = native[native_residue[i].CA].xyz.x; 
    res_CA[i][1] = native[native_residue[i].CA].xyz.y; 
    res_CA[i][2] = native[native_residue[i].CA].xyz.z; 
    res_C[i][0] = native[native_residue[i].C].xyz.x; 
    res_C[i][1] = native[native_residue[i].C].xyz.y; 
    res_C[i][2] = native[native_residue[i].C].xyz.z; 
   }
   
  for (i=0;i<nresidues-2;i++)
   {
    c_dih_ang(res_C[i], res_N[i+1], res_CA[i+1], res_C[i+1], &cur_phi[i]);
    c_dih_ang(res_N[i+1], res_CA[i+1], res_C[i+1], res_N[i+2], &cur_psi[i]);
    cur_phi[i] *= rad2deg;
    cur_psi[i] *= rad2deg;
//    fprintf(STATUS, "%4d %4d %9.3f %9.3f %9.3f %9.3f %9.3f\n", i+1, i+2,  dih_CA[i], ang_CA[i], ang_CA[i+1], len12_CA[i], len03_CA[i]);
   }
  
  return;
 }

//========================================================================================================
void initialize_secstr()
 {
  FILE *ftor;
  char line[1000];
  char iscorrect[600], secstr_org[600];
  char tmp_str[10];
  int i;
      
  fprintf(STATUS, "Opening the file: %s\n", sec_str_file);
  if((ftor=fopen(sec_str_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", sec_str_file);
    exit(1);
   }
  fgets(line, 500, ftor);
  sscanf(line, "%s", iscorrect); 
  fgets(line, 500, ftor);
  sscanf(line, "%s", secstr_org); 
  fclose(ftor);

  for (i=0;i<nresidues;i++)
   {
    sprintf(tmp_str, "%c", iscorrect[i]);
    if((atoi(tmp_str) >= CUT_SECSTR)&&secstr_org[i]=='H')
      secstr[i] = secstr_org[i];
    else if((atoi(tmp_str) >= CUT_SECSTR)&&secstr_org[i]=='E')
      secstr[i] = secstr_org[i];
    else if((atoi(tmp_str) >= CUT_SECSTR)&&secstr_org[i]=='C')
      secstr[i] = 'L';
    else
      secstr[i] = 'C';
   }
  fprintf(STATUS, "\n");
  //fprintf(STATUS, "secondary structure correct?\n%s\n", iscorrect);
  //fprintf(STATUS, "input secondary structure:\n%s\n", secstr_org);
  fprintf(STATUS, "secondary structure:\n%s\n", secstr);

  return;
 }

