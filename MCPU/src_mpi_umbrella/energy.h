Float FullAtomEnergy();
void ResetEnergies(long int check);
Float HydrogenBonds();
Float torsionenergy();
Float sctenergy();
void aromatic_center(int, struct vector *);
float aromatic_plane(int, int, struct vector *, struct vector *);
Float aromaticenergy();

Float Rgenergy();

void ResetEnergies(long int check) {

  TypeContacts();
  E_pot = FullAtomEnergy();
  E_tor = torsionenergy();
  E_sct = sctenergy();
  E_aro = aromaticenergy();
  E_hbond = HydrogenBonds();
  
  	 if (strcmp(constraint_file, "None")!=0 ){ //AB
  E_constraint =  Compute_constraint_energy(native_residue, native); //AB
  }	else { //AB
  E_constraint =0; //AB
   } //AB

  E = weight_potential*E_pot + weight_clash*nclashes + weight_hbond*E_hbond + TOR_WEIGHT*E_tor + SCT_WEIGHT*E_sct + ARO_WEIGHT*E_aro + E_constraint;   

  if (abs(E - prev_E) > 1.0 && check) {
    fprintf(STATUS, "ResetEnergies(%ld): prev_E(%.5f), New E(%.5f) dE = %.5f\n", check, prev_E, E, E-prev_E);
    fprintf(STATUS, "E_pot %.5f %.5f, E_hbond %.5f %.5f, E_tor %.5f %.5f, E_sct %.5f %.5f, E_aro %.5f %.5f E_constraint %.5f %.5f \n", prev_E_pot, E_pot, prev_E_hbond, E_hbond, prev_E_tor, E_tor, prev_E_sct, E_sct, prev_E_aro, E_aro, prev_E_constraint, E_constraint); 
  }

  prev_E_pot = E_pot;
  prev_E_tor = E_tor;
  prev_E_sct = E_sct;
  prev_E_aro = E_aro;
  prev_E_hbond = E_hbond;
  
  prev_E_constraint = E_constraint; //AB

  prev_E = E;
}

Float FullAtomEnergy() {

  Float e = 0;
  int i, j;
  
  for (i=0; i<MAX_TYPES; i++) {
    for (j=i; j<MAX_TYPES; j++) {
      e += type_contacts[i][j]*potential[i][j];
//      if(((i==bb_O_type)&&(j==bb_N_type))||((i==bb_N_type)&&(j==bb_O_type)))
//	fprintf(STATUS, "%d %d %d %5.2f\n", i, j, type_contacts[i][j], potential[i][j]);
    }
  }
  return e;

}

Float torsionenergy()
 {
  int cur_res, x, y, z, w;
  float locang_int = 30.;
  float loctor_int = 60.;
  
  float e = 0;
  int i;

  check_bb();
  for (i=0;i<nresidues-2;i++)
   {
      cur_res = i;
      x=(int)((a_PCA[i])/locang_int);
      y=(int)((a_bCA[i])/locang_int);
      z=(int)((phim[i])/loctor_int);
      w=(int)((psim[i])/loctor_int);
      e += torsion_E[cur_res][x][y][z][w];
   }
   
  e /= 1000.;
  return e;
 }

Float sctenergy()
 {
  int cur_res;
  int i_ang[4];
  float locang_int = 30.;
  float ang[4];

  float e = 0;
  int i, j;

  for (i=0;i<nresidues-2;i++)
   {
     if(native_residue[i+1].ntorsions==0)
       continue;
      cur_res = i;
      for (j=0;j<4;j++)
        i_ang[j] = 0;
      for (j=0;j<native_residue[i+1].ntorsions;j++)
       {
        ang[j] = native_residue[i+1].tmpchi[j]*rad2deg;
        if(ang[j]<-180.)
          while(ang[j]<-180.)
            ang[j] += 360.;
        else if(ang[j]>=180.)
          while(ang[j]>=180.)
            ang[j] -= 360.;
        if(ang[j]<0)
          ang[j] += +180.+0.00001;
        else
          ang[j] += +180.-0.00001;
        ang[j] += 15.;
        ang[j] = fmod(ang[j], 360.);
//        fprintf(STATUS, "%3d %3d %9.3f\n", i, j, ang[j]);
        i_ang[j]=(int)(ang[j]/locang_int);
       }
      e += sct_E[cur_res][i_ang[0]][i_ang[1]][i_ang[2]][i_ang[3]];
   }

  e /= 1000.;
  return e;
 }

Float aromaticenergy()
 {
  float e = 0;
  int r1, r2;
  int i, j, x;
  struct vector vector_1, vector_2, vector_center;
  struct vector plane_1, plane_2;
  float distance_2, plane_angle;
  float aro_int = 10.;

//  if(mcstep>=1)
//    exit(0);
  for (r1=0;r1<Naromatic;r1++)
    for (r2=r1+1;r2<Naromatic;r2++)
     {
      i = Res_aromatic[r1];
      j = Res_aromatic[r2];
      aromatic_center(i, &vector_1);
      aromatic_center(j, &vector_2);
      MakeVector(vector_1, vector_2, &vector_center);
      distance_2 = D2(vector_1, vector_2);
      if(distance_2 < (AROMATIC_DISTANCE*AROMATIC_DISTANCE))
       {
        plane_angle = aromatic_plane(i, j, &plane_1, &plane_2);
        if(plane_angle > 89.9)
          plane_angle = 89.9;
        x=(int)(plane_angle/aro_int);
        e += aromatic_E[x];
//        fprintf(STATUS, "%3d %3d %3d %3d %8.3f %2d %5d %8.3f\n", r1, r2, i, j, plane_angle, x, aromatic_E[x], e);
       }
     }

  e /= 1000.;
  return e;
 }
  
void aromatic_center(int res_no, struct vector *V)
 {
  if((strcmp(native_residue[res_no].res, "PHE")==0)||(strcmp(native_residue[res_no].res, "RING")==0))
   {
    (*V).x = (native[native_residue[res_no].CG].xyz.x + native[native_residue[res_no].CE1].xyz.x
           + native[native_residue[res_no].CE2].xyz.x)/3.;
    (*V).y = (native[native_residue[res_no].CG].xyz.y + native[native_residue[res_no].CE1].xyz.y
           + native[native_residue[res_no].CE2].xyz.y)/3.;
    (*V).z = (native[native_residue[res_no].CG].xyz.z + native[native_residue[res_no].CE1].xyz.z
           + native[native_residue[res_no].CE2].xyz.z)/3.;
   }
  else if(strcmp(native_residue[res_no].res, "TRP")==0)
   {
    (*V).x = (native[native_residue[res_no].CG].xyz.x + native[native_residue[res_no].CZ2].xyz.x
           + native[native_residue[res_no].CZ3].xyz.x)/3.;
    (*V).y = (native[native_residue[res_no].CG].xyz.y + native[native_residue[res_no].CZ2].xyz.y
           + native[native_residue[res_no].CZ3].xyz.y)/3.;
    (*V).z = (native[native_residue[res_no].CG].xyz.z + native[native_residue[res_no].CZ2].xyz.z
           + native[native_residue[res_no].CZ3].xyz.z)/3.;
   }
  return;
 }

float aromatic_plane(int res_a, int res_b, struct vector *plane_a, struct vector *plane_b)
 {
  struct vector tmp1, tmp2, tmp3, tmp4;
  float angle;
   
  if((strcmp(native_residue[res_a].res, "PHE")==0)||(strcmp(native_residue[res_a].res, "RING")==0))
   {
    MakeVector(native[native_residue[res_a].CG].xyz, native[native_residue[res_a].CE1].xyz, &tmp1);
    MakeVector(native[native_residue[res_a].CG].xyz, native[native_residue[res_a].CE2].xyz, &tmp2);
   }
  else if(strcmp(native_residue[res_a].res, "TRP")==0)
   {
    MakeVector(native[native_residue[res_a].CG].xyz, native[native_residue[res_a].CZ2].xyz, &tmp1);
    MakeVector(native[native_residue[res_a].CG].xyz, native[native_residue[res_a].CZ3].xyz, &tmp2);
   }

  if((strcmp(native_residue[res_b].res, "PHE")==0)||(strcmp(native_residue[res_b].res, "RING")==0))
   {
    MakeVector(native[native_residue[res_b].CG].xyz, native[native_residue[res_b].CE1].xyz, &tmp3);
    MakeVector(native[native_residue[res_b].CG].xyz, native[native_residue[res_b].CE2].xyz, &tmp4);
   }
  else if(strcmp(native_residue[res_b].res, "TRP")==0)
   {
    MakeVector(native[native_residue[res_b].CG].xyz, native[native_residue[res_b].CZ2].xyz, &tmp3);
    MakeVector(native[native_residue[res_b].CG].xyz, native[native_residue[res_b].CZ3].xyz, &tmp4);
   }
  CrossProduct(tmp1, tmp2, plane_a);
  CrossProduct(tmp3, tmp4, plane_b);
  angle = Angle(*plane_a, *plane_b);
  angle *= rad2deg;
  if(angle > 90)
    angle = 180-angle;
	
  return angle;
 }
