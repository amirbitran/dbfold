#define SP3_ANGLE 109.5
#define CA_CB_DISTANCE 1.54

Float Phi(struct residue, struct residue, struct atom *);
Float Psi(struct residue, struct residue, struct atom *);
Float CalculateTorsion(struct atom *, int, int, int, int, Float);
void AddCB(struct atom *, struct residue, struct atom *);
Float RadiusOfGyration(struct atom *protein, int num_atoms);
void CenterProtein(struct atom **, int);
void dRms(struct atom *, struct atom *, int, Float *, Float *);
void SSdRms(struct atom *, struct atom *, int, Float *, Float *, int, int);


void CenterProtein(struct atom **protein, int num_atoms) {
  
  int i;
  struct vector center_of_mass;
  
  Zero(&center_of_mass);
  for (i=0; i<num_atoms; i++)
    Add((*protein)[i].xyz, &center_of_mass);
  Scale(1/(Float) num_atoms, &center_of_mass);
  Inverse(&center_of_mass);
  for (i=0; i<num_atoms; i++) {
    Add(center_of_mass, &(*protein)[i].xyz);
  }

  return;

}

void dRms(struct atom *protein1, struct atom *protein2, int num_atoms, Float *sc_rms, Float *bb_rms) {
  Float scrms_tot, bbrms_tot, d1, d2;
  int sc_count=0, ca_count=0, i, j;

  scrms_tot = bbrms_tot = 0;
  for (i = 0; i < num_atoms; ++i)
    for (j = i + 1; j < num_atoms; ++j) {
      if ((protein1[i].is_sidechain) && (protein1[j].is_sidechain)) {
	d1 = sqrt(D2(protein1[i].xyz, protein1[j].xyz));
	d2 = sqrt(D2(protein2[i].xyz, protein2[j].xyz));
	scrms_tot += (d1 - d2)*(d1 - d2);                             
	++sc_count;
      }
      else if (!strcmp(protein1[i].atomname, "CA") && !strcmp(protein1[j].atomname, "CA")) {
	d1 = sqrt(D2(protein1[i].xyz, protein1[j].xyz));
	d2 = sqrt(D2(protein2[i].xyz, protein2[j].xyz));
	bbrms_tot += (d1 - d2)*(d1 - d2);
	++ca_count;
      }
    }

  *sc_rms = sqrt(scrms_tot/(Float)sc_count);
  *bb_rms = sqrt(bbrms_tot/(Float)ca_count);
}

void align_drms(struct atom *protein1, struct residue *residue1, struct atom *protein2, struct residue *residue2, int *map1, int *map2, int align_len, Float *rms) {

  Float bbrms_tot, d1, d2;
  int ca_count=0, i, j;
  
  bbrms_tot = 0;
  for (i=0; i<align_len; ++i)
    for (j=i+1; j<align_len; ++j) {
	d1 = sqrt(D2(protein1[residue1[map1[i]].CA].xyz, protein1[residue1[map1[j]].CA].xyz));
	d2 = sqrt(D2(protein2[residue2[map2[i]].CA].xyz, protein2[residue2[map2[j]].CA].xyz));
	bbrms_tot += (d1 - d2)*(d1 - d2);
	++ca_count;
    }
  *rms = sqrt(bbrms_tot/(Float)ca_count);
}

void SSdRms(struct atom *protein1, struct atom *protein2, int num_atoms, Float *sc_rms, Float *bb_rms, int first_res, int last_res) {
  Float scrms_tot, bbrms_tot, d1, d2;
  int sc_count=0, ca_count=0, i, j;

  scrms_tot = bbrms_tot = 0;
  for (i = 0; i < num_atoms; ++i)
    for (j = i + 1; j < num_atoms; ++j) {
      if ((protein1[i].res_num >= first_res) && (protein1[i].res_num <= last_res) && (protein1[j].res_num >= first_res) && (protein1[j].res_num <= last_res)) {
	if ((protein1[i].is_sidechain) && (protein1[j].is_sidechain)) {
	  d1 = sqrt(D2(protein1[i].xyz, protein1[j].xyz));
	  d2 = sqrt(D2(protein2[i].xyz, protein2[j].xyz));
	  scrms_tot += (d1 - d2)*(d1 - d2);                             
	  ++sc_count; 
	}
	else if (!strcmp(protein1[i].atomname, "CA") && !strcmp(protein1[j].atomname, "CA")) {
	  d1 = sqrt(D2(protein1[i].xyz, protein1[j].xyz));
	  d2 = sqrt(D2(protein2[i].xyz, protein2[j].xyz));
	  bbrms_tot += (d1 - d2)*(d1 - d2);
	  ++ca_count;
	}
      }
    }
  
  *sc_rms = sqrt(scrms_tot/(Float)sc_count);
  *bb_rms = sqrt(bbrms_tot/(Float)ca_count);
}

  
Float RadiusOfGyration(struct atom *protein, int num_atoms) {

  int i, num_ca_atoms;
  Float radius=0;
  struct vector center_of_mass;
  
  Zero(&center_of_mass);
  num_ca_atoms=0;
  for (i=0; i<num_atoms; i++)
    if (!strcmp(protein[i].atomname, "CA")) {
      Add(protein[i].xyz, &center_of_mass);
      num_ca_atoms++;
    }
  Scale(1/(Float) num_ca_atoms, &center_of_mass);

  for (i=0; i<num_atoms; i++)
    if (!strcmp(protein[i].atomname, "CA"))
      radius += D2(protein[i].xyz, center_of_mass);
  radius /= (Float) num_ca_atoms;
  
  return sqrt(radius);
  
}

Float Phi(struct residue current_residue, struct residue prev_residue, struct atom *protein) {

  return CalculateTorsion(protein, prev_residue.C, current_residue.N, current_residue.CA, current_residue.C, 0);


}

Float Psi(struct residue current_residue, struct residue next_residue, struct atom *protein) {

  return CalculateTorsion(protein, next_residue.N, current_residue.C, current_residue.CA, current_residue.N, 0);

}

Float CalculateTorsion(struct atom *protein, int a, int b, int c, int d, Float offset) {
  /* returns degrees */ 

 
  /*    d                               d        d                             */
  /*     \                             /           \                           */
  /*      c -- b          bc--d       cb            cb       d--cb       cb    */
  /*            \         |            |             |           |       ||    */  
  /*     	 a        a            a             a           a       da    */

  /*    torsion          -90          -130         130           90       0    */

  Float dot, norm_c, sign, angle;
  struct vector A1, B1, C1, U1, U2;
  
  Zero(&A1);
  Zero(&B1);
  Zero(&C1);
  Zero(&U1);
  Zero(&U2);
  
  MakeVector(protein[b].xyz, protein[a].xyz, &A1);  /* A1 = a - b */
  MakeVector(protein[c].xyz, protein[d].xyz, &B1);  /* B1 = d - c */
  MakeVector(protein[c].xyz, protein[b].xyz, &C1);  /* C1 = b - c */
  
  dot = Dot(A1, C1);
  norm_c = Norm(C1);
  Copy(C1, &U1);
  Scale(fabs(dot)/(norm_c*norm_c), &U1);
  Inverse(&U1);
  Add(A1, &U1);
  Normalize(&U1);  /* let u1 = A1 - C1*(<A1, C1>/<C1, C1>) ,  so U1 = u1/|u1| */
  dot = Dot(B1, C1);
  Copy(C1, &U2);
  Scale(fabs(dot)/(norm_c*norm_c), &U2);
  Add(B1, &U2);
  Normalize(&U2);  /* let u2 = B1 - C1*(<B1, C1>/<C1, C1>) ,  so U2 = u2/|u2| */

  /* now U1 and U2 are perpendicular to C1, so the angle between them is the dihedral */
  dot = Dot(U1, U2);  

  /* determine orientation */
  Zero(&A1);
  CrossProduct(U1, U2, &A1);
  if (Dot(A1, C1)>=0)
    sign = -1;
  else
    sign = 1;
  
  angle = sign*acos(Tiny(dot))*180.0/PI+offset;
  if (angle>180)
    angle-=360;
  if (angle<-180)
    angle+=360;
  
  return angle;
}


void AddCB(struct atom *protein, struct residue GLY, struct atom *CB) {

  struct vector C1, N1, U1, U2;
  
  Zero(&C1);
  Zero(&N1);
  Zero(&U1);
  Zero(&U2);
  Zero(&(CB->xyz));
  MakeVector(protein[GLY.C].xyz, protein[GLY.CA].xyz, &C1);
  Normalize(&C1);
  MakeVector(protein[GLY.N].xyz, protein[GLY.CA].xyz, &N1);
  Normalize(&N1);
  CrossProduct(C1, N1, &U1);
  Normalize(&U1);
  Scale(CA_CB_DISTANCE*sin(0.5*PI*SP3_ANGLE/PI), &U1);
  Add(C1, &U2);
  Add(N1, &U2);
  Normalize(&U2);
  Inverse(&U2);
  Scale(CA_CB_DISTANCE*cos(0.5*PI*SP3_ANGLE/PI), &U2);
  Add(U1, &(CB->xyz));
  Add(U2, &(CB->xyz));
  Add(protein[GLY.CA].xyz, &(CB->xyz));

  return;
}



