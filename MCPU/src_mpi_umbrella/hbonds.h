float HydrogenBonds();
float FoldHydrogenBonds();
void InitializeHydrogenBonding();
int MatchAtomname(char *);
long int CheckHBond(int, int);

float HydrogenBonds() {

  int i, j, res_n, res_o;
  float e;

  e=0;
  for (res_n =1; res_n<nresidues-1; res_n++)
    for (res_o =1; res_o<nresidues-1; res_o++) {
      i = native_residue[res_n].N;
      j = native_residue[res_o].O;
      if (data[i][j].check_hbond) 
       {
	data[i][j].hbond=data[j][i].hbond=CheckHBond(i, j);
	if(d_memory < HB_INNER*HB_INNER)
   	  data[i][j].closehb=data[j][i].closehb=0;
	else if((d_memory > HB_INNER*HB_INNER) && (d_memory < HB_CUTOFF*HB_CUTOFF))
   	  data[i][j].closehb=data[j][i].closehb=1;
	else if(d_memory > HB_CUTOFF*HB_CUTOFF)
   	  data[i][j].closehb=data[j][i].closehb=3;
       }
      if (data[i][j].check_hbond) 
      if (data[i][j].hbond!=NO_HBOND) {
        if(abs(native[i].res_num-native[j].res_num)>4)
	 {
	  if(d_memory < HB_INNER*HB_INNER)
	   {
	    e+=beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	   }
	  else
	    e+=HB_PENALTY*beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	 }
       else
	 {
	  if(d_memory < HB_INNER*HB_INNER)
	    e+=seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	  else
	    e+=HB_PENALTY*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	 }
      }      
    }
    e=e/1000.0*RDTHREE_CON;
  return e;

}

float FoldHydrogenBonds() {

  int i, j, k, res_n;
  int d;
  float e;

  e = 0;

  for (res_n =1; res_n<nresidues-1; res_n++) {
    i = native_residue[res_n].N;
    temp_atom = &native[i];
    temp_xyz_int = temp_atom->xyz_int;
    temp_cell3 = prev_native[i].matrix;
    temp_cell_array = temp_atom->matrix->neighbors;
    for (d=0; d<27; d++) {
      temp_cell = temp_cell_array[d];
      if (temp_cell->natoms) {
	k=0;
	O = temp_cell->natoms;
	A = temp_cell->atom_list;
	while (k<O) { /* this is O not zero */
	  j = A[k];
			  
	  if (data[i][j].check_hbond) {
	    data[i][j].hbond=data[j][i].hbond=CheckHBond(i, j);
	    if (d_memory < HB_INNER*HB_INNER)
	      data[i][j].closehb=data[j][i].closehb=0;
	    else if ((d_memory > HB_INNER*HB_INNER) && (d_memory < HB_CUTOFF*HB_CUTOFF))
	      data[i][j].closehb=data[j][i].closehb=1;
	    else if (d_memory > HB_CUTOFF*HB_CUTOFF)
	      data[i][j].closehb=data[j][i].closehb=3;
        
	    if (data[i][j].hbond!=NO_HBOND) {
	      if(abs(native[i].res_num-native[j].res_num)>4) {
		if(d_memory < HB_INNER*HB_INNER) {
		  e+=beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		} else
		  e+=HB_PENALTY*beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	      }
	      else {
		if (d_memory < HB_INNER*HB_INNER) {
		  e+=seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		} else
		  e+=HB_PENALTY*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
	      }
	    }      
	  }
	  k++;
	}
      }
    }

    if (temp_cell_array[13] != temp_cell3->neighbors[13]) {
      for (d=0; d<27; d++) {
	k = 0;
	while (k<27 && (temp_cell_array[k] != temp_cell3->neighbors[d]))
	  k++;
	if (k==27) { // if an old neighbor cell is no longer a neighbor cell...
	  temp_cell=temp_cell3->neighbors[d];
	  if (temp_cell->natoms) {
	    k=0;
	    O = temp_cell->natoms;
	    A = temp_cell->atom_list;
	    while (k<O) {
	      j = A[k];
			  
	      if (data[i][j].check_hbond) {
		data[i][j].hbond=data[j][i].hbond=CheckHBond(i, j);
		if(d_memory < HB_INNER*HB_INNER)
		  data[i][j].closehb=data[j][i].closehb=0;
		else if((d_memory > HB_INNER*HB_INNER) && (d_memory < HB_CUTOFF*HB_CUTOFF))
		  data[i][j].closehb=data[j][i].closehb=1;
		else if(d_memory > HB_CUTOFF*HB_CUTOFF)
		  data[i][j].closehb=data[j][i].closehb=3;
        
		if (data[i][j].hbond!=NO_HBOND) {
		  if(abs(native[i].res_num-native[j].res_num)>4) {
		    if(d_memory < HB_INNER*HB_INNER) {
		      e+=beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		    } else
		      e+=HB_PENALTY*beta_favor*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		  } else {
		    if (d_memory < HB_INNER*HB_INNER) {
		      e+=seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		    } else
		      e+=HB_PENALTY*seq_hb[helix_sheet][GetAminoNumber(native[i].res)][GetAminoNumber(native[j].res)]*hbond_E[data[i][j].hbond];
		  }
		}      
	      }

	      k++;
	    }
	  }
	}	
      }
    }    
  }
  e=e/1000.0*RDTHREE_CON;
  return e;

}


void InitializeHydrogenBonding() {

  int i, j, k, l, match, skip;
  int i7, i6, i5, i4, i3, i2, i1;
  int jj1, jj2, jj3, jj4, jj5, jj6, jj7;
  int D, D2, D3, D4, D5, D6, D7, A, A2, A3, A4, A5, A6, A7, A8;
  int E;
  FILE *DATA;
  FILE *fseq_hb;
  char line[250], file[250];
  char line_seq_hb[250], file_seq_hb[250];
  char res[4], res2[4], donor[4], donor2[4], donor3[4], donor4[4], donor5[4], donor6[4], donor7[4];
  char acceptor[4], acceptor2[4], acceptor3[4], acceptor4[4], acceptor5[4], acceptor6[4], acceptor7[4], acceptor8[4];
  int offset[15];
  int tmp_type, tmp_a, tmp_b;
  float tmp_val;
	
  float min_D, max_D, D_int;

  if((hbonds = (struct hydrogen_bond **) calloc(natoms, sizeof(struct hydrogen_bond *))) == NULL)
   {
    fprintf(STATUS, "ERROR: Failed to generate hbonds\n");
    exit(1);
   }
  for (i=0; i<natoms; i++)
    hbonds[i] = (struct hydrogen_bond *) calloc(natoms, sizeof(struct hydrogen_bond));

  strcpy(file_seq_hb, "../config_files/seq_dep_hb_mu_low.energy");
  fprintf(STATUS, "Opening the file: %s\n", file_seq_hb);
  if ((fseq_hb = fopen(file_seq_hb, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", file_seq_hb);
    exit(1);
   fclose(fseq_hb);
   }
  while (fgets(line_seq_hb, 250, fseq_hb)!=NULL)
   {
    sscanf(line_seq_hb, "%d %d %d %f %*s%*s", &tmp_type, &tmp_a, &tmp_b, &tmp_val);
    if (SEQ_DEP_HB == 0)
     {
      seq_hb[tmp_type][tmp_a][tmp_b] = 1.;
     }
    else
     {
      seq_hb[tmp_type][tmp_a][tmp_b] = (-1.)*tmp_val;
     }
   }
  
  if ((DATA = fopen(hydrogen_bonding_data, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", hydrogen_bonding_data);
    exit(1);
   }
  
  while (fgets(line, 250, DATA)!=NULL) {
    if (line[0]!='#') {
      sscanf(line, "%s %s %d %s %d %s %d %s %d %s %d %s %d %s %d %s %s %d %s %d %s %d %s %d %s %d %s %d %s %d %s %d %d %f %f %f %s", 
          res, donor, &offset[0], donor2, &offset[1], donor3, &offset[2], donor4, &offset[7], donor5, &offset[9], donor6, &offset[10], donor7, &offset[13], 
          res2, acceptor, &offset[3], acceptor2, &offset[4], acceptor3, &offset[5], 
          acceptor4, &offset[6], acceptor5, &offset[8], acceptor6, &offset[11], acceptor7, &offset[12], acceptor8, &offset[14], 
          &skip, &min_D, &max_D, &D_int, file);
      
     for (i=0; i<nresidues; i++) 
       for (j=0; j<nresidues; j++) {
/* *******
	 match = 0;
	 if (strcmp(res, "XXX") && strcmp(res2, "XXX")) {
	   if(!strcmp(native_residue[i].res, res) && !strcmp(native[j].res, res2)) 
	     match = 1;
	 }
	 else if (!strcmp(res, "XXX") && strcmp(res2, "XXX")) {
	   if (!strcmp(native[j].res, res2))
	     match = 1;
	 }
	 else if (strcmp(res, "XXX") && !strcmp(res2, "XXX")) {
	   if (!strcmp(native[i].res, res))
	     match = 1;
	 }   
	 else 
    *******************/
	   match = 1;

	 l=0;
	 for (k=0; k<3 && !l; k++)
	   l = (i+offset[k] < 0 || i+offset[k] >= nresidues);
         if(!l)
	   l = (i+offset[7] < 0 || i+offset[7] >= nresidues);
         if(!l)
	   l = (i+offset[9] < 0 || i+offset[9] >= nresidues);
         if(!l)
	   l = (i+offset[10] < 0 || i+offset[10] >= nresidues);
         if(!l)
	   l = (i+offset[13] < 0 || i+offset[13] >= nresidues);
	 for (k=3; k<7 && !l; k++)
	   l = (j+offset[k] < 0 || j+offset[k] >= nresidues);
         if(!l)
	   l = (j+offset[8] < 0 || j+offset[8] >= nresidues);
         if(!l)
	   l = (j+offset[11] < 0 || j+offset[11] >= nresidues);
         if(!l)
	   l = (j+offset[12] < 0 || j+offset[12] >= nresidues);
         if(!l)
	   l = (j+offset[14] < 0 || j+offset[14] >= nresidues);

	 if (!l && match && fabs(i-j)>=skip) {
	   D = native_residue[i+offset[0]].atomnumber[MatchAtomname(donor)];
	   D2 = native_residue[i+offset[1]].atomnumber[MatchAtomname(donor2)];
	   D3 = native_residue[i+offset[2]].atomnumber[MatchAtomname(donor3)];
	   D4 = native_residue[i+offset[7]].atomnumber[MatchAtomname(donor4)];
           D5 = native_residue[i+offset[9]].atomnumber[MatchAtomname(donor5)];
           D6 = native_residue[i+offset[10]].atomnumber[MatchAtomname(donor6)];
           D7 = native_residue[i+offset[13]].atomnumber[MatchAtomname(donor7)];
	   A = native_residue[j+offset[3]].atomnumber[MatchAtomname(acceptor)];
	   A2 = native_residue[j+offset[4]].atomnumber[MatchAtomname(acceptor2)];
	   A3 = native_residue[j+offset[5]].atomnumber[MatchAtomname(acceptor3)];
	   A4 = native_residue[j+offset[6]].atomnumber[MatchAtomname(acceptor4)];
           A5 = native_residue[j+offset[8]].atomnumber[MatchAtomname(acceptor5)];
           A6 = native_residue[j+offset[11]].atomnumber[MatchAtomname(acceptor6)];
           A7 = native_residue[j+offset[12]].atomnumber[MatchAtomname(acceptor7)];
           A8 = native_residue[j+offset[14]].atomnumber[MatchAtomname(acceptor8)];
				  
	   hbonds[A][D].donor = hbonds[D][A].donor = &(native[D]);
	   hbonds[A][D].donor2 = hbonds[D][A].donor2 = &(native[D2]);
	   hbonds[A][D].donor3 = hbonds[D][A].donor3 = &(native[D3]);
           hbonds[A][D].donor4 = hbonds[D][A].donor4 = &(native[D4]);
           hbonds[A][D].donor5 = hbonds[D][A].donor5 = &(native[D5]);
           hbonds[A][D].donor6 = hbonds[D][A].donor6 = &(native[D6]);
           hbonds[A][D].donor7 = hbonds[D][A].donor7 = &(native[D7]);
	   hbonds[A][D].acceptor = hbonds[D][A].acceptor = &(native[A]);
	   hbonds[A][D].acceptor2 = hbonds[D][A].acceptor2 = &(native[A2]);
	   hbonds[A][D].acceptor3 = hbonds[D][A].acceptor3 = &(native[A3]);
	   hbonds[A][D].acceptor4 = hbonds[D][A].acceptor4 = &(native[A4]);
	   hbonds[A][D].acceptor5 = hbonds[D][A].acceptor5 = &(native[A5]);
           hbonds[A][D].acceptor6 = hbonds[D][A].acceptor6 = &(native[A6]);
           hbonds[A][D].acceptor7 = hbonds[D][A].acceptor7 = &(native[A7]);
           hbonds[A][D].acceptor8 = hbonds[D][A].acceptor8 = &(native[A8]);
				 
	   hbonds[A][D].max_D2 = hbonds[D][A].max_D2 = max_D*max_D;
	   hbonds[A][D].min_D2 = hbonds[D][A].min_D2 = min_D*min_D;
	   hbonds[A][D].D_int =  hbonds[D][A].D_int = D_int;
	   data[A][D].check_hbond = data[D][A].check_hbond = 1;
//           fprintf(STATUS, "%d %d %d %d %d %d %s %d %d %d\n", D, D2, D3, A, A2, A3, donor, MatchAtomname(donor), i, j);
	 }
       }
    }

  }
  fclose(DATA);

//  n_hbond_int = (int)((max_D-min_D)/D_int)+1;
  i7 = 3*LTOR*LTOR*LTOR*LTOR*LANG*LANG;
  i6 = LTOR*LTOR*LTOR*LTOR*LANG*LANG;
  i5 = LTOR*LTOR*LTOR*LANG*LANG;
  i4 = LTOR*LTOR*LANG*LANG;
  i3 = LTOR*LANG*LANG;
  i2 = LANG*LANG;
  i1 = LANG;
  if((hbond_E = (short *) calloc(i7, sizeof(short))) == NULL)
   {
    fprintf(STATUS, "ERROR: Failed to generate hbond_E\n");
    exit(1);
   }
  NO_HBOND = -1;
  for (jj1=0; jj1<3; jj1++)
    for (jj2=0; jj2<LTOR; jj2++)
    for (jj3=0; jj3<LTOR; jj3++)
    for (jj4=0; jj4<LTOR; jj4++)
    for (jj5=0; jj5<LTOR; jj5++)
    for (jj6=0; jj6<LANG; jj6++)
    for (jj7=0; jj7<LANG; jj7++)
      hbond_E[jj1*i6+jj2*i5+jj3*i4+jj4*i3+jj5*i2+jj6*i1+jj7] = 0;

  if((DATA = fopen(file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", file);
    exit(1);
   }

  fprintf(STATUS, "Opening the file: %s\n", file);
  while (fgets(line, 250, DATA)!=NULL) {
    sscanf(line, "%d %d %d %d %d %d %d %d %*s%*s", &jj1, &jj2, &jj3, &jj4, &jj5, &jj6, &jj7, &E);
    hbond_E[jj1*i6+jj2*i5+jj3*i4+jj4*i3+jj5*i2+jj6*i1+jj7]=E;
//    fprintf(STATUS, "%2d %2d %2d %2d %7d %5d\n", jj1, jj2, jj3, jj4, jj1*i3+jj2*i2+jj3*i1+jj4, E);
//  fprintf(STATUS, "%d %d %d %d %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n", r1*20*i4+r2*i4+jj1*i3+jj2*i2+jj3*i1+jj4, 
//E, r1, r2, d1, d2, d3, d4, as, i1, jj1, jj2, jj3, jj4, r1*20*i4, r2*i4, jj1*i3, jj2*i2, jj3*i1, jj4); 
  }
  fclose(DATA);
  
//  min_E=0;
//  for (i=0; i<20*20*i4; i++)
//    if (hbond_E[i]<min_E)
//      min_E = hbond_E[i];
//  for (i=0; i<20*20*i4; i++)
//    hbond_E[i]*=(hydrogen_bond/min_E);
//  fprintf(STATUS, "hey %d %d %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n", r1, r2, d1, d2, d3, d4, i1, jj1, jj2, jj3, jj4, r1*20*i4, r2*i4, jj1*i3, jj2*i2, jj3*i1, jj4); 

  return;
  
}

int MatchAtomname(char *name) {

  int i;

  for (i=0; i<NUM_ATOMNAMES; i++)
    if (!strcmp(name, ATOMNAME[i]))
      return i;
  
  return -999;

}

long int CheckHBond(int A, int B) {
  
  struct hydrogen_bond *cur_hbond;
  struct vector V3, H;
  struct vector tmp1, tmp2, tmp3, tmp4;
  struct vector plane1, plane2;
  struct vector bisect1, bisect2;
  float d_HO;
  float d_CA_n0, d_CA_0p, d_CA_np, d_CA_00;
  float ang_Dphi, ang_Dpsi, ang_Aphi, ang_Apsi, ang_PH, ang_bH;
  float ang_PCA, ang_bCA, ang_PCAant, ang_bCAant;
  float ang_CACA;
  float min1, min2, min3;

  int res_dif;
//  int helix_sheet = 0;
  int jj1, jj2, jj3, jj4, jj5, jj6, jj7;
  long int xxx;

  cur_hbond = &(hbonds[A][B]);

  res_dif = abs((cur_hbond->donor)->res_num - (cur_hbond->acceptor)->res_num);

  if(((cur_hbond->donor)->res_num == 0) ||((cur_hbond->acceptor)->res_num == 0))
    return NO_HBOND;
  if(((cur_hbond->donor)->res_num == nresidues-1) ||((cur_hbond->acceptor)->res_num == nresidues-1))
    return NO_HBOND;

  //estimate position of hydrogen atom, H		  
  MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor2)->xyz, &H);
  MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor3)->xyz, &V3);
  Add(V3, &H);
  Normalize(&H);
  Inverse(&H);
  Add((cur_hbond->donor)->xyz, &H);

  //determine if h...O distance is less than cutoff for hbonds
  d_HO = D2(H, (cur_hbond->acceptor)->xyz);
  d_memory = d_HO;
  if (d_HO > HB_CUTOFF*HB_CUTOFF)
    return NO_HBOND;
  
  //distances between CA atoms, p = previous, 0 = same, n = next residue
  d_CA_n0 = D2((cur_hbond->donor5)->xyz, (cur_hbond->acceptor3)->xyz);
  d_CA_0p = D2((cur_hbond->donor2)->xyz, (cur_hbond->acceptor6)->xyz);
  d_CA_np = D2((cur_hbond->donor5)->xyz, (cur_hbond->acceptor6)->xyz);
  d_CA_00 = D2((cur_hbond->donor2)->xyz, (cur_hbond->acceptor3)->xyz);
  min1 = min(d_CA_n0, d_CA_np);
  min2 = min(d_CA_0p, d_CA_00);
  min3 = min(min1, min2);
//  fprintf(STATUS, "%d %f %f %f %f %f\n", res_dif, sqrt(d_CA_n0), sqrt(d_CA_0p), sqrt(d_CA_np), sqrt(d_CA_00), sqrt(min(min1, min2)));
  if(res_dif==4)
   {
    helix_sheet = 0;
    if((min1>5.8*5.8)||(min2>5.8*5.8))
      return NO_HBOND;
    if(min3>5.5*5.5)
      return NO_HBOND;
   }
  if(res_dif>4)
   {
    if((min1>6.0*6.0)||(min2>6.0*6.0))
      return NO_HBOND;
    if(min3>5.4*5.4)
      return NO_HBOND;
    //if secondary structure is predicted as helix, don't form h-bond when res_dif>4
    if((secstr[(cur_hbond->donor)->res_num]=='H') || (secstr[(cur_hbond->acceptor)->res_num]=='H'))
      return NO_HBOND;
   }

   MakeVector((cur_hbond->donor3)->xyz, (cur_hbond->donor)->xyz, &tmp1);
   MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor2)->xyz, &tmp2);
   MakeVector((cur_hbond->donor2)->xyz, (cur_hbond->donor4)->xyz, &tmp3);
   ang_Dphi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
   ang_Dphi *= rad2deg;
   MakeVector((cur_hbond->donor6)->xyz, (cur_hbond->donor5)->xyz, &tmp1);
   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor3)->xyz, &tmp2);
   MakeVector((cur_hbond->donor3)->xyz, (cur_hbond->donor)->xyz, &tmp3);
   ang_Dpsi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
   ang_Dpsi *= rad2deg;

   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor4)->xyz, &tmp1);
   MakeVector((cur_hbond->acceptor4)->xyz, (cur_hbond->acceptor6)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor6)->xyz, (cur_hbond->acceptor7)->xyz, &tmp3);
   ang_Aphi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
   ang_Aphi *= rad2deg;
   MakeVector((cur_hbond->acceptor5)->xyz, (cur_hbond->acceptor3)->xyz, &tmp1);
   MakeVector((cur_hbond->acceptor3)->xyz, (cur_hbond->acceptor2)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor4)->xyz, &tmp3);
   ang_Apsi = struct_calc_dih_ang(tmp1, tmp2, tmp3);
   ang_Apsi *= rad2deg;
   ang_Dphi += 180.;
   ang_Dpsi += 180.;
   ang_Aphi += 180.;
   ang_Apsi += 180.;
   if(res_dif==4)
   //if we are in an alpha-helical conformation, but beta-sheet is predicted in sec_str file, no h-bond
    {
     if((secstr[(cur_hbond->donor)->res_num]=='E') || (secstr[(cur_hbond->acceptor)->res_num]=='E'))
      {
       if((ang_Dphi<180.) && (ang_Dpsi<180.))
         return NO_HBOND;
       if((ang_Aphi<180.) && (ang_Apsi<180.))
         return NO_HBOND;
      }
     if((secstr[(cur_hbond->donor)->res_num]=='L') || (secstr[(cur_hbond->acceptor)->res_num]=='L'))
      {
       if((ang_Dphi<180.) && (ang_Dpsi<180.))
         return NO_HBOND;
       if((ang_Aphi<180.) && (ang_Apsi<180.))
         return NO_HBOND;
      }
    }
   if(res_dif>4)
    {
     // don't form long-range h-bonds if we occupy quadrants on the right of ramachandran plot
     if(ang_Dphi>150.)
       return NO_HBOND;
     if((ang_Dpsi>30.) && (ang_Dpsi<210.))
       return NO_HBOND;
     if(ang_Aphi>150.)
       return NO_HBOND;
     if((ang_Apsi>30.) && (ang_Apsi<210.))
       return NO_HBOND;
     // don't form long-range h-bonds if L conformation is predicted
     if(secstr[(cur_hbond->donor)->res_num]=='L')
       return NO_HBOND;
     if(secstr[(cur_hbond->acceptor)->res_num]=='L')
       return NO_HBOND;
    }
  
   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor7)->xyz, &tmp1);
   MakeVector((cur_hbond->acceptor8)->xyz, (cur_hbond->acceptor6)->xyz, &tmp3);
   Normalize(&tmp1);
   Normalize(&tmp3);
   ang_CACA = Angle(tmp1, tmp3);
   ang_CACA *= rad2deg;
   if(res_dif>4)
    {
     if(ang_CACA < 90)
       helix_sheet = 1;
     else
       helix_sheet = 2;
    }
	  
   MakeVector((cur_hbond->donor2)->xyz, (cur_hbond->donor)->xyz, &tmp1);
   MakeVector((cur_hbond->donor2)->xyz, (cur_hbond->donor4)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor3)->xyz, (cur_hbond->acceptor5)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor3)->xyz, (cur_hbond->acceptor2)->xyz, &tmp4);
   CrossProduct(tmp1, tmp2, &plane1);
   CrossProduct(tmp3, tmp4, &plane2);
   ang_PCA = Angle(plane1, plane2);
   ang_PCA *= rad2deg;

   MakeVector((cur_hbond->donor2)->xyz, (cur_hbond->donor)->xyz, &tmp1);
   MakeVector((cur_hbond->donor2)->xyz, (cur_hbond->donor4)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor3)->xyz, (cur_hbond->acceptor5)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor3)->xyz, (cur_hbond->acceptor2)->xyz, &tmp4);
   Normalize(&tmp1);
   Normalize(&tmp2);
   Normalize(&tmp3);
   Normalize(&tmp4);
   bisect(tmp1, tmp2, &bisect1);
   bisect(tmp3, tmp4, &bisect2);
   ang_bCA = Angle(bisect1, bisect2);
   ang_bCA *= rad2deg;
	  
   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor6)->xyz, &tmp1);
   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor3)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor6)->xyz, (cur_hbond->acceptor4)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor6)->xyz, (cur_hbond->acceptor7)->xyz, &tmp4);
   CrossProduct(tmp1, tmp2, &plane1);
   CrossProduct(tmp3, tmp4, &plane2);
   ang_PCAant = Angle(plane1, plane2);
   ang_PCAant *= rad2deg;

   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor6)->xyz, &tmp1);
   MakeVector((cur_hbond->donor5)->xyz, (cur_hbond->donor3)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor6)->xyz, (cur_hbond->acceptor4)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor6)->xyz, (cur_hbond->acceptor7)->xyz, &tmp4);
   Normalize(&tmp1);
   Normalize(&tmp2);
   Normalize(&tmp3);
   Normalize(&tmp4);
   bisect(tmp1, tmp2, &bisect1);
   bisect(tmp3, tmp4, &bisect2);
   ang_bCAant = Angle(bisect1, bisect2);
   ang_bCAant *= rad2deg;
   ang_Dphi = ang_PCA;
   ang_Dpsi = ang_bCA;
   ang_Aphi = ang_PCAant;
   ang_Apsi = ang_bCAant;
	  
   MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor3)->xyz, &tmp1);
   MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor2)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor3)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor4)->xyz, &tmp4);
   CrossProduct(tmp1, tmp2, &plane1);
   CrossProduct(tmp3, tmp4, &plane2);
   ang_PH = Angle(plane1, plane2);
   ang_PH *= rad2deg;

   MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor3)->xyz, &tmp1);
   MakeVector((cur_hbond->donor)->xyz, (cur_hbond->donor2)->xyz, &tmp2);
   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor3)->xyz, &tmp3);
   MakeVector((cur_hbond->acceptor2)->xyz, (cur_hbond->acceptor4)->xyz, &tmp4);
   Normalize(&tmp1);
   Normalize(&tmp2);
   Normalize(&tmp3);
   Normalize(&tmp4);
   bisect(tmp1, tmp2, &bisect1);
   bisect(tmp3, tmp4, &bisect2);
   ang_bH = Angle(bisect1, bisect2);
   ang_bH *= rad2deg;
//   fprintf(STATUS, "%d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", helix_sheet, ang_Dphi, ang_Dpsi, ang_Aphi, ang_Apsi, ang_PH, ang_bH);
//   fprintf(STATUS, "%d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", helix_sheet, ang_Dphi-180, ang_Dpsi-180, ang_Aphi-180, ang_Apsi-180, ang_PH, ang_bH);
	  
   jj1 = helix_sheet;
   jj2 = (int) (ang_Dphi/ang3_int);
   jj3 = (int) (ang_Dpsi/ang3_int);
   jj4 = (int) (ang_Aphi/ang3_int);
   jj5 = (int) (ang_Apsi/ang3_int);
   jj6 = (int) (ang_PH/ang_int);
   jj7 = (int) (ang_bH/ang_int);

  xxx = jj7 + jj6*LANG + jj5*LANG*LANG + jj4*LTOR*LANG*LANG + jj3*LTOR*LTOR*LANG*LANG + jj2*LTOR*LTOR*LTOR*LANG*LANG + jj1*LTOR*LTOR*LTOR*LTOR*LANG*LANG;
//  fprintf(STATUS, "%ld\n", xxx);
  return jj7 + jj6*LANG + jj5*LANG*LANG + jj4*LTOR*LANG*LANG + jj3*LTOR*LTOR*LANG*LANG + jj2*LTOR*LTOR*LTOR*LANG*LANG + jj1*LTOR*LTOR*LTOR*LTOR*LANG*LANG;
}
