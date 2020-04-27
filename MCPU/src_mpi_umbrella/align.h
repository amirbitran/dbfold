void SetupAlignmentStructure(void);
void SetupAlignmentPotential(void);

void SetAlignHardCore() {
  int i, j, k;
  Float temp;
  Float rad1, rad2;

  if (!USE_GO_POTENTIAL) {
    align_hard_core = (Float **) calloc(MAX_TYPES, sizeof(Float *));
    for (i=0; i<MAX_TYPES; i++)
      align_hard_core[i] = (Float *) calloc(MAX_TYPES, sizeof(Float));
    
    for (i=0; i<MAX_TYPES; i++)
      for (j=0; j<MAX_TYPES; j++) {   
	for (k=0; k<natom_type_list; ++k)
	  if (i == atom_type_list[k].type_num)
	    break;
	rad1 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
	for (k = 0; k < natom_type_list; ++k)
	  if (j == atom_type_list[k].type_num)
	    break;
	rad2 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
	
	temp = ALPHA*(rad1 + rad2);
	align_hard_core[i][j] = temp*temp;
      }
  } else {     
    align_hard_core = (Float **) calloc(struct_natoms, sizeof(Float *));
    for (i=0; i<struct_natoms; i++)
      align_hard_core[i] = (Float *) calloc(struct_natoms, sizeof(Float));
    
    for (i=0; i<struct_natoms; i++)
      for (j=0; j<struct_natoms; j++) {
	rad1 = radii[struct_native[i].atomtype];
	rad2 = radii[struct_native[j].atomtype];

	temp = ALPHA*(rad1 + rad2);
	align_hard_core[i][j] = temp*temp;
      }
  }

  return;
}

void SetAlignContactDistance() {
  Float rad1, rad2;
  Float temp;
  int i, j, k;

  if (!USE_GO_POTENTIAL) {
    align_con_dist = (struct align_cutoff **) calloc(MAX_TYPES, sizeof(struct align_cutoff *));
    for (i=0; i<MAX_TYPES; i++)
      align_con_dist[i] = (struct align_cutoff *) calloc(MAX_TYPES, sizeof(struct align_cutoff));
    
    for (i=0; i<MAX_TYPES; i++)
      for (j=0; j<MAX_TYPES; j++) {
	for (k = 0; k < natom_type_list; ++k)
	  if (i == atom_type_list[k].type_num)
	    break;
	rad1 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
	for (k = 0; k < natom_type_list; ++k)
	  if (j == atom_type_list[k].type_num)
	    break;
	rad2 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
    
	temp = LAMBDA*ALPHA*(rad1 + rad2);
	align_con_dist[i][j].b = temp*temp;
	align_con_dist[i][j].a = 0;
      }  
 
    if (MAX_TYPES==84) {
      align_con_dist[79][82].b = align_con_dist[82][79].b = 3.25*3.25;
      align_con_dist[79][83].b = align_con_dist[83][79].b = 3.25*3.25;
      align_con_dist[79][82].a = align_con_dist[82][79].a = 2.75*2.75;
      align_con_dist[79][83].a = align_con_dist[83][79].a = 2.75*2.75;
      for (j = 79; j <= 83; ++j)
	for (i = 0; i < MAX_TYPES; ++i)
	  if ((i < 79) || (i > 83))
	    align_con_dist[i][j].b =  align_con_dist[j][i].b = 0;
    }
    if (MAX_TYPES==28) {
      align_con_dist[23][26].b = align_con_dist[26][23].b = 3.25*3.25;
      align_con_dist[23][27].b = align_con_dist[27][23].b = 3.25*3.25;
      align_con_dist[23][26].a = align_con_dist[26][23].a = 2.75*2.75;
      align_con_dist[23][27].a = align_con_dist[27][23].a = 2.75*2.75;
      for (j = 23; j <= 27; ++j)
	for (i = 0; i < MAX_TYPES; ++i)
	  if ((i < 23) || (i > 27))
	    align_con_dist[i][j].b =  align_con_dist[j][i].b = 0;
    } 
  }
  else {
    align_con_dist = (struct align_cutoff **) calloc(struct_natoms, sizeof(struct align_cutoff *));
    for (i=0; i<struct_natoms; i++)
      align_con_dist[i] = (struct align_cutoff *) calloc(struct_natoms, sizeof(struct align_cutoff));
 
    for (i=0; i<struct_natoms; i++)
      for (j=0; j<struct_natoms; j++) {
	rad1 = radii[struct_native[i].atomtype];
	rad2 = radii[struct_native[j].atomtype];

	temp = LAMBDA*ALPHA*(rad1 + rad2);
	align_con_dist[i][j].b = temp*temp;
	align_con_dist[i][j].a = 0;
      }
  }
  
  return;
}

void ReadHelicityData() {
  int i;
  char line[10000], line1[10000];
  
  if((DATA = fopen(helicity_data, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", helicity_data);
    exit(1);
   }
  fgets(line, 10000, DATA);
  fgets(line, 10000, DATA);
  fgets(line1, 10000, DATA);
  
  helix = (int *) calloc(strlen(line), sizeof(int));

  for (i=0; i<strlen(line); i++) {
    if (line[i]=='H')
      helix[i] = (int) (line1[i])-'0';
    else
      helix[i]=0;
  }
  fclose(DATA);
  
  return;
}

void ReadAlignment() {
  int seq1, seq2, struct1, struct2, i;
  int nalign_seq, nalign_struct;

  map_to_seq = (int *) calloc(1000, sizeof(int));
  map_to_struct = (int *) calloc(1000, sizeof(int));
  nalign_seq = nalign_struct = 0;
  nseg = 0;
    seq1 = 0;
    seq2 = nresidues - 1;
    struct1 = 0;
    struct2 = nresidues - 1;
    for (i = seq1; i <= seq2; ++i)
      map_to_seq[nalign_seq++] = i;
    for (i = struct1; i <= struct2; ++i)
      map_to_struct[nalign_struct++] = i;
    str_segment[nseg].a = struct1;
    str_segment[nseg].b = struct2;
    seq_segment[nseg].a = seq1;
    seq_segment[nseg].b = seq2;
    ++nseg;
  nalign = nalign_seq;
}

void AlignCheckForContacts(short a, short b) {
  Float distance;

  distance = D2(struct_native[a].xyz, struct_native[b].xyz);

  if (struct_data[a][b].check_clashes && distance < align_hard_core[struct_native[a].smogtype][struct_native[b].smogtype]){
    struct_data[a][b].clashes=struct_data[b][a].clashes=1;
    struct_nclashes++;
  }
  
  if (struct_data[a][b].check_contacts)
    if ((distance <= align_con_dist[struct_native[a].smogtype][struct_native[b].smogtype].b) && (distance >= align_con_dist[struct_native[a].smogtype][struct_native[b].smogtype].a)) {
      struct_data[a][b].contacts=struct_data[b][a].contacts=1;
      struct_ncontacts++;      
    }
  
  return;
}

void AlignContacts() {
  int i, j;
  /* Re-initializes the contact matrices */
  /* Counts clashes, and native/non-native contacts */
  /* Values for ncontacts, nnon_native_contacts will not be correct */
  /*        until GetNativeContacts has been called once */  

  struct_nclashes=0;
  struct_ncontacts=0; 

  for (i=0; i<struct_natoms; i++) {
    for (j=i+1; j<struct_natoms; j++) {
      struct_data[i][j].clashes=struct_data[j][i].clashes=0;
      struct_data[i][j].contacts=struct_data[j][i].contacts=0;
      if (struct_data[i][j].check_contacts || struct_data[i][j].check_clashes)
	AlignCheckForContacts(i, j);
    }
  }
  return;
}

void SetupAlignmentStructure() {
  int i;

  struct_native = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));
  ReadNative(structure_file, struct_native, &struct_natoms);
  /* despite name, struct_natoms is an int */
 
  /* this is so unclear:
   * both functions below use struct_native struct_natoms and a bunch of globals
   */
  SetAlignHardCore();
  SetAlignContactDistance();

  struct_nresidues = 0;
  for (i=0; i<struct_natoms; i++) 
    if (!strncmp(struct_native[i].atomname, "CA", 2)) 
      struct_nresidues++;
  struct_residue = (struct residue *) calloc(struct_nresidues, sizeof(struct residue));
  GetResidueInfo(struct_native, struct_residue, struct_nresidues, struct_natoms);  
  struct_data = (struct contact_data **) calloc(struct_natoms, sizeof(struct contact_data *));
  for (i=0; i<struct_natoms; i++)
    struct_data[i] = (struct contact_data *) calloc(struct_natoms, sizeof(struct contact_data));

  GetPhiPsi(struct_native, struct_residue, struct_nresidues);
  CheckCorrelation(struct_data, struct_native, struct_residue, struct_natoms);

  seq_to_struct = (int *) calloc(nresidues, sizeof(int));
  struct_to_seq = (int *) calloc(struct_nresidues, sizeof(int));
  
  for (i=0; i <nresidues; ++i)
    seq_to_struct[i] = -1;
  for (i=0; i <struct_nresidues; ++i)
    struct_to_seq[i] = -1;
  
  for (i=0; i <nalign; ++i) {
    seq_to_struct[map_to_seq[i]]=map_to_struct[i];
    struct_to_seq[map_to_struct[i]] = map_to_seq[i];
  }

  AlignContacts();
  fprintf(STATUS, "---TEMPLATE---\n");
  fprintf(STATUS, "  # of clashes:\t\t%d\n  # of contacts:\t%d\n\n", struct_nclashes, struct_ncontacts);
    
}

void SetupAlignmentPotential() {
  int ngo_seq, ngo_str;
  int i, j, k, l, m, str_bbA[4], str_bbB[4], seq_bbA[4], seq_bbB[4];
  int **done_pairs;
  struct residue str_resA, str_resB, seq_resA, seq_resB;
  int strA, strB, seqA, seqB;
  short same_segment;

  /* setup aligned backbone potential */
  ngo_seq = 0;
  ngo_str = 0;
  for (i = 0; i < nalign; ++i)
    for (j = i+1; j < nalign; ++j) {
      str_resA = struct_residue[map_to_struct[i]];
      str_resB = struct_residue[map_to_struct[j]];
      seq_resA = native_residue[map_to_seq[i]];
      seq_resB = native_residue[map_to_seq[j]];
      str_bbA[0] = str_resA.CA; str_bbB[0] = str_resB.CA;
      str_bbA[1] = str_resA.N ; str_bbB[1] = str_resB.N ;
      str_bbA[2] = str_resA.C ; str_bbB[2] = str_resB.C ;
      str_bbA[3] = str_resA.O ; str_bbB[3] = str_resB.O ;
      seq_bbA[0] = seq_resA.CA; seq_bbB[0] = seq_resB.CA;
      seq_bbA[1] = seq_resA.N ; seq_bbB[1] = seq_resB.N ;
      seq_bbA[2] = seq_resA.C ; seq_bbB[2] = seq_resB.C ;
      seq_bbA[3] = seq_resA.O ; seq_bbB[3] = seq_resB.O ;

      /* determine whether i and j are in same aligned segment in structure */

      same_segment = 0;
      for (m = 0; m < nseg; ++m)
	if ((map_to_struct[i] >= str_segment[m].a) && (map_to_struct[i] <= str_segment[m].b) && (map_to_struct[j] >= str_segment[m].a) && (map_to_struct[j] <= str_segment[m].b))
	  same_segment = 1;

      /* generate potential: */
      /*   if aligned residues are in contact, assign attraction  */
      /*   if not, then assign repulsion only if they are in same segment */
      /*   otherwise, assign non specific energy */

      for (k = 0; k < 4; ++k)
	for (l = 0; l < 4; ++l)
	  if (struct_data[str_bbA[k]][str_bbB[l]].contacts) {
	    potential[seq_bbA[k]][seq_bbB[l]] = NATIVE_ATTRACTION;
	    potential[seq_bbB[l]][seq_bbA[k]] = NATIVE_ATTRACTION;
	    ++ngo_str;
	    ++ngo_seq;
	  }
	  else if (same_segment){
	    potential[seq_bbA[k]][seq_bbB[l]] = NON_NATIVE_REPULSION;
	    potential[seq_bbB[l]][seq_bbA[k]] = NON_NATIVE_REPULSION;
	  }
	  else {
	    potential[seq_bbA[k]][seq_bbB[l]] = NON_SPECIFIC_ENERGY;
	    potential[seq_bbB[l]][seq_bbA[k]] = NON_SPECIFIC_ENERGY;
	  }
    }
  fprintf(STATUS, "Sequence Go, backbone: %d\n", ngo_seq);
  fprintf(STATUS, "Structure Go, backbone: %d\n", ngo_str);
  
  /* setup alignment side-chain potential */

  done_pairs = (int **) calloc(nresidues, sizeof(int *));
  for (i = 0; i < nresidues; ++i)
    done_pairs[i] = (int *) calloc(nresidues, sizeof(int));

  for (i = 0; i < struct_natoms; ++i)
    for (j = i+1; j < struct_natoms; ++j) {
      strA = struct_native[i].res_num;
      strB = struct_native[j].res_num;
      seqA = struct_to_seq[strA];
      seqB = struct_to_seq[strB];
      if ((seqA != -1) && (seqB != -1))
	if (struct_native[i].is_sidechain && struct_native[j].is_sidechain) {
	  same_segment = 0;
	  for (m = 0; m < nseg; ++m)
	    if ((strA >= str_segment[m].a) && (strA <= str_segment[m].b) && (strB >= str_segment[m].a) && (strB <= str_segment[m].b))
	      same_segment = 1;
	  if (struct_data[i][j].contacts) {
	    ++ngo_str;
	    if  (!done_pairs[seqA][seqB]) {
	      done_pairs[seqA][seqB] = 1;
	      done_pairs[seqB][seqA] = 1;
	      for (k = 0; k < natoms; ++k)
		for (l = k+1; l < natoms; ++l)
		  if ((native[k].res_num == seqA) && (native[l].res_num == seqB))
		    if (native[k].is_sidechain && native[l].is_sidechain) {
		      potential[k][l] = NATIVE_ATTRACTION;
		      potential[l][k] = NATIVE_ATTRACTION;
		      ++ngo_seq;
		    }
	    }
	  }
	  else {
	    if  ((!done_pairs[seqA][seqB]) && (same_segment)) {
	      done_pairs[seqA][seqB] = 1;
	      done_pairs[seqB][seqA] = 1;
	      for (k = 0; k < natoms; ++k)
		for (l = k+1; l < natoms; ++l)
		  if ((native[k].res_num == seqA) && (native[l].res_num == seqB))
		    if (native[k].is_sidechain && native[l].is_sidechain) {
		      potential[k][l] = NON_NATIVE_REPULSION;
		      potential[l][k] = NON_NATIVE_REPULSION;
		    }
	    }
	  }
	}
    }
  

  fprintf(STATUS, "Sequence Go, all: %d\n", ngo_seq);
  fprintf(STATUS, "Structure Go, all: %d\n", ngo_str);

  for (i = 0; i < nresidues; ++i)
    free(done_pairs[i]);
  free(done_pairs);
}

