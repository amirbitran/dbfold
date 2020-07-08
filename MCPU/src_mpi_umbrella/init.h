void SetProgramOptions();
void InitializeProtein();
void DetermineTriplets();
void SetRadii();
void SetHardCore();
void SetContactDistance();
void ReadNative();
void GetResidueInfo(struct atom *, struct residue *, int, int);
void GetPhiPsi(struct atom *, struct residue *, int);
void GetChi();
void InitializeData();
void ReadSidechainTorsionData();
void InitializeSidechainRotationData();
void InitializeBackboneRotationData();
void CheckCorrelation(struct contact_data **, struct atom *, struct residue *, int);
void ReadPotential();
void ReadDistPotential();
int SkipSelf(int, int, struct atom *, struct residue *);
int SkipNeighbors(int, int, struct atom *, struct residue *);
void TurnOffNativeClashes(int);
void ReadAlignment();
void SetupAlignmentStructure(); 
void SetupAlignmentPotential(); 
void ReadAvgChis(void);
void ReadHelicityData();
void InitializeHydrogenBonding();
int MatchAtomname(char *);

/*================================================*/
/*           main initialization routine          */
/*================================================*/

void ReadTypesFile(void) {
  int type_num, i;
  char res_name[4], atom_name[4];
  FILE *type_file;

  if ((type_file = fopen(atom_type_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", atom_type_file);
    exit(1);
   }
  if (type_file == NULL) {
    fprintf(STATUS, "Error: Type file non-existant\n");
    exit(-1);
  }
  i = -1;
  while (fscanf(type_file, "%s %s %d", atom_name, res_name, &type_num)!= EOF) {
    ++i;
    strcpy(atom_type_list[i].atom_name, atom_name);
    strcpy(atom_type_list[i].res_name, res_name);
    atom_type_list[i].type_num = type_num;
    if (!strcmp(res_name, "XXX")) {
      if (!strcmp(atom_name, "N"))
	bb_N_type = type_num;
      if (!strcmp(atom_name, "O"))
	bb_O_type = type_num;
      if (!strcmp(atom_name, "OXT"))
	bb_OXT_type = type_num;
    }
  }
  natom_type_list = i+1;
  MAX_TYPES = atom_type_list[i].type_num + 1;
  fclose(type_file);
}

void SetupMuPotential (void) {
  int i, j, k, l;
  Float n_nat, avg;
  int **nat_cons, **nnat_cons;

  nat_cons = (int **) calloc(MAX_TYPES, sizeof(int *));
  nnat_cons = (int **) calloc(MAX_TYPES, sizeof(int *));
  for (i = 0; i < MAX_TYPES; ++i) {
    nat_cons[i] = (int *) calloc(MAX_TYPES, sizeof(int));
    nnat_cons[i] = (int *) calloc(MAX_TYPES, sizeof(int));
  }

  for (k = 0; k < natoms; ++k) 
    for (l = k+1; l < natoms; ++l) {
      if (data[k][l].contacts) {
	if (native[k].smogtype <= native[l].smogtype)
	  nat_cons[native[k].smogtype][native[l].smogtype]++;
	else
	  nat_cons[native[l].smogtype][native[k].smogtype]++;
      }
      else {
	if (native[k].smogtype <= native[l].smogtype)
	  nnat_cons[native[k].smogtype][native[l].smogtype]++;
	else
	  nnat_cons[native[l].smogtype][native[k].smogtype]++;
      }
    }

  avg = 0;
  n_nat = 0;
  for (i = 0; i < MAX_TYPES; ++i)
    for (j = i; j < MAX_TYPES; ++j) {
      if (nat_cons[i][j] != 0) {
	avg += nnat_cons[i][j]/nat_cons[i][j];
	n_nat = n_nat + 1;
      }
    }
  mu = (avg/(float) n_nat)/(1 + avg/(float) n_nat);
  fprintf(STATUS, "mu: %f\n", mu);

  for (i = 0; i < MAX_TYPES; ++i)
    for (j = i; j < MAX_TYPES; ++j)
      if ((nat_cons[i][j] !=0) || (nnat_cons[i][j]!=0)) {
	if (mu == 1) {
	  if (nat_cons[i][j] != 0)
	    potential[i][j] = -1;
	  else
	    potential[i][j] = 1;
	}
	else if (mu == 0) {
	  if (nnat_cons[i][j] != 0)
	    potential[i][j] = 1;
	  else
	    potential[i][j] = -1;
	}
	else   /*AB: Here is the mu potential, ex. eq. (2) in Yang..Shakhnovich et. al PNAS 2007*/
	  potential[i][j] = ((1-mu)*nnat_cons[i][j] - mu*nat_cons[i][j])/(mu*nat_cons[i][j] + (1-mu)*nnat_cons[i][j]); 
	potential[j][i] = potential[i][j];
      }

  for (i = 0; i < MAX_TYPES; ++i) {
    free(nat_cons[i]);
    free(nnat_cons[i]);
  }
  free(nat_cons);
  free(nnat_cons);
  
}


void InitializeProtein() {
  int i;
  /* Reset variables */
  /**See atom.h to learn about how atom structures are set up/
  fprintf(STATUS, "---MODEL---\n");
  fprintf(STATUS, "  file:\t\t%s   \n", native_file);  /*Set up arrays for native file, minimized file, and file with min RMSD*/
  native = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));  /*Allocates memory for an array of atom_structures, each of which corresponds to data for one atom*/
  native_Emin = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));
  native_RMSDmin = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));
  prev_native = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));
  orig_native = (struct atom *) calloc(MAX_ATOMS, sizeof(struct atom));
  natoms = 0;

  fprintf(STATUS, "Initialized molecular data structures.\n");
 
  /* Initialize static data */

  amino_acids = (struct amino *) calloc(20, sizeof(struct amino));
  fprintf(STATUS, "Made amino_acids structure\n");

  ReadTypesFile(); /*Left off here*/
  fprintf(STATUS, "Read Types File\n");
  ReadHelicityData();
  fprintf(STATUS, "Read Helicity Data\n");
  
  /* Read data from native_file into native */
  ReadNative(native_file, native, &natoms);  //AB: This reads native_file into the structure native
  fprintf(STATUS, "Read Native\n");

  /* Count nresidues */
  nresidues = 0;
  for (i=0; i<natoms; i++)
    if (!strncmp(native[i].atomname, "CA", 2)) nresidues++;

  if (nresidues!=(native[natoms-1].res_num+1))
    fprintf(STATUS, "FILE: %s   MISSING RESIDUES!!!\n", native_file);

  fprintf(STATUS, "  pdb length:\t\t%d\n  # of CA's:\t\t%d\n\n", nresidues, native[natoms-1].res_num+1);
  
  buf_in = (float *) calloc(3*natoms, sizeof(float));
  buf_out = (float *) calloc(3*natoms, sizeof(float));

  ReadAlignment();  
  fprintf(STATUS, "Read Alignment\n");
  get_template();  
  fprintf(STATUS, "Read Template Information\n");

  if (USE_GO_POTENTIAL) MAX_TYPES = natoms;

  SetRadii();
  SetHardCore();
  SetContactDistance();  

  /* Allocate potential-related structures */
  
  potential = (Float **) calloc(MAX_TYPES, sizeof(Float *));
  for (i=0; i<MAX_TYPES; i++) {
    potential[i] = (Float *) calloc(MAX_TYPES, sizeof(Float));
  }
  
  CenterProtein(&native, natoms);
  for (i=0;i<natoms; i++) 
    FindLatticeCoordinates(&native[i]); /* this also computes the integer value of the coordinates */
    
  /* Get residue info */

  native_residue = (struct residue *) calloc(nresidues, sizeof(struct residue)); //I believe native_residue stores the indices for the different atoms (where they are located within native), and other info
  cur_rotamers = (int *) calloc(nresidues, sizeof(int));
  GetResidueInfo(native, native_residue, nresidues, natoms);   //This adds the necessary info into native_residue
 
  /* Allocate memory for data structures */
  
  InitializeData(); 

  /* Determine phi-psi angles */
  
  GetPhiPsi(native, native_residue, nresidues);

  /* Set up correlation matrix */
  
  CheckCorrelation(data, native, native_residue, natoms);

  /* Set up side-chain rotation structure */

  if (USE_SIDECHAINS) {
    /* Initialize sidechain rotation data */
    ReadSidechainTorsionData();
    InitializeSidechainRotationData();
    /* initialize sidechain torsions */ 
    GetChi();
    if (USE_ROTAMERS) {
      ReadAvgChis();
      /* this will generate 1 side-chain move in the absence of backbone moves: */
      SIDECHAIN_MOVES = 1; 
    }
  }
  else
    SIDECHAIN_MOVES = 0;

  /* initialize moved atoms data structures */

  ab = (struct pair *) calloc(natoms*natoms, sizeof(struct pair));
  cd = (struct pair *) calloc(natoms*natoms, sizeof(struct pair));


  initialize_torsion();
  initialize_sct();
  initialize_aromatic();
  initialize_secstr();
  read_cluster();
  if (weight_hbond){
    InitializeHydrogenBonding();
    hbond_pair = (struct pair *) calloc(natoms*natoms, sizeof(struct pair));
  }

  /* Determine residue triplets */

  DetermineTriplets();
  
  /* Initialize backbone rotation data */
  InitializeBackboneRotationData();

  /* Get contacts */

  Contacts();

  /* Setup potential and related structures */

  SetupAlignmentStructure();

  if (USE_GO_POTENTIAL){
    mu = 1;
    SetupAlignmentPotential();
  }
  else if (READ_POTENTIAL)
    ReadPotential();
  else
    SetupMuPotential();
      
//  if (!USE_GO_POTENTIAL) {
//    potential[bb_O_type][bb_N_type] = hydrogen_bond;
//    potential[bb_N_type][bb_O_type] = hydrogen_bond;
//    potential[bb_OXT_type][bb_N_type] = hydrogen_bond;
//    potential[bb_N_type][bb_OXT_type] = hydrogen_bond;
//  }
  
  TypeContacts();
  native_E = FullAtomEnergy();

//  free(amino_acids);
  
  return;

}

void TurnOffNativeClashes(int ReportBack) {
  int i, j;

  for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++)
      if (data[i][j].clashes) {
	if (ReportBack)
	  fprintf(STATUS, "native clash\t%d %s - %s\t%d %s - %s \t %.3f %.3f\n", i, native[i].res, native[i].atomname, j, native[j].res, native[j].atomname, sqrt(D2(native[i].xyz, native[j].xyz)), sqrt(hard_core[native[i].smogtype][native[j].smogtype])/100);
	data[i][j].clashes = data[j][i].clashes = 0;
	nclashes--;
      }
}

/*================================================*/
/*     routines for initializing static data      */
/*================================================*/

void SetHardCore() {
  int i, j, k;
  Float temp;
  Float rad1, rad2;

  hard_core = (long int **) calloc(MAX_TYPES, sizeof(long int *));
  for (i=0; i<MAX_TYPES; i++)
    hard_core[i] = (long int *) calloc(MAX_TYPES, sizeof(long int));

  /* smogtype is always used for atom sizes */
  for (i=0; i<MAX_TYPES; i++)
    for (j=0; j<MAX_TYPES; j++) {      
      if (!USE_GO_POTENTIAL) {
	for (k = 0; k < natom_type_list; ++k)
	  if (i == atom_type_list[k].type_num)
	    break;
	rad1 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
	for (k = 0; k < natom_type_list; ++k)
	  if (j == atom_type_list[k].type_num)
	    break;
	rad2 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
      }
      else {
	rad1 = radii[native[i].atomtype];
	rad2 = radii[native[j].atomtype];
      }
      temp = ALPHA*(rad1 + rad2);
      hard_core[i][j] = (long int) (temp*temp*INT_PRECISION*INT_PRECISION);
    }

  return;

}

void SetContactDistance() {
  Float rad1, rad2;
  Float temp;
  Float temp00;
  int i, j, k;

  contact_distance = (struct cutoff **) calloc(MAX_TYPES, sizeof(struct cutoff *));
  for (i=0; i<MAX_TYPES; i++)
    contact_distance[i] = (struct cutoff *) calloc(MAX_TYPES, sizeof(struct cutoff));

  for (i=0; i<MAX_TYPES; i++)
    for (j=0; j<MAX_TYPES; j++) {
      if (!USE_GO_POTENTIAL) {
	for (k = 0; k < natom_type_list; ++k)
	  if (i == atom_type_list[k].type_num)
	    break;
	rad1 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
	for (k = 0; k < natom_type_list; ++k)
	  if (j == atom_type_list[k].type_num)
	    break;
	rad2 = radii[TypeAtom(atom_type_list[k].atom_name, atom_type_list[k].res_name)];
      }
      else {
	rad1 = radii[native[i].atomtype];
	rad2 = radii[native[j].atomtype];
      }
      temp = LAMBDA*ALPHA*(rad1 + rad2);
      contact_distance[i][j].b = (long int) (temp*temp*INT_PRECISION*INT_PRECISION);
      contact_distance[i][j].a = 0;
      beta = 0.;
      temp00 = beta*(rad1 + rad2);
        contact_distance[i][j].a = (long int) (temp00*temp00*INT_PRECISION*INT_PRECISION);
			
    }

  /* sets up proper distances for h-bonding */
  if (!USE_GO_POTENTIAL) {
//      contact_distance[bb_O_type][bb_N_type].b = 3.25*3.25*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_N_type][bb_O_type].b = 3.25*3.25*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_OXT_type][bb_N_type].b = 3.25*3.25*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_N_type][bb_OXT_type].b = 3.25*3.25*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_O_type][bb_N_type].a = 2.75*2.75*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_N_type][bb_O_type].a = 2.75*2.75*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_OXT_type][bb_N_type].a = 2.75*2.75*INT_PRECISION*INT_PRECISION;
//      contact_distance[bb_N_type][bb_OXT_type].a = 2.75*2.75*INT_PRECISION*INT_PRECISION;
      /* turns off all backbone-sidechain contacts */
      for (j = bb_N_type; j <= bb_OXT_type; ++j)
	for (i = 0; i < MAX_TYPES; ++i)
	  if ((i < bb_N_type) || (i > bb_OXT_type))
	   {
	    contact_distance[i][j].a =  contact_distance[j][i].a = 0;
	    contact_distance[i][j].b =  contact_distance[j][i].b = 0;
	   }
  }
  return;
}

void SetRadii() {
  
  radii = (Float *) calloc(13, sizeof(Float));

  /* carbons */
  radii[0] = 1.61;
  radii[1] = 1.76;
  radii[2] = radii[3] = radii[4] = 1.88; 
  /* nitrogens */
  radii[5] = radii[6] = radii[7] = radii[8] = 1.64;
  /* oxygens */
  radii[9] = 1.42;
  radii[10] = 1.46;
  /* sulfurs */
  radii[11] = radii[12] = 1.77;
  
  return;

}


/*=============================================*/
/*   routines for initializing and reading     */
/*               the native pdb                */
/*=============================================*/

void ReadNative(char *file_name, struct atom *protein, int *Natoms) {
  FILE *the_file;
  char line[250];

  *Natoms = 0;
  if ((the_file = fopen(file_name, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", file_name);
    exit(1);
   }
  while (fgets(line, 100, the_file) != NULL) 
    ParsePDBLine(line, protein, Natoms);
  fclose(the_file);
  
  return;  
}

void GetResidueInfo(struct atom *Chain, struct residue *Residue, int Nres, int Natoms) {
  int i, j;
  
  for (i=0; i<Natoms; i++) 
    for (j=0; j<40; j++) 
      Residue[Chain[i].res_num].atomnumber[j] = -1;
  
  for (i=0; i<Natoms; i++) {
    j = MatchAtomname(Chain[i].atomname);
    Residue[Chain[i].res_num].atomnumber[j] = i;
    if (!strncmp(Chain[i].atomname, "CA", 2)) {
      Residue[Chain[i].res_num].CA = i;
      strcpy(Residue[Chain[i].res_num].res, Chain[i].res);
      Residue[Chain[i].res_num].amino_num = GetAminoNumber(Chain[i].res);
      Residue[Chain[i].res_num].psi=0;
      Residue[Chain[i].res_num].phi=0;
      Residue[Chain[i].res_num].chi[0]=0;
      Residue[Chain[i].res_num].chi[1]=0;
      Residue[Chain[i].res_num].chi[2]=0;
      Residue[Chain[i].res_num].chi[3]=0;
      Residue[Chain[i].res_num].is_core=Chain[i].is_core;
      Residue[Chain[i].res_num].is_designed=Chain[i].is_designed;
    }
    else if (!strcmp(Chain[i].atomname, "N")) 
      Residue[Chain[i].res_num].N = i;
    else if (!strcmp(Chain[i].atomname, "C")) 
      Residue[Chain[i].res_num].C = i;
    else if (!strcmp(Chain[i].atomname, "O")) 
      Residue[Chain[i].res_num].O = i;
    else if (!strcmp(Chain[i].atomname, "CB")) 
      Residue[Chain[i].res_num].CB = i;
    else if (!strcmp(Chain[i].atomname, "CG"))
      Residue[Chain[i].res_num].CG = i;
    else if (!strcmp(Chain[i].atomname, "CE1"))
      Residue[Chain[i].res_num].CE1 = i;
    else if (!strcmp(Chain[i].atomname, "CE2"))
      Residue[Chain[i].res_num].CE2 = i;
    else if (!strcmp(Chain[i].atomname, "CZ2"))
      Residue[Chain[i].res_num].CZ2 = i;
    else if (!strcmp(Chain[i].atomname, "CZ3"))
      Residue[Chain[i].res_num].CZ3 = i;
  }
  
  /* Handle CB atoms for glycine */
  
  for (i=0; i<Nres; i++) {
    if (!strcmp(Residue[i].res, "GLY")) {
      Residue[i].CB = -999;
    }
    /*  AddCB(native, native_residue[i]); */
  }
  
  return;
  
}

  
void GetPhiPsi(struct atom *Chain, struct residue *Residue, int Nres) { /* Phi/Psi angles are stored in radians */
  int i;
  
  for (i=0; i<Nres; i++) {
    if (i!=0)
      Residue[i].phi = PI/180.0*Phi(Residue[i], Residue[i-1], Chain);
    else
      Residue[i].phi = -999;
    if (i!=Nres-1)
      Residue[i].psi = PI/180.0*Psi(Residue[i], Residue[i+1], Chain);
    else
      Residue[i].psi = -999;
  }
  
  return;

}

void GetChi() { /* Chi angles are stored in radians */
  int i, j;
  /* this routine also resets the native chis to the current values */
  
  for (i=0; i<nresidues; i++) 
    for (j=0; j<amino_acids[native_residue[i].amino_num].ntorsions; j++) {
      native_residue[i].chi[j] = PI/180.0*CalculateTorsion(native, sidechain_torsion[i][j][0], sidechain_torsion[i][j][1], sidechain_torsion[i][j][2], sidechain_torsion[i][j][3], 0); 
      native_residue[i].native_chi[j] = native_residue[i].chi[j];
      native_residue[i].tmpchi[j] = native_residue[i].chi[j];
    }
  
  return;
}

void ReadAvgChis (void) {
  int i;
  char name[4], line[200];
  Float X, Y, Z, W;
  float value;
  float sX, sY, sZ, sW;

  rotamer_angles = calloc(nresidues, sizeof(struct angles));
  
  if ((DATA = fopen(rotamer_data_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", rotamer_data_file);
    exit(1);
   }
  while (fgets(line, 150, DATA)!=NULL) {
    sscanf(line, "%s %*d%*d%*d%*d %*d%*d %f%*f%*f%*f  %f%f %f%f %f%f %f%f", 
           name, &value, &X, &sX, &Y, &sY, &Z, &sZ, &W, &sW);
    for (i = 0; i < nresidues; ++i) {
      if (strcmp(name, native_residue[i].res)==0)
       {
	rotamer_angles[i].chis[no_chi_list[GetAminoNumber(name)]][0] = X;
	rotamer_angles[i].chis[no_chi_list[GetAminoNumber(name)]][1] = Y;
	rotamer_angles[i].chis[no_chi_list[GetAminoNumber(name)]][2] = Z;
	rotamer_angles[i].chis[no_chi_list[GetAminoNumber(name)]][3] = W;
       }	
    }
   deviation_ang[GetAminoNumber(name)][no_chi_list[GetAminoNumber(name)]][0] = sX;
   deviation_ang[GetAminoNumber(name)][no_chi_list[GetAminoNumber(name)]][1] = sY;
   deviation_ang[GetAminoNumber(name)][no_chi_list[GetAminoNumber(name)]][2] = sZ;
   deviation_ang[GetAminoNumber(name)][no_chi_list[GetAminoNumber(name)]][3] = sW;
   prob_ang[GetAminoNumber(name)][no_chi_list[GetAminoNumber(name)]] = value;
   no_chi_list[GetAminoNumber(name)]++;
  }   
  fclose(DATA); 
// for (i = 0; i < nresidues; ++i) 
//   fprintf(STATUS, "%3d %s %2d %3d %8.3f %8.3f %8.3f\n", 
//     i, native_residue[i].res, native_residue[i].amino_num, 
//     no_chi_list[native_residue[i].amino_num], rotamer_angles[i].chis[2][0], 
//     deviation_ang[native_residue[i].amino_num][2][0], prob_ang[native_residue[i].amino_num][2]);
// exit(0);
}


void ReadSidechainTorsionData() {
  int temp_x, i, j, k;
  short temp;
  char  AA[4], BB[4], CC[4], DD[4], symbol[2], name[4], line[250];

  if ((DATA = fopen(amino_data_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", amino_data_file);
    exit(1);
   }
  
  temp =0;
  while (fgets(line, 150, DATA)!=NULL) {
    if (!strncmp(line, "*", 1))
      temp++;
    else if (strncmp(line, "!!", 2)) {
      if (temp == 0)
	sscanf(line, "%*s %*s");
      else if (temp == 1) {
	sscanf(line, "%d %s %*d %d %*s %s", &i, name, &j, symbol);
	strcpy(amino_acids[i].name, name);
	strcpy(amino_acids[i].symbol, symbol);
	amino_acids[i].ntorsions = j;
	if (!strcmp(name, "PRO"))
	  amino_acids[i].nrotamers = 2;
	else if (!strcmp(name, "TYR") || !strcmp(name, "HIS") || !strcmp(name, "PHE"))
	  amino_acids[i].nrotamers = 6;
	else if (!strcmp(name, "GLY"))
	  amino_acids[i].nrotamers = 0;
	else
	  amino_acids[i].nrotamers = (int) three[amino_acids[i].ntorsions];
      }
      else if (temp ==2) {
	sscanf(line, "%s %d %s %s %s %s %*f", name, &temp_x, AA, BB, CC, DD);
	for (i=0; i<20; i++) 
	  if (!strcmp(name, amino_acids[i].name))
	    break;
	strcpy(amino_acids[i].torsion[temp_x][0], AA);
	strcpy(amino_acids[i].torsion[temp_x][1], BB);
	strcpy(amino_acids[i].torsion[temp_x][2], CC);
	strcpy(amino_acids[i].torsion[temp_x][3], DD);
	
      }
      else if (temp == 3) {
	strncpy(name, line, 3);
	for (i=0; i<20; i++)
	  if (!strncmp(name, amino_acids[i].name, 3))
	    break;
	strtok(line, " \t");
	j = atoi(strtok(NULL, " \t\n"));
	amino_acids[i].rotate_natoms[j] = atoi(strtok(NULL, " \t\n"));
	for (k=0; k<amino_acids[i].rotate_natoms[j]; k++) 
	  strcpy(amino_acids[i].rotate_atom[j][k], strtok(NULL, " \t\n"));
      }
    }
  }
  fclose(DATA);

  return;
}

/*=============================================================*/
/*                Read full atom potential                     */
/*=============================================================*/

void ReadPotential(){
  FILE *pot_file;
  int i, j;
  float val;
  /* read max types potential */
  
  if ((pot_file = fopen(potential_file, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", potential_file);
    exit(1);
   }
  while (fscanf(pot_file, "%d %d %f", &i, &j, &val)!= EOF) {
    potential[i][j] = potential[j][i] = val;
  }
    
  fclose(pot_file);
}

/*=============================================================*/
/*                   allocate data structures                  */
/*=============================================================*/

void InitializeData() {
  int i;
#if DEBUG
  debug_contacts = (unsigned char **) calloc(natoms, sizeof(unsigned char *));
  for (i=0; i<natoms; i++)
    debug_contacts[i] = (unsigned char *) calloc(natoms, sizeof(unsigned char));
  debug_dcontacts = (unsigned char **) calloc(natoms, sizeof(unsigned char *));
  for (i=0; i<natoms; i++)
    debug_dcontacts[i] = (unsigned char *) calloc(natoms, sizeof(unsigned char));
  debug_clashes = (unsigned char **) calloc(natoms, sizeof(unsigned char *));
  for (i=0; i<natoms; i++)
    debug_clashes[i] = (unsigned char *) calloc(natoms, sizeof(unsigned char));
#endif

  data = (struct contact_data **) calloc(natoms, sizeof(struct contact_data *));
  for (i=0; i<natoms; i++)
    data[i] = (struct contact_data *) calloc(natoms, sizeof(struct contact_data));
  
  type_contacts = (short **) calloc(MAX_TYPES, sizeof(short *));
  for (i=0; i<MAX_TYPES; i++)
    type_contacts[i] = (short *) calloc(MAX_TYPES, sizeof(short));

  is_rotated = (unsigned char *) calloc(natoms, sizeof(unsigned char));
  for (i=0; i<natoms; i++)
    is_rotated[i]=0;

  return;
  
}

/*=============================================================*/
/*           routines for initializing the move set            */
/*=============================================================*/


void DetermineTriplets() {
  int i, j, k;

  residue_triplets = (struct triplet *) calloc((nresidues-5)*81, sizeof(struct triplet));
  /* this is way more memory than needed, given that the exact number of triplets is 16n - 55 */
  total_triplets = 0;
  for (i = 0; i <nresidues; i++) {
    residue_triplets[total_triplets].a = i;
    residue_triplets[total_triplets].b = -1;
    residue_triplets[total_triplets++].c = -1;
  }
  TOTAL_SINGLE_LOOP_MOVES = total_triplets;
  
  for (i = 0; i < nresidues; i++)
    for (j = i+1; j < nresidues; j++)
      if (j<i+6) {
	residue_triplets[total_triplets].a = i;
	residue_triplets[total_triplets].b = j;
	residue_triplets[total_triplets++].c = -1;
      }
  TOTAL_DOUBLE_LOOP_MOVES = total_triplets-TOTAL_SINGLE_LOOP_MOVES;
  
  for (i = 0; i < nresidues; i++)
    for (j = i+1; j < nresidues; j++)
      for (k = j+1; k < nresidues; k++)
	if (k<i+6 && j<i+6) { 
	  residue_triplets[total_triplets].a = i;
	  residue_triplets[total_triplets].b = j;
	  residue_triplets[total_triplets++].c = k;
	}
  TOTAL_TRIPLE_LOOP_MOVES = total_triplets-TOTAL_SINGLE_LOOP_MOVES-TOTAL_DOUBLE_LOOP_MOVES;

  return;
  
}

void InitializeBackboneRotationData() {
  int i, j, k;
  
  /* rotate_atom[0=psi, 1=phi][which residue][list of atoms] */

  yang_rotated_atoms = (short *) calloc(natoms, sizeof(short));
  yang_not_rotated = (char *) calloc(natoms, sizeof(char));
  
  rotate_natoms = (short **) calloc(2, sizeof(short *));
  rotate_atom = (short ***) calloc(2, sizeof(short **));
  not_rotated = (char ***) calloc(2, sizeof(char **));
  for (i=0; i<2; i++) {
    rotate_atom[i] = (short **) calloc(nresidues, sizeof(short *));
    not_rotated[i] = (char **) calloc(nresidues, sizeof(char *));
    rotate_natoms[i] = (short *) calloc(nresidues, sizeof(short));
    for (j=0; j<nresidues; j++) {
      rotate_atom[i][j] = (short *) calloc(natoms, sizeof(short));
      not_rotated[i][j] = (char *) calloc(natoms, sizeof(char));
    }
  }
  
  /* rotate the short end of the chain for each residue */
  /* and determine which atoms were rotated, for either phi or psi rotation */

  for (i=0; i<nresidues; i++) {
    
    rotate_natoms[0][i]=0;
    rotate_natoms[1][i]=0;
    
    if (i > nresidues/2.0) { 
      for (j=0; j<natoms; j++)
	if (native[j].res_num > i) {
	  rotate_atom[0][i][rotate_natoms[0][i]++] = j;
	  rotate_atom[1][i][rotate_natoms[1][i]++] = j;
	}	
	else if (native[j].res_num == i) {
	  if (j != native_residue[i].N && j!= native_residue[i].CA)
	    rotate_atom[1][i][rotate_natoms[1][i]++] = j;
	  if (j == native_residue[i].O)
	    rotate_atom[0][i][rotate_natoms[0][i]++] = j;
	}
    }
    else {
      for (j=0; j<natoms; j++)
	if (native[j].res_num < i) {
	  rotate_atom[0][i][rotate_natoms[0][i]++] = j;
	  rotate_atom[1][i][rotate_natoms[1][i]++] = j;
	}
	else if (native[j].res_num == i) {
	  if (j != native_residue[i].O && j != native_residue[i].C && j != native_residue[i].CA)
	    rotate_atom[0][i][rotate_natoms[0][i]++] = j;
	}
    } 
  }

  /* set up the not_rotated array, with 1's at every unrotated atom */
  
  for (i=0; i<nresidues; i++) {
    for (j=0; j<rotate_natoms[0][i]; j++)
      not_rotated[0][i][rotate_atom[0][i][j]]=1;
    for (j=0; j<rotate_natoms[1][i]; j++)
      not_rotated[1][i][rotate_atom[1][i][j]]=1;
    for (j=0; j<natoms; j++) {
      not_rotated[0][i][j]=!not_rotated[0][i][j];
      not_rotated[1][i][j]=!not_rotated[1][i][j];
    }
  }

  /* store atoms rotated by loop moves */

  loop_rotate_natoms = (short **) calloc(total_triplets, sizeof(short *));
  loop_int_rotate_natoms = (short **) calloc(total_triplets, sizeof(short *));
  loop_rotate_atoms = (short ***) calloc(total_triplets, sizeof(short **));
  loop_int_rotate_atoms = (short ***) calloc(total_triplets, sizeof(short **));
  loop_not_rotated = (char ***) calloc(total_triplets, sizeof(char **));
  for (i=0; i<total_triplets; i++){
    /* there are up to 6 bonds to rotate for each loop move */
    loop_rotate_natoms[i] = (short *) calloc(6, sizeof(short));
    loop_int_rotate_natoms[i] = (short *) calloc(6, sizeof(short));
    loop_rotate_atoms[i] = (short **) calloc(6, sizeof(short *));
    loop_int_rotate_atoms[i] = (short **) calloc(6, sizeof(short *));
    loop_not_rotated[i] = (char **) calloc(6, sizeof(char *));

    /* all loop moves have at least 2 bonds */
    loop_rotate_atoms[i][0] = (short *) calloc(natoms, sizeof(short));
    loop_rotate_atoms[i][1] = (short *) calloc(natoms, sizeof(short));
    loop_int_rotate_atoms[i][0] = (short *) calloc(natoms, sizeof(short));  
    loop_int_rotate_atoms[i][1] = (short *) calloc(natoms, sizeof(short));   
    loop_not_rotated[i][0] = (char *) calloc(natoms, sizeof(char)); 
    loop_not_rotated[i][1] = (char *) calloc(natoms, sizeof(char)); 

    if (residue_triplets[i].b >=0) {
      /* double and triple loop moves */
      loop_rotate_atoms[i][2] = (short *) calloc(natoms, sizeof(short));
      loop_rotate_atoms[i][3] = (short *) calloc(natoms, sizeof(short)); 
      loop_int_rotate_atoms[i][2] = (short *) calloc(natoms, sizeof(short)); 
      loop_int_rotate_atoms[i][3] = (short *) calloc(natoms, sizeof(short)); 
      loop_not_rotated[i][2] = (char *) calloc(natoms, sizeof(char)); 
      loop_not_rotated[i][3] = (char *) calloc(natoms, sizeof(char)); 
    }
    else {
      /* single loop moves */
      loop_rotate_atoms[i][2] = (short *) calloc(1, sizeof(short));  /* why are these allocated? */
      loop_rotate_atoms[i][3] = (short *) calloc(1, sizeof(short)); 
      loop_int_rotate_atoms[i][2] = (short *) calloc(1, sizeof(short)); 
      loop_int_rotate_atoms[i][3] = (short *) calloc(1, sizeof(short)); 
      loop_not_rotated[i][2] = (char *) calloc(natoms, sizeof(char)); 
      loop_not_rotated[i][3] = (char *) calloc(natoms, sizeof(char)); 
    }

    if (residue_triplets[i].c >=0) {
      /* triple loop moves */
      loop_rotate_atoms[i][4] = (short *) calloc(natoms, sizeof(short)); 
      loop_rotate_atoms[i][5] = (short *) calloc(natoms, sizeof(short)); 
      loop_int_rotate_atoms[i][4] = (short *) calloc(natoms, sizeof(short)); 
      loop_int_rotate_atoms[i][5] = (short *) calloc(natoms, sizeof(short)); 
      loop_not_rotated[i][4] = (char *) calloc(natoms, sizeof(char)); 
      loop_not_rotated[i][5] = (char *) calloc(natoms, sizeof(char)); 
    }
    else {
      /* single and double loop moves */
      loop_rotate_atoms[i][4] = (short *) calloc(1, sizeof(short)); 
      loop_rotate_atoms[i][5] = (short *) calloc(1, sizeof(short)); 
      loop_int_rotate_atoms[i][4] = (short *) calloc(1, sizeof(short)); 
      loop_int_rotate_atoms[i][5] = (short *) calloc(1, sizeof(short)); 
      loop_not_rotated[i][4] = (char *) calloc(natoms, sizeof(char)); 
      loop_not_rotated[i][5] = (char *) calloc(natoms, sizeof(char)); 
    }
  }
  
  for (i=0; i<total_triplets; i++) {

    if (residue_triplets[i].a > nresidues/2.0) {  /* phi then psi */

      loop_rotate_natoms[i][0] = rotate_natoms[1][residue_triplets[i].a]; /* phi */
      for (j=0; j < loop_rotate_natoms[i][0]; j++) 
	loop_rotate_atoms[i][0][j] = rotate_atom[1][residue_triplets[i].a][j];
      loop_rotate_natoms[i][1] = rotate_natoms[0][residue_triplets[i].a]; /* psi */
      for (j=0; j < loop_rotate_natoms[i][1]; j++) 
	loop_rotate_atoms[i][1][j] = rotate_atom[0][residue_triplets[i].a][j];

      if (residue_triplets[i].b >=0) {
	loop_rotate_natoms[i][2] = rotate_natoms[1][residue_triplets[i].b]; /* phi */
	for (j=0; j < loop_rotate_natoms[i][2]; j++) 
	  loop_rotate_atoms[i][2][j] = rotate_atom[1][residue_triplets[i].b][j];
	loop_rotate_natoms[i][3] = rotate_natoms[0][residue_triplets[i].b]; /* psi */
	for (j=0; j < loop_rotate_natoms[i][3]; j++) 
	  loop_rotate_atoms[i][3][j] = rotate_atom[0][residue_triplets[i].b][j];
      }
      else {
	loop_rotate_natoms[i][2] = 0;
	loop_rotate_natoms[i][3] = 0;
      }
      if (residue_triplets[i].c >=0) {
	loop_rotate_natoms[i][4] = rotate_natoms[1][residue_triplets[i].c]; /* phi */
	for (j=0; j < loop_rotate_natoms[i][4]; j++) 
	  loop_rotate_atoms[i][4][j] = rotate_atom[1][residue_triplets[i].c][j];
	loop_rotate_natoms[i][5] = rotate_natoms[0][residue_triplets[i].c]; /* psi */
	for (j=0; j < loop_rotate_natoms[i][5]; j++) 
	  loop_rotate_atoms[i][5][j] = rotate_atom[0][residue_triplets[i].c][j];
      }
      else {
	loop_rotate_natoms[i][4] = 0;
	loop_rotate_natoms[i][5] = 0;
      }
    }
    else { /* psi then phi */
      j=0;
      if (residue_triplets[i].c >=0) {
	for (k=0; k<natoms; k++)
	  if (native[k].res_num < residue_triplets[i].c) {
	    loop_rotate_atoms[i][j][loop_rotate_natoms[i][j]++] = k;   /* psi */
	    loop_rotate_atoms[i][j+1][loop_rotate_natoms[i][j+1]++] = k;  /* phi */
	  }
	  else if (native[k].res_num == residue_triplets[i].c) {
	    if (k != native_residue[residue_triplets[i].c].O && k != native_residue[residue_triplets[i].c].C && k != native_residue[residue_triplets[i].c].CA)
	      loop_rotate_atoms[i][j][loop_rotate_natoms[i][j]++] = k; /* psi */
	  }
	j+=2;
      }
      else {
	loop_rotate_natoms[i][4] = 0;
	loop_rotate_natoms[i][5] = 0;
      }
      if (residue_triplets[i].b >=0) {
 	for (k=0; k<natoms; k++)
	  if (native[k].res_num < residue_triplets[i].b) {
	    loop_rotate_atoms[i][j][loop_rotate_natoms[i][j]++] = k;  /* psi */
	    loop_rotate_atoms[i][j+1][loop_rotate_natoms[i][j+1]++] = k; /* phi */
	  }
	  else if (native[k].res_num == residue_triplets[i].b) {
	    if (k != native_residue[residue_triplets[i].b].O && k != native_residue[residue_triplets[i].b].C && k != native_residue[residue_triplets[i].b].CA)
	      loop_rotate_atoms[i][j][loop_rotate_natoms[i][j]++] = k; /* psi */
	  }
	j+=2;
      }
      else {
	loop_rotate_natoms[i][2] = 0;
	loop_rotate_natoms[i][3] = 0;
      }
      loop_rotate_natoms[i][j] = rotate_natoms[0][residue_triplets[i].a]; /* psi */
      for (k=0; k < loop_rotate_natoms[i][j]; k++) 
	loop_rotate_atoms[i][j][k] = rotate_atom[0][residue_triplets[i].a][k];
      loop_rotate_natoms[i][j+1] = rotate_natoms[1][residue_triplets[i].a]; /* phi */
      for (k=0; k < loop_rotate_natoms[i][j+1]; k++) 
	loop_rotate_atoms[i][j+1][k] = rotate_atom[1][residue_triplets[i].a][k];
    }
  }

  for (i=0; i<total_triplets; i++) 
    for (j=0; j<6; j++) {
      for (k=0; k<loop_rotate_natoms[i][j]; k++)
	loop_not_rotated[i][j][loop_rotate_atoms[i][j][k]]=1;
      for (k=0; k<natoms; k++)
	loop_not_rotated[i][j][k]=!loop_not_rotated[i][j][k];
    }
  
  for (i=0; i<total_triplets; i++)
    for (j=5; j>0; j--) {
      if (loop_rotate_natoms[i][j] && loop_rotate_natoms[i][j-1]) {
	for (k=0; k<loop_rotate_natoms[i][j-1]; k++)
	  if (loop_not_rotated[i][j][loop_rotate_atoms[i][j-1][k]])
	    loop_int_rotate_atoms[i][j-1][loop_int_rotate_natoms[i][j-1]++]=loop_rotate_atoms[i][j-1][k];
      }
    }
      
  return;

}

void InitializeSidechainRotationData() {
  int i, j, k, l, m, n;

  sidechain_torsion = (short ***) calloc(nresidues, sizeof(short **));
  for (i=0; i<nresidues; i++) {
    sidechain_torsion[i] = (short **) calloc(amino_acids[native_residue[i].amino_num].ntorsions, sizeof(short *));
    for (j=0; j<amino_acids[native_residue[i].amino_num].ntorsions; j++)
      sidechain_torsion[i][j] = (short *) calloc(4, sizeof(short));
  }

  sct_E = (short *****) calloc(nresidues, sizeof(short ****));
  for (i=0; i<nresidues; i++)
   {
    sct_E[i] = (short ****) calloc(12, sizeof(short ***));
    for (j=0; j<12; j++)
     {
      sct_E[i][j] = (short ***) calloc(12, sizeof(short **));
      for (k=0; k<12; k++)
       {
        sct_E[i][j][k] = (short **) calloc(12, sizeof(short *));
        for (l=0; l<12; l++)
         {
          sct_E[i][j][k][l] = (short *) calloc(12, sizeof(short));
         }
       }
     }
   }

  for (i=0; i<nresidues; i++) {
    native_residue[i].ntorsions = amino_acids[native_residue[i].amino_num].ntorsions;
    native_residue[i].nrotamers = amino_acids[native_residue[i].amino_num].nrotamers;

    /* copy rotamer angles from amino_acid structure into native_residue structure */
    for (j=0; j<4; j++)
      for (k=0; k<4; k++)
	for (l=0; l<4; l++)
	  for (m=0; m<4; m++)
	    for (n=0; n<4; n++)
	      native_residue[i].avg_angle[j][k][l][m][n] = amino_acids[native_residue[i].amino_num].avg_angle[j][k][l][m][n];

    /* rot_position gives the 4-digit base 3 representation of the rotamer */

    native_residue[i].rot_position = (short **) calloc(native_residue[i].nrotamers, sizeof(short *));
    for (j=0; j<native_residue[i].nrotamers; j++)
      native_residue[i].rot_position[j] = (short *) calloc(4, sizeof(short));
        
    if (native_residue[i].nrotamers>1)
      for (k=0; k<native_residue[i].nrotamers; k++) {
	j=native_residue[i].ntorsions;
	l = k;
	do {
	  native_residue[i].rot_position[k][j-1] = l/three[j-1];
	  l -= three[j-1]*native_residue[i].rot_position[k][j-1];
	  native_residue[i].rot_position[k][j-1] +=1;
	  j--;
	} while (j>0); 
      }
  }

  /* these structures keep track of which sidechain atoms rotate at each torsion */

  rotate_sidechain_atom = (short ***) calloc(nresidues, sizeof(short **));
  sidechain_not_rotated = (char ***) calloc(nresidues, sizeof(char **));
  for (i=0; i<nresidues; i++) {
    rotate_sidechain_atom[i] = (short **) calloc(native_residue[i].ntorsions, sizeof(short *));
    sidechain_not_rotated[i] = (char **) calloc(native_residue[i].ntorsions, sizeof(char *));
    for (j=0; j<native_residue[i].ntorsions; j++) {
      rotate_sidechain_atom[i][j] = (short *) calloc(amino_acids[native_residue[i].amino_num].rotate_natoms[j], sizeof(short *));
      sidechain_not_rotated[i][j] = (char *) calloc(natoms, sizeof(char));
    }
  }
  
  rotate_sidechain_natoms = (short **) calloc(nresidues, sizeof(short *));
  for (i=0; i<nresidues; i++) 
    rotate_sidechain_natoms[i] = (short *) calloc(native_residue[i].ntorsions, sizeof(short));
  for (i=0; i<nresidues; i++) {
    for (j=0; j<native_residue[i].ntorsions; j++)
      rotate_sidechain_natoms[i][j] = amino_acids[native_residue[i].amino_num].rotate_natoms[j];
    for (j=0; j<native_residue[i].ntorsions; j++)
      for (k=0; k<4; k++) {
	for (l=0; l<natoms; l++) 
	  if (native[l].res_num == i && !strcmp(native[l].atomname, amino_acids[native_residue[i].amino_num].torsion[j][k])) {
	    /* sidechain_torsion records the 4 atoms that define each torsion */
	    sidechain_torsion[i][j][k] = l; 
	    break;
	  }
	if (l==natoms)
         {
	  fprintf(STATUS, "WARNING -- atom %s, residue %d  %s not found!\n", amino_acids[native_residue[i].amino_num].torsion[j][k], i, native_residue[i].res);
          exit(1);
         }
      }
    for (j=0; j<native_residue[i].ntorsions; j++)
      for (k=0; k<rotate_sidechain_natoms[i][j]; k++) {
	for (l=0; l<natoms; l++) 
	  if (native[l].res_num == i && !strcmp(native[l].atomname, amino_acids[native_residue[i].amino_num].rotate_atom[j][k])) {
	    rotate_sidechain_atom[i][j][k] = l;
	    break;
	  } 
	if (l==natoms)
         {
	  fprintf(STATUS, "WARNING -- atom %s, residue %d  %s not found!\n", amino_acids[native_residue[i].amino_num].rotate_atom[j][k], i, native_residue[i].res);
          exit(1);
         }
      }
  }
  
  for (i=0; i<nresidues; i++)
    for (j=0; j< native_residue[i].ntorsions; j++) {
      for (k=0; k<natoms; k++)
	sidechain_not_rotated[i][j][k]=1;
      for (k=0; k<rotate_sidechain_natoms[i][j]; k++)
	sidechain_not_rotated[i][j][rotate_sidechain_atom[i][j][k]]=0;
    }
  
  return;
}

/*=============================================================*/
/*               set up correlation data                       */
/*=============================================================*/

int SkipSelf(int s, int b, struct atom *Protein, struct residue *Residue) {
  /* returns 1 if sidechain-backbone atom pair is separated by less than 3 bonds, or 0 if not */

  if ((b == Residue[Protein[b].res_num].C || b == Residue[Protein[b].res_num].N || b == Residue[Protein[b].res_num].CA)  && s == Residue[Protein[s].res_num].CB)
    return 1;
  else if (Residue[Protein[s].res_num].amino_num==14)
    return 1;
  else if (b == Residue[Protein[b].res_num].CA && Protein[s].atomname[1] == 'G') 
    return 1;
  else
    return 0;

}
  
int SkipNeighbors(int i, int j, struct atom *Protein, struct residue *Residue) {

  int first, second;
  
  if (Protein[i].res_num < Protein[j].res_num) {
    first = i;
    second = j;
  }
  else {
    second = i;
    first = j;
  }
  if (!(first == Residue[Protein[first].res_num].C && !strcmp(Protein[second].atomname, "CD") && Residue[Protein[second].res_num].amino_num == 14) && !(first == Residue[Protein[first].res_num].CA && !strcmp(Protein[second].atomname, "CD") && Residue[Protein[second].res_num].amino_num == 14)) 
    if ((first == Residue[Protein[first].res_num].N) || (second != Residue[Protein[second].res_num].CA && second != Residue[Protein[second].res_num].N) || Protein[first].is_sidechain || Protein[second].is_sidechain)
      return 0;
  
  return 1;

}

int Disulfide(int a, int b, struct atom *Protein, struct residue *Residue) {

  if ((Residue[Protein[a].res_num].amino_num == 4) && (Residue[Protein[b].res_num].amino_num == 4) && !strcmp(Protein[a].atomname, "SG") && !strcmp(Protein[b].atomname, "SG"))
    return 1;
  else
    return 0;

}

void CheckCorrelation(struct contact_data **Data, struct atom *Protein, struct residue *Residue, int Natoms) {
  int i, j;
  /*====================================================================================*/
  /*  non-local residues     check clashes and contacts of all non-disulfide pairs      */
  /*  non-local disulfides   check only contacts of S-S atom pairs                      */
  /*  self                   check clashes of only non-bonded sidechain-backbone pairs  */
  /*  i-i+1                  check clashes of non-bonded pairs                          */
  /*  i-i+2                  check clashes of all pairs                                 */
  /*  relevant bb atoms      check clashes of all pairs                                 */
  /*====================================================================================*/

  for (i=0; i<Natoms; i++)
    for (j=i+1; j<Natoms; j++) {
      
      /* non-local residues */
      
      if (fabs(Protein[i].res_num-Protein[j].res_num)>SKIP_LOCAL_CONTACT_RANGE) {
	if (!Disulfide(i, j, Protein, Residue)) { 
	  Data[i][j].check_clashes=1;
	  Data[j][i].check_clashes=1;
	  Data[i][j].check_contacts=1;
	  Data[j][i].check_contacts=1;
	}
     	else {
	  /* disulfide */
	  Data[i][j].disulfide=1;
	  Data[j][i].disulfide=1;
	  Data[i][j].check_clashes=0;
	  Data[j][i].check_clashes=0;
	  Data[i][j].check_contacts=1;
	  Data[j][i].check_contacts=1;
	}
      }

      /* self */
      
      else if (Protein[i].res_num == Protein[j].res_num) {
	if (Protein[i].is_sidechain && !Protein[j].is_sidechain) {  /* sidechain - backbone */
 	  if (!SkipSelf(i, j, Protein, Residue)) {
	    Data[i][j].check_clashes=1;
	    Data[j][i].check_clashes=1;
	    Data[i][j].check_contacts=0;
	    Data[j][i].check_contacts=0;
	  }
	}
	else if (!Protein[i].is_sidechain && Protein[j].is_sidechain) {  /* backbone - sidechain */
	  if (!SkipSelf(j, i, Protein, Residue)) {
	    Data[i][j].check_clashes=1;
	    Data[j][i].check_clashes=1;
	    Data[i][j].check_contacts=0;
	    Data[j][i].check_contacts=0;
	  }
	}
      }
     
      /* i-i+1 */

      else if (fabs(Protein[i].res_num-Protein[j].res_num)==1) { 
	if (!SkipNeighbors(i, j, Protein, Residue))  {
	  Data[i][j].check_clashes=1;
	  Data[j][i].check_clashes=1;
	  Data[i][j].check_contacts=0;
	  Data[j][i].check_contacts=0;
	}
      }
      
      /* all other local residues */

      else if (fabs(Protein[i].res_num-Protein[j].res_num)<=SKIP_LOCAL_CONTACT_RANGE && fabs(Protein[i].res_num-Protein[j].res_num)>=2) {  
	Data[i][j].check_clashes=1;
	Data[j][i].check_clashes=1;
	Data[i][j].check_contacts=0;
	Data[j][i].check_contacts=0;
      }

      else {
	fprintf(STATUS, "Uncategorizable pair\n");
	exit(0);
      }

      /* Skip relevant backbone contacts */

      if (!IsSidechainAtom(Protein[i].atomname) && !IsSidechainAtom(Protein[j].atomname) && fabs(Protein[i].res_num - Protein[j].res_num) <= SKIP_BB_CONTACT_RANGE) {
	Data[i][j].check_contacts=0;
	Data[j][i].check_contacts=0;
      }

      /* up to here, backbone contacts below the SKIP_BB_CONTACT_RANGE are turned off */
      
    }
  
  return;
}


void SetProgramOptions(int argc, char *argv[]) {
  
  char line[500];		/* increased to 500 from 150 */
  char token[50];
  char name[500];		/* increased to 500 from 50 */
  float value;
  int find_yang_move = 0;
  int find_yang_scale = 0;

  char cfg_file[200];
  int l=0, ls, MPI_STOP=0;

  /* Acquire cfg_file name */
  memset(cfg_file, '\0', 200);
  if (myrank == 0) {
    strcpy(cfg_file, argv[1]);
    if (argc != 2) {
      printf("ERROR!!! Usage is like this: ./fold_potential config_file, argc : %d\n", argc);
      for (l=0; l<argc; l++){
        printf("argc : %3d, argv : %s\n", l, argv[l]);
      }
      MPI_STOP = 1;
    }
  }
	
  ierr = MPI_Bcast(&MPI_STOP, 1, MPI_INT, 0, mpi_world_comm);

  if (MPI_STOP == 1) {
    MPI_Finalize();
    exit(1);
  } else {
    ierr = MPI_Bcast(cfg_file, 200, MPI_CHAR, 0, mpi_world_comm);
  }

  /* open cfg file */
  
  if ((DATA = fopen(cfg_file, "r"))==NULL) {
    printf("ERROR: Can't open the file: %s!\n", cfg_file);
    exit(1);
  }

  /* acquire params from cfg file */
  /* increased from 150 to 500 */
  while (fgets(line, 500, DATA) != NULL) {
    /* printf("%s", line);*/
    if (strncmp(line, "!", 1)) {
      sscanf(line, "%s %f", token, &value);
      if (!strcmp(token, "NATIVE_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(native_file, name);
      }
      else if (!strcmp(token, "STRUCTURE_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(structure_file, name);
      }
      else if (!strcmp(token, "NATIVE_DIRECTORY")){
	sscanf(line, "%*s %s", name);
	strcpy(native_directory, name);     
	  }
      else if (!strcmp(token, "TEMPLATE_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(template_file, name);
      }
      else if (!strcmp(token, "SUBSTRUCTURE_FILE")){
      sscanf(line, "%*s %s", name);
      strcpy(substructure_path, name);
      }
      else if (!strcmp(token, "ALIGNMENT_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(alignment_file, name);
      }
      else if (!strcmp(token, "TRIPLET_ENERGY_FILE")) {
        sscanf(line, "%*s %s", name);
        strcpy(triplet_file, name);
      }
      else if (!strcmp(token, "SIDECHAIN_TORSION_FILE")) {
        sscanf(line, "%*s %s", name);
        strcpy(sctorsion_file, name);
      }
      else if (!strcmp(token, "SECONDARY_STRUCTURE_FILE")) {
        sscanf(line, "%*s %s", name);
        strcpy(sec_str_file, name);
      }
      else if (!strcmp(token, "AMINO_DATA_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(amino_data_file, name);
      }
      else if (!strcmp(token, "ATOM_TYPE_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(atom_type_file, name);
      }
      else if (!strcmp(token, "ROTAMER_DATA_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(rotamer_data_file, name);
      }
      else if (!strcmp(token, "PDB_OUT_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(pdb_out_file, name);
	strcpy(std_prefix, name);
 	ls = strlen(name);
        /*sprintf(pdb_out_file+ls, "_%5.3f", MC_TEMP);*/
      }
      else if (!strcmp(token, "POTENTIAL_DATA")) {
	sscanf(line, "%*s %s", name);
	strcpy(potential_file, name);
      } 
      else if (!strcmp(token, "HELICITY_DATA")) {
	sscanf(line, "%*s %s", name);
	strcpy(helicity_data, name);
      } 
      else if (!strcmp(token, "HYDROGEN_BONDING_DATA")) {
	sscanf(line, "%*s %s", name);
	strcpy(hydrogen_bonding_data, name);
      } 
      else if (!strcmp(token, "AROMATIC_FILE")) {
	sscanf(line, "%*s %s", name);
	strcpy(aromatic_file, name);
      } 
      else if (!strcmp(token, "PROTEIN_NAME")) {
	sscanf(line, "%*s %s", name);
	strcpy(PROTEIN_NAME, name);
      }
      else if (!strcmp(token, "PRINT_PDB"))
	PRINT_PDB = (int) value;
      else if (!strcmp(token, "MC_STEPS"))
	MC_STEPS = (long int) value;
      else if (!strcmp(token, "MC_ANNEAL_STEPS"))
	MC_ANNEAL_STEPS = (long int) value;
      else if (!strcmp(token, "MC_STEP_SIZE"))
	STEP_SIZE =  value*PI/180.0;
      else if (!strcmp(token, "SIDECHAIN_NOISE"))
	SIDECHAIN_NOISE = value*PI/180.0;
      else if (!strcmp(token, "SIDECHAIN_MOVES"))
        SIDECHAIN_MOVES = (int) value;
      else if (!strcmp(token, "MC_PRINT_STEPS"))
	MC_PRINT_STEPS = (int) value;
      else if (!strcmp(token, "MC_PDB_PRINT_STEPS"))
	MC_PDB_PRINT_STEPS = (int) value;
      else if (!strcmp(token, "MC_TEMP_MIN"))
	MC_TEMP_MIN = value;
      else if (!strcmp(token, "TEMP_STEP"))
	TEMP_STEP = value;
      else if (!strcmp(token, "NODES_PER_TEMP"))
	NODES_PER_TEMP = (int) value;
      else if (!strcmp(token, "ALPHA"))
	ALPHA = value;      
      else if (!strcmp(token, "LAMBDA"))
	LAMBDA = value;      
      else if (!strcmp(token, "NATIVE_ATTRACTION"))
	NATIVE_ATTRACTION = value;      
      else if (!strcmp(token, "NON_NATIVE_REPULSION"))
	NON_NATIVE_REPULSION = value;
      else if (!strcmp(token, "USE_ROT_PROB"))
        USE_ROT_PROB = value;       
      else if (!strcmp(token, "SEQ_DEP_HB"))
        SEQ_DEP_HB = value;
      else if (!strcmp(token, "CLASH_WEIGHT"))
	weight_clash = value;       
      else if (!strcmp(token, "HYDROGEN_BOND"))
	hydrogen_bond = value;             
      else if (!strcmp(token, "RMS_WEIGHT"))
	weight_rms = value;              
      else if (!strcmp(token, "CONSTRAINT_RMSD"))
	rmsd_constraint = value;              
      else if (!strcmp(token, "NON_SPECIFIC_ENERGY"))
	NON_SPECIFIC_ENERGY = value;         
      else if (!strcmp(token, "USE_GLOBAL_BB_MOVES"))
	USE_GLOBAL_BB_MOVES = value;
	  else if (!strcmp(token, "CONSTRAINT_FILE")){ //AB
	sscanf(line, "%*s %s", name);
	strcpy(constraint_file, name);
	}
	  else if (!strcmp(token, "K_CONSTRAINT")) //AB
	k_constraint = value;
      else if (!strcmp(token, "YANG_MOVE"))
       {
	YANG_MOVE = value;
	find_yang_move = 1;
       }
      else if (!strcmp(token, "YANG_SCALE"))
       {
	YANG_SCALE = value;
	find_yang_scale = 1;
       }
      else if (!strcmp(token, "SKIP_LOCAL_CONTACT_RANGE"))
	SKIP_LOCAL_CONTACT_RANGE = (int) value;      
      else if (!strcmp(token, "SKIP_BB_CONTACT_RANGE"))
        SKIP_BB_CONTACT_RANGE = (int) value;  
      else if  (!strcmp(token, "MY_RANK_OFFSET"))
      	my_rank_offset = (int) value;
      else if (!strcmp(token, "USE_SIDECHAINS"))
        USE_SIDECHAINS = (int) value;      
      else if (!strcmp(token, "NO_NEW_CLASHES"))
        NO_NEW_CLASHES = (int) value;      
      else if (!strcmp(token, "USE_ROTAMERS"))
        USE_ROTAMERS = (int) value;      
      else if (!strcmp(token, "READ_POTENTIAL"))
	READ_POTENTIAL = (int) value;      
      else if (!strcmp(token, "USE_GO_POTENTIAL"))
	USE_GO_POTENTIAL = (int) value;      
      else if (!strcmp(token, "MC_REPLICA_STEPS"))
        MC_REPLICA_STEPS = (int) value;
      else if (!strcmp(token, "MAX_EXCHANGE"))
        MAX_EXCHANGE = (int) value;
      else if (!strcmp(token, "UMBRELLA"))
        umbrella = (int) value;
      else if (!strcmp(token, "K_BIAS"))
        k_bias =  value;   
      else if (!strcmp(token, "NUMBER_OF_CONTACTS_MAX"))
        number_of_contacts_max =  (int) value;   
      else if (!strcmp(token, "CONTACTS_STEP"))
        contacts_step = (int) value;  
      else if (!strcmp(token, "MIN_SEQ_SEP"))
      	min_seq_sep = (int) value; 
      else if (!strcmp(token, "CONTACT_CALPHA_CUTOFF"))
      	contact_calpha_cutoff = value; //note that by default, this value is set to 7 in backbone.h, unless you enter something in cfg file
      else if (!strcmp(token, "USE_CLUSTER")) //With what probability should code use a knowledge based move?
      	USE_CLUSTER = value; //note that by default, this value is set to 0 in backbone.h, unless you enter something in cfg file
      else if (!strcmp(token, "MAX_CLUSTERSTEP")) //As of what MC should knowledge based moves be stopped?
      	MAX_CLUSTERSTEP = value; //note that by default, this value is set to 0 in backbone.h, unless you enter something in cfg file
      else if (!strcmp(token, "CLUSTER_MOVE")) //How often do you call the loopmove function? This does not mean you do a knowledge move
      	CLUSTER_MOVE = value; //note that by default, this value is set to 0.33 in backbone.h, unless you enter something in cfg file
	  else {
		printf("config file option not found: %s\n", token);
	//exit(0);
      }
    }
  }
  
  fclose(DATA);			/* close cfg file */

	/*Here's what we gonna do
	If we have entered in a native directory, then the files in this directory will override the native file that is entered here
	The files in the native directory need to be in the form my_rank.pdb/
	*/
  
  //int mod_rank=(myrank+my_rank_offset) % NODES_PER_TEMP;
  //char temp1[500];
  
  if (strcmp(native_directory, "None")!=0) {
  	sprintf(native_file, "%s/%d.pdb", native_directory, myrank+my_rank_offset);
  	//printf("%s", native_directory);
  	
  	/*
  	Now, suppose we have not specified one starting file for each core
  	Then we repeat starting files periodically every nodes_per_temp
  	
  	So if nodes_per_temp=10, then cores 0, 10, 20, etc will all use 0.pdb
  	*/

  	FILE *test_file;
  	if ((test_file = fopen(native_file, "r"))==NULL) {  
  		int mod_rank = (myrank+my_rank_offset) % NODES_PER_TEMP;
  		//printf("The current value of mod_rank is %i", mod_rank);
  		sprintf(native_file, "%s/%d.pdb", native_directory, mod_rank);
  		//printf("The current value of native_file is %s", native_file);
  		
  	}
  }

  if (nprocs % NODES_PER_TEMP!=0) {
  	printf("ERROR! Number of cores must be divisible by NODES_PER_TEMP");
  	exit(0);
  }
  
  /*
  Alright bro; here's the deal
  The temperatures and setpoints will be assigned as follows:
  
  Suppose nodes_per_temp=5, contacts_max=300, contacts_step=10, mc_temp_min=0.1, temp_step=0.1, and n_procs=20
  Then the temperautre assignments will be as follows
  
  	Cnode[0]=300, Tnode[0]=0.100000                                                                 
	Cnode[1]=290, Tnode[1]=0.100000                                                                 
	Cnode[2]=280, Tnode[2]=0.100000                                                                 
	Cnode[3]=270, Tnode[3]=0.100000                                                                 
	Cnode[4]=260, Tnode[4]=0.100000                                                                 
	Cnode[5]=300, Tnode[5]=0.200000  
	Cnode[6]=290, Tnode[6]=0.20000
	....
	...
	...
	Cnode[15]=300, Tnode[15]=0.400000                                                               
	Cnode[16]=290, Tnode[16]=0.400000                                                               
	Cnode[17]=280, Tnode[17]=0.400000                                                               
	Cnode[18]=270, Tnode[18]=0.400000                                                               
	Cnode[19]=260, Tnode[19]=0.400000   
  
  That is, at every node, we decrease the contacts setpoint by contacts_step
  But every nodes_per_temp nodes, we increment the temperature by temp_step, and reset the contacts_setpoint to contacts_max
  */
  
  
  /* establish some ladder/replica exchange items */
  /*Note for C newbies such as AB: The call
  
  (float *) calloc(nprocs, sizeof(float));
  
  Means that we are allocating an array of floats with nprocs elements
  The function calloc allocates the necessary amount of memory (which is given by
  nprocs times the memory per array element) and sets those elements to zero
  The (float*) in front of it indicates that we are asking for a pointer to an
  array of floats
  */
  Tnode = (float *) calloc(nprocs, sizeof(float));
  Enode = (float *) calloc(nprocs, sizeof(float));
  Cnode = (int *) calloc(nprocs, sizeof(int)); //number of contacts setpoints
  Nnode = (int *) calloc(nprocs, sizeof(int));
  
  float current_temp=MC_TEMP_MIN; 
  int current_setpoint;
  for (l=0; l<nprocs; l++) {
  	current_setpoint=number_of_contacts_max-l%NODES_PER_TEMP*contacts_step;
  	Cnode[l]=current_setpoint;
    Tnode[l]=current_temp;
    //current_setpoint=current_setpoint-contacts_step;
    if ((l+1) % NODES_PER_TEMP==0){
    	current_temp=current_temp+TEMP_STEP;
    }
  }
  

  MC_TEMP=Tnode[myrank];   //temperature must always be externally modified if we want ot change the my_rank_offset
  number_of_contacts_setpoint=Cnode[myrank]; ///same with setpoint


  /* SET name of PDB and log file  */
  
  if (umbrella==1){
  	sprintf(std_file, "%s_T_%5.3f_%d.log", std_prefix, MC_TEMP, number_of_contacts_setpoint);
  	sprintf(pdb_out_file+ls, "_%5.3f_%d", MC_TEMP, number_of_contacts_setpoint);
  	}
  else {
  	k_bias=0;  //regardless of what we had inputted for k_bias, if umbrella is off, we don't want any bias!!
  	if (NODES_PER_TEMP==1 && my_rank_offset==0){
  		sprintf(std_file, "%s_T_%5.3f.log", std_prefix, MC_TEMP);
  		sprintf(pdb_out_file+ls, "_%5.3f", MC_TEMP);
  	}
  	else {
  		sprintf(std_file,"%s_T_%5.3f_%d.log",std_prefix, MC_TEMP, myrank+my_rank_offset);
  		sprintf(pdb_out_file+ls,"_%5.3f_%d", MC_TEMP, myrank+my_rank_offset);
  	}
  }

  /* OPEN log file  for business*/
  /*sprintf(std_file, "%s_%5.3f.log", std_prefix, MC_TEMP);*/
  STATUS = fopen(std_file, "w");	/* for some reason, code gives file handle the name STATUS */
  //fprintf(STATUS,"The value of temp1 is %s and the value of mod_rank is %i", temp1, mod_rank);
  //fflush(STATUS);
  
  



  replica_index = (int *) calloc(nprocs, sizeof(int));
  
  
  //AB attempted the following
  //accepted_replica = (int **) calloc(nprocs, sizeof(int*)); //AB made this 2D array 
  //rejected_replica = (int **) calloc(nprocs, sizeof(int*));
  //int i;
  //for (i=0; i<nprocs; i++) {
        //accepted_replica[i]=(int*) calloc(nprocs, sizeof(int));
        //rejected_replica[i]=(int*) calloc(nprocs, sizeof(int));
        //for (l=0; l<nprocs; l++) {
        //	accepted_replica[i][l] = rejected_replica[i][l] = 0;
 		//}
    //}
  
  //Previously:
  accepted_replica = (int *) calloc(nprocs, sizeof(int));
  rejected_replica = (int *) calloc(nprocs, sizeof(int));
  
  for (l=0; l<nprocs; l++) {
  	accepted_replica[l] = rejected_replica[l] = 0;
   }


  fprintf(STATUS, "The native directory is %s \n", native_directory);
  fprintf(STATUS, "Temperature/setpoint range for replica exchange!\n");
  
  
  for (l=0; l<nprocs; l++) {
    fprintf(STATUS, "%4d : %5.3f %d \n", l, Tnode[l], Cnode[l]);
  }
  
  
  
  fprintf(STATUS, "myrank : %4d, cfg_file : %s\n k_bias : %5.3f\n setpoint : %d\n ", myrank, cfg_file, k_bias, number_of_contacts_setpoint);
  fflush(STATUS);

  if (find_yang_move == 0) {
    fprintf(STATUS, "There is nothing on YANG_MOVE!\n");
    exit(1);
  }
  if (find_yang_scale == 0) {
    fprintf(STATUS, "There is nothing on YANG_SCALE!\n");
    exit(1);
  }

  /* lattice parameters */
  if (DISTANCE_DEPENDENCE)
    LATTICE_SIZE = 1.0/5.25;
  else
    LATTICE_SIZE = 1/(LAMBDA*ALPHA*(1.88+1.88));
  MATRIX_SIZE = 20;
  HALF_MATRIX_SIZE = MATRIX_SIZE/2;


	/*Read the constraints--AB*/
   if (strcmp(constraint_file, "None")!=0 ){
   fprintf(STATUS, "There are constraints to read!\n");
   fflush(STATUS);
   Read_constraints();
   fprintf(STATUS, "Constraints read\n");
   fflush(STATUS);
   }

  return;

}



