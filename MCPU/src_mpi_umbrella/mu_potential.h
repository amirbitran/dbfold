void OutputPotential(void);

void OutputNab (char *file_name) {
  int i, j, k, l;
  int **nat_cons, **nnat_cons;
  FILE *Nab_file;

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

  if((Nab_file = fopen(file_name, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", file_name);
    exit(1);
   }
  for (i = 0; i < MAX_TYPES; ++i)
    for (j = i; j < MAX_TYPES; ++j) {
      fprintf(Nab_file, "%d %d %d %d\n", i, j, nat_cons[i][j], nnat_cons[i][j]); 
    }
  fclose(Nab_file);

  for (i = 0; i < MAX_TYPES; ++i) {
    free(nat_cons[i]);
    free(nnat_cons[i]);
  }
  free(nat_cons);
  free(nnat_cons);
  
}
void OutputPotential (void) {
  FILE *temp_file;
  int i, j;
  char potential_name[200], Nab_file_name[200];

  sprintf(Nab_file_name, "%s_Nab", PROTEIN_NAME);
  OutputNab(Nab_file_name);

  sprintf(potential_name, "%s_potential_%.3f", PROTEIN_NAME, mu);
  if((temp_file = fopen(potential_name, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", potential_name);
    exit(1);
   }
  for (i = 0; i < MAX_TYPES; ++i)
    for (j = i; j < MAX_TYPES; ++j) {
      fprintf(temp_file, "%d %d %.3f\n", i, j, potential[i][j]); 
    }
  fclose(temp_file);
}
