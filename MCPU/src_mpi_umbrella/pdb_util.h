void ParsePDBLine(char *, struct atom *, int *);
int TypeAtom(char *, char *);
int IsSidechainAtom(char *);
int IsDesignedResidue(int);
int IsCoreResidue(int);
int GetAminoNumber(char *);
void PrintPDB(char *);
void PrintPDB_Emin(char *);
void PrintPDB_RMSDmin(char *);
int SMoGType(char *, char *);
int GetSMoGType(char *);
void PrintReplica(char *);
int GetReplica(char *);

void PrintPDB(char *filename) {
  int i;

  if((DATA = fopen(filename, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", filename);
    exit(1);
   }

  fprintf(DATA, "REMARK E %8.3f E_pot %8.3f E_hb %8.3f E_tor %8.3f E_sct %8.3f E_aro %8.3f\n", E, E_pot, E_hbond, E_tor, E_sct, E_aro);

  for (i=0; i<natoms; i++)
    if (strcmp(native[i].atomname, "CB") || strcmp(native[i].res, "GLY"))
      fprintf(DATA, "ATOM%7d  %-3s %3s  %4d%12.3f%8.3f%8.3f  1.00  0.00\n", i+1, native[i].atomname, native[i].res, native[i].res_num+1, native[i].xyz.x, native[i].xyz.y, native[i].xyz.z);
  fprintf(DATA, "TER\n");
  fprintf(DATA, "END\n");

  fclose(DATA);
  
  return;
}
void PrintReplica(char *filename) {
  int i;

  if((DATA = fopen(filename, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", filename);
    exit(1);
   }

  for (i=0; i<natoms; i++)
    if (strcmp(native[i].atomname, "CB") || strcmp(native[i].res, "GLY"))
      fprintf(DATA, "%7d  %f %f %f\n", i, native[i].xyz.x, native[i].xyz.y, native[i].xyz.z);
  fprintf(DATA, "TER\n");

  fclose(DATA);
  
  return;
}
int GetReplica(char *filename) {
  int i;
  char end_line[100];

  if((DATA = fopen(filename, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", filename);
    exit(1);
   }

  for (i=0; i<natoms; i++)
      fscanf(DATA, "%*d  %f%f%f", &native[i].xyz.x, &native[i].xyz.y, &native[i].xyz.z);
  fscanf(DATA, "%s", end_line);
  if(strcmp(end_line, "TER")==0)
   {
    fclose(DATA);
    return 0;
   }
  else
   {
    fclose(DATA);
    return 1;
   }
  
  return 0;
}

void PrintPDB_Emin(char *filename) {
  int i;

  if((DATA = fopen(filename, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", filename);
    exit(1);
   }

  for (i=0; i<natoms; i++)
    if (strcmp(native[i].atomname, "CB") || strcmp(native[i].res, "GLY"))
      fprintf(DATA, "ATOM%7d  %-3s %3s  %4d%12.3f%8.3f%8.3f  1.00  0.00\n", i+1, native[i].atomname, native[i].res, native[i].res_num+1, native_Emin[i].xyz.x, native_Emin[i].xyz.y, native_Emin[i].xyz.z);
  fprintf(DATA, "TER\n");
  fprintf(DATA, "END\n");

  fclose(DATA);

  return;
}

void PrintPDB_RMSDmin(char *filename) {
  int i;

  if((DATA = fopen(filename, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", filename);
    exit(1);
   }

  for (i=0; i<natoms; i++)
    if (strcmp(native[i].atomname, "CB") || strcmp(native[i].res, "GLY"))
      fprintf(DATA, "ATOM%7d  %-3s %3s  %4d%12.3f%8.3f%8.3f  1.00  0.00\n", i+1, native[i].atomname, native[i].res, native[i].res_num+1, native_RMSDmin[i].xyz.x, native_RMSDmin[i].xyz.y, native_RMSDmin[i].xyz.z);
  fprintf(DATA, "TER\n");
  fprintf(DATA, "END\n");

  fclose(DATA);

  return;
}

int SMoGType(char *atomname, char *residue) {
  
  if (!strncmp(atomname, "H", 1))
    return(GetSMoGType("H"));
  else if (!(strcmp(atomname, "C")))
    return(GetSMoGType("CC"));
  else if (!(strcmp(atomname, "CA")))
    return(GetSMoGType("CA"));
  else if (!(strcmp(atomname, "CB"))) {
    if ((!(strcmp(residue, "SER")))||(!(strcmp(residue, "THR"))))
      return(GetSMoGType("CP"));
    else
      return(GetSMoGType("C3"));
  }
  else if (!(strncmp(atomname, "CG", 2))) {
    if((!(strcmp(residue, "ASP")))||(!(strcmp(residue, "ASN"))))
      return(GetSMoGType("CC"));
    else if((!(strcmp(residue, "PHE")))||(!(strcmp(residue, "TYR")))||
	    (!(strcmp(residue, "HIS")))||(!(strcmp(residue, "TRP"))))
      return(GetSMoGType("C2"));
    else 
      return(GetSMoGType("C3"));
  }
  else if (!(strncmp(atomname, "CD", 2))) {
    if((!(strcmp(residue, "GLU")))||(!(strcmp(residue, "GLN"))))
      return(GetSMoGType("CC"));
    else if (!(strcmp(residue, "ARG")))
      return(GetSMoGType("CP"));
    else if (!(strcmp(residue, "PRO")))
      return(GetSMoGType("CA"));
    else if((!(strcmp(residue, "PHE")))||(!(strcmp(residue, "TYR")))||
	    (!(strcmp(residue, "HIS")))||(!(strcmp(residue, "TRP"))))
      return(GetSMoGType("C2"));
    else
      return(GetSMoGType("C3"));
  }
  else if (!(strncmp(atomname, "CE", 2))) {
    if ((!(strcmp(residue, "LYS")))||(!(strcmp(residue, "HIS"))))
      return(GetSMoGType("CP"));
    else
      return(GetSMoGType("C2"));
  }
  else if (!(strncmp(atomname, "CZ", 2)))   {
    if(!(strcmp(residue, "ARG")))
      return(GetSMoGType("CC"));
    else if (!(strcmp(residue, "TYR")))
      return(GetSMoGType("CP"));
    else
      return(GetSMoGType("C2"));
  }
  else if (!(strncmp(atomname, "CH", 2)))
    return(GetSMoGType("C2"));
  else if (!(strcmp(atomname, "N")))
    return(GetSMoGType("NM"));
  else if (!(strncmp(atomname, "ND", 2))) {
    if(!(strcmp(residue, "HIS")))
      return(GetSMoGType("NC"));
    else
      return(GetSMoGType("ND"));
  }
  else if (!(strncmp(atomname, "NE", 2))) {
    if(!(strcmp(residue, "HIS")))
      return(GetSMoGType("NC"));
    else
      return(GetSMoGType("ND"));
  }		
  else if (!(strcmp(atomname, "NZ")))
    return(GetSMoGType("NC"));
  else if (!(strncmp(atomname, "NH", 2)))
    return(GetSMoGType("NC"));
  else if (!(strcmp(atomname, "O")))
    if (!(strncmp(residue, "HOH", 3)))
      return(GetSMoGType("OD"));
    else
      return(GetSMoGType("OB"));
  else if (!(strncmp(atomname, "OG", 2)))
    return(GetSMoGType("OD"));
  else if (!(strncmp(atomname, "OD", 2))) {
    if(!(strcmp(residue, "ASP")))
      return(GetSMoGType("OC"));
    else
      return(GetSMoGType("OB"));
  }
  else if (!(strncmp(atomname, "OE", 2))) {
    if(!(strcmp(residue, "GLU")))
      return(GetSMoGType("OC"));
    else
      return(GetSMoGType("OB"));
  }
  else if (!(strncmp(atomname, "OH", 2)))
    return(GetSMoGType("OD"));
  else if (!(strncmp(atomname, "S", 1)))
    return(GetSMoGType("SP"));
  else if (!(strncmp(atomname, "OXT", 3)))
    return(GetSMoGType("OC"));
  else if (!(strncmp(atomname, "OCT", 3)))
    return(GetSMoGType("OC"));
  else if (!strncmp(atomname, "ME", 2) || !strncmp(atomname, "ZN", 2) || !strncmp(atomname, "FE", 2) || !strncmp(atomname, "MG", 2) || !strncmp(atomname, "MN", 2) || !strncmp(atomname, "CU", 2) )
    return(GetSMoGType("ME"));
  else {
    return(GetSMoGType("H"));
  }
  
  return -500;
}

int GetSMoGType(char *c) {
  int i;

  for (i=0; i<NSMOGTYPES; i++) 
    if (!strcmp(at_types[i], c))
      return i;
  fprintf(STATUS, "atom type %s not found in the parameter table. Check!\n", c);
  exit(0);
}

int TypeMaxAtoms(char *s, char *res) {
  int i;

  for (i = 0; i < natom_type_list; ++i) {
    if(!strcmp(s, atom_type_list[i].atom_name)) {
      if (!strcmp(atom_type_list[i].res_name, "XXX"))
	return(atom_type_list[i].type_num);
      else if (!strcmp(atom_type_list[i].res_name, res))
	return(atom_type_list[i].type_num);
    }
  }
  return(-500);
}

int TypeAtom(char *s, char *res) {
  
  if (!strncmp(s, "C", 1)) {
    if (!strcmp(s, "C") || (!strcmp(s, "CG") && (!strcmp(res, "ASN") || !strcmp(res, "ASP") || !strcmp(res, "HIS") || !strcmp(res, "PHE") || !strcmp(res, "TYR") || !strcmp(res, "TRP"))) || (!strcmp(s, "CD") && (!strcmp(res, "GLN") || !strcmp(res, "GLU"))) || (!strcmp(s, "CZ") && (!strcmp(res, "ARG") || !strcmp(res, "TYR"))) ||  (!strcmp(s, "CD2") && !strcmp(res, "TRP")) || (!strcmp(s, "CE2") && !strcmp(res, "TRP")))
      return 0;
    else if ((!strcmp(s, "CD1") && (!strcmp(res, "PHE") || !strcmp(res, "TYR"))) || (!strcmp(s, "CD2") && (!strcmp(res, "HIS") || !strcmp(res, "PHE") || !strcmp(res, "TYR"))) || (!strcmp(s, "CZ") && !strcmp(res, "PHE")) || (!strcmp(res, "TRP") && (!strncmp(s, "CH", 2) || !strncmp(s, "CZ", 2) || !strcmp(s, "CE3") || !strcmp(s, "CD1"))))
      return 1;
    else if ((!strcmp(s, "CA") && strcmp(res, "GLY")) || (!strcmp(s, "CB") && (!strcmp(res, "ILE") || !strcmp(res, "THR") || !strcmp(res, "VAL"))) || (!strcmp(s, "CG") && !strcmp(res, "LEU")))
      return 2;
    else if ((!strcmp(s, "CB") && !strcmp(res, "ALA")) || (!strcmp(s, "CD1") && (!strcmp(res, "ILE") || !strcmp(res, "LEU"))) || (!strcmp(s, "CD2") && !strcmp(res, "LEU")) || (!strcmp(s, "CG1") && !strcmp(res, "VAL")) || (!strcmp(s, "CG2") && (!strcmp(res, "ILE") || !strcmp(res, "THR") || !strcmp(res, "VAL"))))
      return 4;
    else
      return 3;
  }
  else if (!strncmp(s, "N", 1)) {
    if ((!strcmp(s, "N") && !strcmp(res, "PRO")) || (!strcmp(s, "NE2") && !strcmp(res, "HIS")))
      return 5;
    else if (!strcmp(s, "NZ"))
      return 8;
    else if (!strcmp(s, "N") || (!strcmp(s, "ND1") && !strcmp(res, "HIS")) || (!strcmp(s, "NE") && !strcmp(res, "ARG")) || (!strcmp(s, "NE1") && !strcmp(res, "TRP")))
      return 6;
    else
      return 7;
  }
  else if (!strncmp(s, "O", 1)) {
    if (!strcmp(s, "OH") || !strcmp(s, "OG") || !strcmp(s, "OG1"))
      return 10;
    else
      return 9;
  }
  else if (!strncmp(s, "S", 1)) {
    if (!strcmp(res, "MET"))
      return 11;
    else
      return 12;
  }
  else
    return -500;
  
}

int IsSidechainAtom(char *atomname) {

  if (!strcmp(atomname, "C") || !strcmp(atomname, "N") || !strcmp(atomname, "O") || !strcmp(atomname, "CA") || !strcmp(atomname, "OCT") || !strcmp(atomname, "OXT"))
    return 0;
  else
    return 1;
}

void ParsePDBLine(char *line, struct atom *protein, int *Natoms) {
  char temp[10];
  short primary_rotamer=1;

  static int previous_res = -1;
  
  if (!strncmp(line, "ATOM", 4)) {
    strncpy(protein[*Natoms].res, &(line[17]), 3);
    strcpy(&(protein[*Natoms].res[3]), "\0");
    strncpy(protein[*Natoms].atomname, &(line[12]), 4);
    strcpy(&(protein[*Natoms].atomname[4]), "\0");
    squeeze(protein[*Natoms].res, ' ');	  
    squeeze(protein[*Natoms].atomname, ' '); 
    strncpy(temp, &(line[22]), 4);
    strcpy(&temp[4], "\0");
    protein[*Natoms].res_num = atoi(temp);
    protein[*Natoms].is_core = IsCoreResidue(protein[*Natoms].res_num);
    protein[*Natoms].is_designed = IsDesignedResidue(protein[*Natoms].res_num);
    strncpy(temp, &(line[30]), 8);
    strcpy(&temp[8], "\0");
    protein[*Natoms].xyz.x=atof(temp);
    strncpy(temp, &(line[38]), 8);
    strcpy(&temp[8], "\0");
    protein[*Natoms].xyz.y=atof(temp);
    strncpy(temp, &(line[46]), 8);
    strcpy(&temp[8], "\0");
    protein[*Natoms].xyz.z=atof(temp);
    protein[*Natoms].is_sidechain = IsSidechainAtom(protein[*Natoms].atomname);
    protein[*Natoms].atomtype = TypeAtom(protein[*Natoms].atomname, protein[*Natoms].res);
    if (!USE_GO_POTENTIAL)
      protein[*Natoms].smogtype = TypeMaxAtoms(protein[*Natoms].atomname, protein[*Natoms].res);
    else 
      protein[*Natoms].smogtype = *Natoms;
    if (!isspace(line[16])) {
      fprintf(STATUS, "WARNING: there is an alternate location indicator:\n");
      fprintf(STATUS, "\t%s\n", line);
      if (line[16]!='A')
	primary_rotamer = 0;
    }
    if ((protein[*Natoms].atomtype!=-500) && (protein[*Natoms].smogtype!=-500) && (USE_SIDECHAINS || !protein[*Natoms].is_sidechain) && (primary_rotamer)){ 
      (*Natoms)++;
      if (*Natoms==1)
	first_atom_res = protein[0].res_num;
      protein[*Natoms-1].res_num-=first_atom_res;

      if (protein[*Natoms-1].res_num != previous_res)
       {
	res_atomno[protein[*Natoms-1].res_num] = *Natoms -1;
	previous_res = protein[*Natoms-1].res_num;
       }
    }
  }

  return;
}

int IsCoreResidue(int res) {

  switch(res) {
  case 3: case 5: case 7: case 13: case 15: case 17: case 21: case 23: case 26: case 30: case 36: case 41: case 43: case 45: case 50: case 56: case 61: case 67: case 69:
    return 1;
  default:
    return 0;
  }
  
  return -1;

}

int IsDesignedResidue(int res) {
  
  switch(res) {
  case 1: case 3: case 4: case 5: case 7: case 13: case 15: case 17: case 21: case 22: case 23: case 26: case 27: case 29: case 30: case 31: case 33: case 34: case 36: case 40: case 41: case 42: case 43: case 44: case 45: case 48: case 49: case 50: case 54: case 55: case 56: case 58: case 59: case 61: case 65: case 67: case 68: case 69: case 71: case 72: 
    return 1;
  default:  
    return 0;
  }
  return -1;

}

int GetAminoNumber(char *name) {
  
  if (!strcmp(name, "ALA"))
    return 0;
  else if (!strcmp(name, "ARG"))
    return 1;
  else if (!strcmp(name, "ASN"))
    return 2;
  else if (!strcmp(name, "ASP"))
    return 3;
  else if (!strcmp(name, "CYS"))
    return 4;
  else if (!strcmp(name, "GLN"))
    return 5;
  else if (!strcmp(name, "GLU"))
    return 6;
  else if (!strcmp(name, "GLY"))
    return 7;
  else if (!strcmp(name, "HIS"))
    return 8;
  else if (!strcmp(name, "ILE"))
    return 9;
  else if (!strcmp(name, "LEU"))
    return 10;
  else if (!strcmp(name, "LYS"))
    return 11;
  else if (!strcmp(name, "MET"))
    return 12;
  else if (!strcmp(name, "PHE"))
    return 13;
  else if (!strcmp(name, "PRO"))
    return 14;
  else if (!strcmp(name, "SER"))
    return 15;
  else if (!strcmp(name, "THR"))
    return 16;
  else if (!strcmp(name, "TRP"))
    return 17;
  else if (!strcmp(name, "TYR"))
    return 18;
  else if (!strcmp(name, "VAL"))
    return 19;
  else
    return 0;

};

