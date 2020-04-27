struct low_resol {
  struct vector xyz;
};

struct atom {
  struct vector xyz;
  struct int_vector xyz_int;
  char res[5];        /* residue name */
  short res_num;        /* residue number */
  char atomname[5];
  short atomtype;
  short smogtype;
  short is_core;
  short is_designed;
  short is_sidechain;
  short sec_structure;
  unsigned char X, Y, Z;
  struct cell *matrix;
};

struct residue {
  char res[4];
  int amino_num;
  short is_core;
  short is_designed;
  Float psi;
  Float phi;
  int ntorsions;
  int nrotamers;
  Float avg_angle[4][4][4][4][4];
  short **rot_position;
  Float chi[4];
  Float native_chi[4];
  Float tmpchi[4];
  int CA, N, C, O, CB, CG, CE1, CE2, CZ2, CZ3;
  int atomnumber[40];
};


void CopyAtom(struct atom, struct atom *);

/* copy from protein1 to protein2 */
void CopyAtom(struct atom protein1, struct atom *protein2) {
  
  (*protein2).X = protein1.X;
  (*protein2).Y = protein1.Y;
  (*protein2).Z = protein1.Z;
  (*protein2).xyz_int.x = protein1.xyz_int.x;
  (*protein2).xyz_int.y = protein1.xyz_int.y;
  (*protein2).xyz_int.z = protein1.xyz_int.z;
  (*protein2).xyz.x = protein1.xyz.x;
  (*protein2).xyz.y = protein1.xyz.y;
  (*protein2).xyz.z = protein1.xyz.z;
  (*protein2).res_num = protein1.res_num;
  (*protein2).atomtype = protein1.atomtype;
  (*protein2).smogtype = protein1.smogtype;	
  (*protein2).is_core = protein1.is_core;
  (*protein2).is_designed = protein1.is_designed;
  (*protein2).is_sidechain = protein1.is_sidechain;
  strcpy((*protein2).res, protein1.res);
  strcpy((*protein2).atomname, protein1.atomname);
  (*protein2).matrix = protein1.matrix;
  (*protein2).sec_structure = protein1.sec_structure;
  return;

}
