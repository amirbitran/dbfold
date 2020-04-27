
struct amino {
  char name[4];
  char symbol[2];
  char torsion[4][4][4];
  Float avg_angle[4][4][4][4][4];
  int ntorsions;
  int nrotamers;
  int rotate_natoms[4];
  char rotate_atom[4][10][4];
};
