/* MATRIX_SIZE is set to be 20 in init.h (it's a u. char, not a macro) */

struct cell {
  char X, Y, Z;
  struct cell *neighbors[27];
  short atom_list[MAX_CELL_ATOMS];  //AB: Each cell has this array called atom_list which has MAX_CELL_ATOMS = 100 elements
  unsigned char natoms;
};

void FindLatticeCoordinates(struct atom *);
void UpdateLattice(short, short *);
void InitializeMatrix();
unsigned char PerBound(signed char);
void CopyLatticeCoordinates(struct atom, struct atom *);

void CopyLatticeCoordinates(struct atom A, struct atom *B) {
  
  (*B).X = A.X;
  (*B).Y = A.Y;
  (*B).Z = A.Z;
  (*B).xyz_int.x = A.xyz_int.x;
  (*B).xyz_int.y = A.xyz_int.y;
  (*B).xyz_int.z = A.xyz_int.z;
  (*B).matrix = A.matrix;

  return;

}

void UpdateLattice(short rotate_natoms, short *rotate_atom) {
  int j;

  for (j=0; j<rotate_natoms; j++) {
    temp_atom = &native[rotate_atom[j]];
    temp_atom->X = PerBound((int)(floor(temp_atom->xyz.x*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
    temp_atom->Y = PerBound((int)(floor(temp_atom->xyz.y*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
    temp_atom->Z = PerBound((int)(floor(temp_atom->xyz.z*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
    temp_atom->matrix = &(the_matrix[temp_atom->X][temp_atom->Y][temp_atom->Z]);
    temp_atom->xyz_int.x = (long int) (temp_atom->xyz.x*INT_PRECISION);
    temp_atom->xyz_int.y = (long int) (temp_atom->xyz.y*INT_PRECISION);
    temp_atom->xyz_int.z = (long int) (temp_atom->xyz.z*INT_PRECISION);
  }
 
  return;
  
}

  /* computes which lattice cell the given atom falls into */
  /* also, finds the integer representation of the coordinates */
void FindLatticeCoordinates (struct atom *ATOM) {

  (*ATOM).X = PerBound((int)(floor((*ATOM).xyz.x*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
  (*ATOM).Y = PerBound((int)(floor((*ATOM).xyz.y*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
  (*ATOM).Z = PerBound((int)(floor((*ATOM).xyz.z*LATTICE_SIZE)) + HALF_MATRIX_SIZE);
  (*ATOM).xyz_int.x = (long int) ((*ATOM).xyz.x*INT_PRECISION);
  (*ATOM).xyz_int.y = (long int) ((*ATOM).xyz.y*INT_PRECISION);
  (*ATOM).xyz_int.z = (long int) ((*ATOM).xyz.z*INT_PRECISION);
  (*ATOM).matrix = &(the_matrix[(*ATOM).X][(*ATOM).Y][(*ATOM).Z]);
  
  return;
  
}

void InitializeMatrix() {

/*The Matrix is a  3D array, each dimension is given by MATRIX_SIZE = 20 (# of AAs)*/
  short i, j, k, a, b, c;

  the_matrix = (struct cell ***) calloc(MATRIX_SIZE, sizeof(struct cell **));
  for (i=0; i<MATRIX_SIZE; i++) {
    the_matrix[i] = (struct cell **) calloc(MATRIX_SIZE, sizeof(struct cell *));
    for (j=0; j<MATRIX_SIZE; j++)
      the_matrix[i][j] = (struct cell *) calloc(MATRIX_SIZE, sizeof(struct cell));
  }  
  for (i=0; i<MATRIX_SIZE; i++)
    for (j=0; j<MATRIX_SIZE; j++)
      for (k=0; k<MATRIX_SIZE; k++) {
	the_matrix[i][j][k].natoms = 0;
	the_matrix[i][j][k].X = i;
	the_matrix[i][j][k].Y = j;
	the_matrix[i][j][k].Z = k;
	for (a=-1; a<2; a++)
	  for (b=-1; b<2; b++)
	    for (c=-1; c<2; c++) {
	      the_matrix[i][j][k].neighbors[(a+1)*9+(b+1)*3+(c+1)] = &(the_matrix[PerBound(i+a)][PerBound(j+b)][PerBound(k+c)]);
	    }
      }

  return;

}

unsigned char PerBound(signed char X) {
  
  signed char value;
  
  if (X < 0) 
   {
    value = X+MATRIX_SIZE;
    while (value < 0)
      value += MATRIX_SIZE;
    return value;
   }
  else if (X >=MATRIX_SIZE) 
   {
    value = X-MATRIX_SIZE;
    while (value >= MATRIX_SIZE)
      value -= MATRIX_SIZE;
    return value;
   }
  else
    return X;
  
}
