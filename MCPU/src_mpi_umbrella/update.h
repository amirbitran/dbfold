void Restore();
void Update();

void Update() {  /* This function is called AFTER a move has been made. Counterintuitively, the 
values for prev_E are set as the CURRETNT values of the energy. But that's because we make another move before
this function is called again, and so after that other move is made, the prev_E valious really are previous

On the first step of the simulation, the function Reset_energies is used to compute the initial energy, and the values for prev_E are set to that initial energy
*/
  int i, j, k, temp1;
  short M, N, temp;

  for (i=0; i<all_rotated_natoms; i++) {
    N = all_rotated_atoms[i];
    temp_atom = &native[N];
    temp_prev_atom = &prev_native[N];
    // identify matrix index j of the rotated atom N in the previous step
    j=0;
    while(N!=temp_prev_atom->matrix->atom_list[j])
      j++;
    // remove the rotated atom N from the previous matrix
    for (k=j; k<(temp_prev_atom->matrix->natoms-1); k++)
      temp_prev_atom->matrix->atom_list[k] = temp_prev_atom->matrix->atom_list[k+1];
    temp_prev_atom->matrix->natoms--;
    
    if (temp_atom->matrix->natoms) {
      j=0;
      while (j < temp_atom->matrix->natoms && N > temp_atom->matrix->atom_list[j])
	j++;
      if (j!=temp_atom->matrix->natoms) {
	temp = temp_atom->matrix->atom_list[j];
	temp_atom->matrix->atom_list[j] = N;
	for (k=j+1; k< (temp_atom->matrix->natoms); k++) {
	  temp1 = temp_atom->matrix->atom_list[k]; 
	  temp_atom->matrix->atom_list[k] = temp; 
	  temp = temp1;
	}
	temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = temp;
      }
      else
	temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = N;

      if(temp_atom->matrix->natoms>=MAX_CELL_ATOMS){ //We have some sort of problem if the number of atoms in the matrix's atom_list exceeds or is equal to MAX_CELL_ATOMS (100)
        fprintf(STATUS, "Lattice Error: Update(), num. of atoms exceeds maximum num. %d of the cell", MAX_CELL_ATOMS);
        fprintf(STATUS, "Lattice Error: atom %4s %4d %4s\n", temp_atom->atomname, temp_atom->res_num, temp_atom->res);
        fprintf(STATUS, "Lattice Error: atom lists\n");
        for (k=0;k<temp_atom->matrix->natoms;k++){ //Loop through all entries in the atom_list for the matrix affiliated with the current atom
          M=temp_atom->matrix->atom_list[k];  //The value for kth entry in atom list, indexes some other atom in the native structure
          fprintf(STATUS, "Lattice Error: %4d %4s %4d %4s\n", M, native[M].atomname, native[M].res_num, native[M].res);
        }
        exit(1);
      }
    }
    else
      temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = N;

    Copy(temp_atom->xyz, &temp_prev_atom->xyz);
    CopyLatticeCoordinates(native[N], &prev_native[N]);     
  }

  E_pot += dE_pot;
  E_tor += dE_tor;
  E_sct += dE_sct;
  E_aro += dE_aro;
  E_hbond += dE_hbond;
  E_constraint+=dE_constraint;
  E += dE;
  prev_E = E;
  prev_E_pot = E_pot;
  prev_E_tor = E_tor;
  prev_E_sct = E_sct;
  prev_E_aro = E_aro;
  prev_E_hbond = E_hbond;
  prev_E_constraint = E_constraint;
  //fprintf(stdout, "Update : dE_sct : %.5f, E :%.5f, dE : %.5f, prev_E : %.5f\n", dE_sct, E, dE, prev_E);
  //fprintf(stdout, "%9.5f %9.5f\n", E, dE);

  ncontacts+=delta_contacts;

//  if ((sidechain_step!=0)||(USE_ROTAMERS))
  if (sidechain_step!=0)
    for (i=0; i<native_residue[mc.sel_res_num].ntorsions; i++) 
      native_residue[mc.sel_res_num].chi[i] += mc.delta_angle[i];
  else
    for (i=0; i<mc.loop_size; i++) {
      native_residue[mc.selected[i]].phi += mc.delta_phi_angle[i];
      native_residue[mc.selected[i]].psi += mc.delta_psi_angle[i];
    }

  for (i=0; i<total_pairs; i++) {
    M = ab[i].a;
    N = ab[i].b;
    data[M][N].contacts = data[N][M].contacts = data[M][N].delta_contacts;
    data[M][N].delta_contacts=0;
    if (!mc_flags.init)
      data[M][N].delta_clashes=0;
    
  }

  if (!mc_flags.init)
   {
    for (i=0; i< total_hbond_pairs; i++) {
      M = hbond_pair[i].a;
      N = hbond_pair[i].b;
      data[M][N].hbond = data[N][M].hbond = data[M][N].delta_hbond;
      data[M][N].closehb = data[N][M].closehb = data[M][N].delta_closehb;
      data[M][N].delta_hbond=data[N][M].delta_hbond=NO_HBOND;
      data[M][N].delta_closehb=data[N][M].delta_closehb=3;
    }
   }
  
  if (mc_flags.init) {
    nclashes+=delta_nclashes;
    for (i=0; i<total_pairs2; i++) {
      M = cd[i].a;
      N = cd[i].b;
      data[M][N].clashes=data[N][M].clashes=data[M][N].delta_clashes;
      data[M][N].delta_clashes=0;
    }
  }
  
  return;  
}

void Restore() {
  int i;
  short M, N;

  //fprintf(stdout, "Restore(): mc.sel_res_num %d, sidechain_step %1d, sidemovedone %1d\n", mc.sel_res_num, sidechain_step, sidemovedone);

  if ((sidechain_step != 0) && (sidemovedone != 0))
    for (i=0; i<native_residue[mc.sel_res_num].ntorsions; i++)
      native_residue[mc.sel_res_num].tmpchi[i] -= mc.delta_angle[i];

  for (i=0; i<all_rotated_natoms; i++) {
    N = all_rotated_atoms[i];
    temp_atom = &native[N];
    temp_prev_atom = &prev_native[N];
    Copy(temp_prev_atom->xyz, &temp_atom->xyz);
    CopyLatticeCoordinates(prev_native[N], &native[N]);
  }

  if (USE_ROTAMERS)
    cur_rotamers[mc.sel_res_num] = old_rotamer;
 
  for (i=0; i<total_pairs; i++) {
    M = ab[i].a;
    N = ab[i].b;
    data[M][N].delta_contacts = 0;
    data[M][N].delta_clashes = 0;
  }
  
  if (mc_flags.init) {
    for (i=0; i<total_pairs2; i++) {
      M = cd[i].a;
      N = cd[i].b;
      data[M][N].delta_clashes=0;
    }
  }

  else {
    for (i=0; i<total_hbond_pairs; i++){ 
      M = hbond_pair[i].a;
      N = hbond_pair[i].b;
      data[M][N].delta_hbond = data[N][M].delta_hbond = NO_HBOND;
      data[M][N].delta_closehb = data[N][M].delta_closehb = 3;
    }
  }

  E = prev_E;
  E_pot = prev_E_pot;
  E_tor = prev_E_tor;
  E_sct = prev_E_sct;
  E_aro = prev_E_aro;
  E_hbond = prev_E_hbond;
  E_constraint = prev_E_constraint;

  return;
}

