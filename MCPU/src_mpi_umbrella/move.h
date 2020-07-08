void BackboneMove(Float);
void LoopBackboneMove(Float);
void LocalBackboneMove(Float);
void SidechainMove();
void MakeMove(Float, Float);



/*==================================================*/
/*                  main program                    */
/*==================================================*/


void MakeMove(Float step_size, Float use_global_bb_moves) {
  
  int reject, del, N, M, i;
  sidechain_step = 0;
  int use_yang = 0;
  int n_soln = 0;
//  float hb_before_move = 0.;
//  float hb_after_move = 0.;
//  int j;
//  float e;
  struct alignment align;

  NFRAG = 1;
  align.seqptr[1].x1 = 1;
  align.structptr[1].x1 = 1;
  align.seqptr[1].x2 = nresidues;
  align.structptr[1].x2 = nresidues;
  backbone_accepted = 0;
  
    /* wmj ************************************* */
  diff_number_of_contacts_current =natives- number_of_contacts_setpoint;
  
  /* end wmj ********************************* */
  

  do {                          /* See the last line of this clause, while (sidechain_step++ < SIDECHAIN_MOVES)
                                 In other words, you make a backbone move first (since sidechain_step above is set to 0
                                 Then you make however many sidechain moves as dictated by variable SIDECHAIN_MOVES
                                 i.e. with SIDECHAIN_MOVES = 1, the loop execs twice...One backbone, one sidechain  */

    reject = 0;
    mc_flags.clashed = 0;
    nomove = 0;
    sidemovedone = 0;
    
    /* make backbone move */
    if (sidechain_step == 0) {

      if ((use_global_bb_moves) && (use_global_bb_moves > 0.001)) {
	fprintf(STATUS, "USE_GLOBAL_BB_MOVES is turned on!!!");
	exit(1);
      }
      
      //Try one of the various possible backbone moves
      if (threefryrand() < CLUSTER_MOVE) { //Cluster_move will often be set to 0 to suppress knowledge moves
        //fprintf(STATUS, "LoopBackbondMove():\n");
        LoopBackboneMove(step_size); //Make a knowledge based backbone move
        use_yang = 0;
      }
      else {
        if (YANG_MOVE && threefryrand() < YANG_MOVE) { //Make a local move, as in Dill paper
	    integloop(step_size, &n_soln); //See loop.h
	    use_yang = 1;
        }
        else {
            //fprintf(STATUS, "LocalBackboneMove():\n");
	    LocalBackboneMove(step_size); //Counterintuitively, this is actually a global move
	    use_yang = 0;
        }
      } /* end (sidechain_step != 0) */

      /* VZ
       * struct_f2, the structure backbone, is updated here
       * because a backbone move was attempted.
       * Acts only if rmsd_constraint is active
       * Note: struct_f1 and struct_f2 are declared in backbone.h,
       *   and struct_f1 is set in fold.h
       */
      if (rmsd_constraint > 0) {
        for (i=0; i<nalign; ++i) {
          struct_f2[i+1].CA.x = native[native_residue[map_to_seq[i]].CA].xyz.x;
          struct_f2[i+1].CA.y = native[native_residue[map_to_seq[i]].CA].xyz.y;
          struct_f2[i+1].CA.z = native[native_residue[map_to_seq[i]].CA].xyz.z;
        }
        new_rms = getrms(struct_f1, struct_f2, align);
      }
    }
    else {    /* make sidechain move (sidechain_step != 0) */
      if (total_ntorsions != 0) {
        //fprintf(STATUS, "SidechainMove():\n");
        SidechainMove();
      }
      use_yang = 0;
    }

    if (!mc_flags.init && mc_flags.clashed) {
      reject = 1;
      if ((nomove == 0) && (sidechain_step == 0))
        nrejected++;
    }
    else if ((use_yang == 1) && (n_soln == 0)) {
      reject = 1;
      //fprintf(STATUS, "Yang move rejected \n");
      if (sidechain_step == 0)
        nothers++;
    }
    //else if ((use_yang == 1) && (n_soln != 0)) {
    //fprintf(STATUS, "Yang move accepted \n");

    //}
    /* below block is new code to support RMSD constraint
     * commenting this block fixes bug.
     * --> has to be something messing with logic?
     */
    else if ((sidechain_step == 0) && (rmsd_constraint > 0) && (new_rms > rmsd_constraint)) {
      reject = 1;
      nrejected++;
    }
    else {

      delta_nclashes = 0;
      if (mc_flags.init)  
	for (i=0; i<total_pairs2; i++) 
	  delta_nclashes += data[cd[i].a][cd[i].b].delta_clashes - data[cd[i].a][cd[i].b].clashes;

      dE = 0; dE_pot = 0; dE_hbond = 0; dE_tor = 0; dE_aro = 0; dE_sct = 0; dE_constraint = 0; 

      if (sidechain_step == 0)
	dE_tor = torsionenergy() - prev_E_tor;

      if ((sidechain_step != 0) && sidemovedone)
        dE_sct = sctenergy() - prev_E_sct;

      dE_aro = aromaticenergy() - prev_E_aro;
      delta_contacts = 0;
      for (i=0; i<total_pairs; i++) { /*AB: Loop through all pairs of atoms and compute potential change due to move*/
	N = ab[i].a; M = ab[i].b;  /*N, M are indices for pairs of atoms*/
	del = data[N][M].delta_contacts-data[N][M].contacts;  /*I think del indicates somethign like how many new contacts have formed?*/
	delta_contacts += del;
	dE_pot += ((Float) del)*potential[native[N].smogtype][native[M].smogtype];  /*Compute change in mu potential contribution ...I think the del in front somehow restricts your attention to atoms that have come into contact?*/
      }
      if (weight_hbond) {
        if (sidechain_step == 0)
	  dE_hbond = FoldHydrogenBonds() - prev_E_hbond;
      }

      /* here, code used to test for rmsd */
//	align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &bb_rms);
      
      
      
      if (strcmp(constraint_file, "None")!=0 ){ //AB
  	dE_constraint =  Compute_constraint_energy (native_residue, native) - prev_E_constraint; //AB
   	} //AB: Otherwise it is 0, as previously initialized
      
      dE = weight_potential*dE_pot + weight_clash*delta_nclashes + weight_hbond*dE_hbond + TOR_WEIGHT*dE_tor + ARO_WEIGHT*dE_aro + SCT_WEIGHT*dE_sct + dE_constraint; //AB addded last term
      double arg = dE / MC_TEMP;
    	/* wmj ****************************** */
    	
      new_natives= number_of_calpha_native_contacts (native_residue, native, orig_contactstring);
      diff_number_of_contacts_new =
	new_natives
	- number_of_contacts_setpoint;
      double dbias = k_bias * (diff_number_of_contacts_new*diff_number_of_contacts_new
			       - diff_number_of_contacts_current*diff_number_of_contacts_current);
      arg += dbias;
      /* end wmj ************************** */
      
      
      /************************  Start AB ***********************************/
      double acceptance_crit = exp(-arg);
      
	  //Jacobi fix implemented 7/17 by AB
      if (use_yang == 1){
      acceptance_crit = acceptance_crit * (n_soln/(double)soln_no_before) * (jacobi_after/jacobi_before); //AB changed to obey detailed balance in case local moves are used
      //fprintf(STATUS, "jacobi_after = %7.3f, jacobi_before = %7.3f, n_solutions_after = %i, n_solutions_before = %i, dE/MC_TEMP = %7.3f, acceptance criterion is %7.5f \n", jacobi_after, jacobi_before, n_soln, soln_no_before,arg, acceptance_crit);
      //if (acceptance_crit>=1 || threefryrand() < acceptance_crit){
      //fprintf(STATUS, "Move may have been accepted \n");
      //};
      }
      /*End AB*/
      
      /* Metropolis */
      
      /*AB commented out: */
      //if (arg <= 0 || threefryrand() < exp(-arg)){
      
      /*AB replaced the above with*/
      if (acceptance_crit >=1 || threefryrand() < acceptance_crit){ //Move is accepted!
	//if (use_yang ==1){
	//fprintf(STATUS, "Move accepted \n");
	//}
	
	//end AB
	Update();
	natives=new_natives;
        //fprintf(STATUS, "Updated ...\n");
        //ResetEnergies();
        if ((nomove == 0) && (sidechain_step == 0)) {
	  naccepted++;
	  backbone_accepted = 1;
        }
        else 
	  n_sidechain_accepted++;
//    fprintf(STATUS, "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", E, 0.2*E_pot+0.8*E_hbond, E_pot, E_hbond, dE, dE_pot, dE_hbond);
      }
      else {
	reject = 1;
        if ((nomove==0) && (sidechain_step == 0))
	  nrejected++;
      }
    }

    if (reject) { 
      Restore();
      //fprintf(STATUS, "rejected ...\n");
      //ResetEnergies();
    }

    for (i=0; i<all_rotated_natoms; i++) {
      is_rotated[all_rotated_atoms[i]] = 0;
    }
    
  } while (sidechain_step++ < SIDECHAIN_MOVES);
  
  return;
  
}

/*==================================================*/


void LocalBackboneMove(Float step_size) {
  int i;

  mc.loop_size = 1;
  all_rotated_natoms = 0;
  total_pairs = total_pairs2 = 0;
  total_hbond_pairs = 0;
  nomove = 0;

  mc.sel_res_num = (int)(threefryrand()*nresidues);
  if(is_template[mc.sel_res_num]==1) {
    nomove = 1;
    nothers++;
    return;
  }
  if (native_residue[mc.sel_res_num].amino_num==14) { /* cannot rotate around proline phi */
    if (mc.sel_res_num != nresidues-1) 
      mc.is_phi = 0; 
    else { /* if proline is last residue, no rotation possible */
      nomove = 1;
      nothers++;
      return;
    }
  }
  else if (mc.sel_res_num == 0)
    mc.is_phi = 0;
  else if (mc.sel_res_num == nresidues-1)
    mc.is_phi = 1;
  else
    mc.is_phi = (int) (threefryrand()*2);
  //fprintf(STATUS, "mainchain move at %5d\n", mc.sel_res_num);

  step_size = step_size*GaussianNum();;

  BackboneMove(step_size);
  UpdateLattice(rotate_natoms[mc.is_phi][mc.sel_res_num], rotate_atom[mc.is_phi][mc.sel_res_num]);

  for (i=0;i<rotate_natoms[mc.is_phi][mc.sel_res_num];i++) {
    is_rotated[rotate_atom[mc.is_phi][mc.sel_res_num][i]]=1;
  }

  NewDeltaContacts(rotate_natoms[mc.is_phi][mc.sel_res_num], rotate_atom[mc.is_phi][mc.sel_res_num], not_rotated[mc.is_phi][mc.sel_res_num]);
  return;

}

void LoopBackboneMove(Float absolute_step_size) {
/*Contrary to the naming, this one actually does the knowledge based move*/

  int a, b, c, d, i, j;
  Float step_phi, step_psi;
  Float desire_phi, desire_psi;
  int cluster_bin;
  float use_cluster;

  /* every 1 moves, a new triplet of residues is selected. */
  /* each move, the phi/psi angles of a single residue are changed */
 
  mc.loop_size = 1; 

  /* triplet is selected */
  if (0 == 0) {
    mc.sel_triplet = (int) (threefryrand()*TOTAL_TRIPLE_LOOP_MOVES)+TOTAL_SINGLE_LOOP_MOVES+TOTAL_DOUBLE_LOOP_MOVES;
 
    if (residue_triplets[mc.sel_triplet].a > nresidues/2.0) {
      mc.selected[0] = residue_triplets[mc.sel_triplet].a;
      mc.selected[1] = residue_triplets[mc.sel_triplet].b;
      mc.selected[2] = residue_triplets[mc.sel_triplet].c;
    }
    else {
      if (residue_triplets[mc.sel_triplet].c>=0) {
	mc.selected[0] = residue_triplets[mc.sel_triplet].c;
	mc.selected[1] = residue_triplets[mc.sel_triplet].b;
	mc.selected[2] = residue_triplets[mc.sel_triplet].a;
      }
      else if (residue_triplets[mc.sel_triplet].b>=0) {
	mc.selected[0] = residue_triplets[mc.sel_triplet].b;
	mc.selected[1] = residue_triplets[mc.sel_triplet].a;
	mc.selected[2] = residue_triplets[mc.sel_triplet].c;
      } 
      else {
	mc.selected[0] = residue_triplets[mc.sel_triplet].a;
	mc.selected[1] = residue_triplets[mc.sel_triplet].b;
	mc.selected[2] = residue_triplets[mc.sel_triplet].c;
      }   
    }
  }
  if(is_template[mc.selected[0]]==1)
   {
    nomove = 1;
    nothers++;
    return;
   }      

  all_rotated_natoms = 0;
  total_pairs=total_pairs2=0;
  total_hbond_pairs=0;

  /* move-cycle defines the residue to be changed */
  i = move_cycle%1;
  move_cycle = i;
//  ++move_cycle;

  check_phipsi();
  cluster_bin = (int)(NOCLUSTERS*threefryrand());
//  mc.selected[i] = 92;
//  cluster_bin = 1;
  desire_phi = cluster_phi[native_residue[mc.selected[i]].amino_num][cluster_bin];
  desire_psi = cluster_psi[native_residue[mc.selected[i]].amino_num][cluster_bin];
//fprintf(STATUS, "%d\n", cluster_bin);
  /* rotates phi and psi */
  if(secstr[mc.selected[i]]=='H')
    use_cluster = 0.5;
  else if(secstr[mc.selected[i]]=='E')
    use_cluster = 0.0;
  else if(secstr[mc.selected[i]]=='L')
    use_cluster = 0.0;
  else if(secstr[mc.selected[i]]=='C')
    use_cluster = USE_CLUSTER;
  else
   {
    fprintf(STATUS, "ERROR! secondary structure prediction has a unknown charactrer: %c\n", secstr[mc.selected[i]]);
    exit(1);
   }
  if (mc.selected[i] > nresidues/2.0) {
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1) && native_residue[mc.selected[i]].amino_num!=14) {
	/* can't rotate around a proline phi */
        if (use_cluster > threefryrand()) {
          step_phi = desire_phi - cur_phi[mc.selected[i]-1];
          step_phi += GaussianNum()*CLUSTER_NOISE;
        //fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
        }
        else
          //fprintf(STATUS, "LoopMove called, but no knowledge based move done\n");
          step_phi = GaussianNum()*CLUSTER_NOISE;
	step_phi *= deg2rad;
//	fprintf(STATUS, "%f\n", step_size);
	
	a = native_residue[mc.selected[i]-1].C;
	b = native_residue[mc.selected[i]].N;
	c = native_residue[mc.selected[i]].CA;
	d = native_residue[mc.selected[i]].C;
        DoRotation(a, b, c, d, step_phi, loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0]);
	mc.delta_phi_angle[i]=step_phi;
        for (j=0;j<loop_rotate_natoms[mc.selected[i]][0];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][0][j]]=1;
        }
      }
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1)) {
        if (use_cluster > threefryrand())
         {
          step_psi = desire_psi - cur_psi[mc.selected[i]-1];
          step_psi += GaussianNum()*CLUSTER_NOISE;
          //fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
          //fprintf(STATUS, "Did a knowledge based move!!\n", secstr[mc.selected[i]]);
         }
        else
          //fprintf(STATUS, "LoopMove called, but no knowledge based move done\n");
          step_psi = GaussianNum()*CLUSTER_NOISE;
	step_psi *= deg2rad;
	
	a = native_residue[mc.selected[i]].N;
	b = native_residue[mc.selected[i]].CA;
	c = native_residue[mc.selected[i]].C;
	d = native_residue[mc.selected[i]+1].N;
       DoRotation(a, b, c, d, step_psi, loop_rotate_natoms[mc.selected[i]][1], loop_rotate_atoms[mc.selected[i]][1]);
	mc.delta_psi_angle[i]=step_psi;
        for (j=0;j<loop_rotate_natoms[mc.selected[i]][1];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][1][j]]=2;
        }
      } 
  }
  else {
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1)) {
        if (use_cluster > threefryrand())
         {
          step_psi = desire_psi - cur_psi[mc.selected[i]-1];
          step_psi += GaussianNum()*CLUSTER_NOISE;
          //fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
         }
        else
          step_psi = GaussianNum()*CLUSTER_NOISE;
	step_psi *= -deg2rad;

	a = native_residue[mc.selected[i]].N;
	c = native_residue[mc.selected[i]].CA;
	b = native_residue[mc.selected[i]].C;
	d = native_residue[mc.selected[i]+1].N;
	DoRotation(a, b, c, d, -step_psi, loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0]);
	mc.delta_psi_angle[i]=-step_psi;
        for (j=0;j<loop_rotate_natoms[mc.selected[i]][0];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][0][j]]=1;
        }
      } 
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1) && native_residue[mc.selected[i]].amino_num!=14) {
        if (use_cluster > threefryrand())
         {
          step_phi = desire_phi - cur_phi[mc.selected[i]-1];
          step_phi += GaussianNum()*CLUSTER_NOISE;
         }
        else
          step_phi = GaussianNum()*CLUSTER_NOISE;
	step_phi *= -deg2rad;
	
	a = native_residue[mc.selected[i]-1].C;
	c = native_residue[mc.selected[i]].N;
	b = native_residue[mc.selected[i]].CA;
	d = native_residue[mc.selected[i]].C;
	DoRotation(a, b, c, d, -step_phi, loop_rotate_natoms[mc.selected[i]][1], loop_rotate_atoms[mc.selected[i]][1]);
	mc.delta_phi_angle[i]=-step_phi;
        for (j=0;j<loop_rotate_natoms[mc.selected[i]][1];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][1][j]]=2;
        }
      }
  }
//  check_phipsi();
//   fprintf(STATUS, "%3d %s %3d %3d %7.5f %7.5f %7.5f %7.5f\n", 
//       mc.selected[i], native_residue[mc.selected[i]].res, native_residue[mc.selected[i]].amino_num, cluster_bin, 
//       desire_phi, cur_phi[mc.selected[i]-1], desire_psi, cur_psi[mc.selected[i]-1]);

  UpdateLattice(loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0]);
  all_rotated_natoms = loop_rotate_natoms[mc.selected[i]][0]; 
  all_rotated_atoms = loop_rotate_atoms[mc.selected[i]][0];
  
  NewDeltaContacts(loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0], loop_not_rotated[mc.selected[i]][0]);
  return;

}
 
void BackboneMove(Float step_size) {
  short temp;
  int a, b, c, d;   /*          d    */
                 /*         /     */ 
                 /*    b - c      */
                 /*   /           */
                 /*  a            */

  mc.delta_angle[0] = step_size;

  if (mc.is_phi) {
    a = native_residue[mc.sel_res_num-1].C;
    b = native_residue[mc.sel_res_num].N;
    c = native_residue[mc.sel_res_num].CA;
    d = native_residue[mc.sel_res_num].C;
    mc.delta_phi_angle[0] = step_size; /* just for proper updating purposes */
    mc.delta_psi_angle[0] = 0;
  }
  else {
    a = native_residue[mc.sel_res_num].N;
    b = native_residue[mc.sel_res_num].CA;
    c = native_residue[mc.sel_res_num].C;
    d = native_residue[mc.sel_res_num+1].N;
    mc.delta_psi_angle[0] = step_size;
    mc.delta_phi_angle[0] = 0;
  } 

  if (mc.sel_res_num <= nresidues/2.0) { 
    mc.delta_angle[0] *=-1;
    mc.delta_psi_angle[0] *=-1;
    mc.delta_phi_angle[0] *=-1;
    temp = b;
    b = c;
    c = temp;
  } 

  DoRotation(a, b, c, d, mc.delta_angle[0], rotate_natoms[mc.is_phi][mc.sel_res_num], rotate_atom[mc.is_phi][mc.sel_res_num]);

  all_rotated_natoms = rotate_natoms[mc.is_phi][mc.sel_res_num];
  all_rotated_atoms = rotate_atom[mc.is_phi][mc.sel_res_num];
  
  return;

}

void MakeSidechainMove() {
  int a, b, c, d, i, j;
  int sel_no_chi;
  float p_0to1;
  float cummul_prob;

  mc.sel_rotamer = (int) (threefryrand()*native_residue[mc.sel_res_num].nrotamers);
  if (USE_ROT_PROB == 1) {
    p_0to1 = threefryrand()*100;
    cummul_prob = 0.;
    for (i=0; i<no_chi_list[native_residue[mc.sel_res_num].amino_num]; i++) {
      cummul_prob += prob_ang[native_residue[mc.sel_res_num].amino_num][i]; 
      if (cummul_prob > p_0to1)
	break;
    }
    sel_no_chi = i;
  }
  else
    sel_no_chi = (int) (threefryrand()*no_chi_list[native_residue[mc.sel_res_num].amino_num]);
  //fprintf(STATUS, "sidechain move at %5d\n", mc.sel_res_num);
  
  old_rotamer = cur_rotamers[mc.sel_res_num];
  cur_rotamers[mc.sel_res_num] = mc.sel_rotamer;

  for (i=0; i<native_residue[mc.sel_res_num].ntorsions; i++) {

    a = sidechain_torsion[mc.sel_res_num][i][0];
    b = sidechain_torsion[mc.sel_res_num][i][1];
    c = sidechain_torsion[mc.sel_res_num][i][2];
    d = sidechain_torsion[mc.sel_res_num][i][3];
//    fprintf(STATUS, "%10ld %3d %s chi: %d  no_list: %2d %8.3f %8.3f %8.3f %8.3f %8.3f\n", mcstep, mc.sel_res_num, native_residue[mc.sel_res_num].res, i, 
//           sel_no_chi, rotamer_angles[mc.sel_res_num].chis[sel_no_chi][i], deviation_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi][i], 
//	   p_0to1, cummul_prob, prob_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi]);

    if (USE_ROTAMERS)
      mc.delta_angle[i] = rotamer_angles[mc.sel_res_num].chis[sel_no_chi][i] +
                     	  GaussianNum()*deviation_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi][i]
                  	  - native_residue[mc.sel_res_num].chi[i]; 
    else
      mc.delta_angle[i] = GaussianNum()*SIDECHAIN_NOISE;
    native_residue[mc.sel_res_num].tmpchi[i] = native_residue[mc.sel_res_num].chi[i] + mc.delta_angle[i];
    
    DoRotation(a, b, c, d, mc.delta_angle[i], rotate_sidechain_natoms[mc.sel_res_num][i], rotate_sidechain_atom[mc.sel_res_num][i]);
    for (j=0;j<rotate_sidechain_natoms[mc.sel_res_num][i];j++) {
      is_rotated[rotate_sidechain_atom[mc.sel_res_num][i][j]]=i+1;
    }
  }

  all_rotated_natoms = rotate_sidechain_natoms[mc.sel_res_num][0]; 
  all_rotated_atoms = rotate_sidechain_atom[mc.sel_res_num][0];
  
  return;
}

void SidechainMove() {

  all_rotated_natoms = 0;
  total_pairs = total_pairs2 = 0;
  total_hbond_pairs = 0;

  do {
    mc.sel_res_num = (int) (threefryrand()*nresidues);
  } while (native_residue[mc.sel_res_num].nrotamers <=1 || native_residue[mc.sel_res_num].amino_num==14);

  MakeSidechainMove();
  sidemovedone = 1;
  UpdateLattice(rotate_sidechain_natoms[mc.sel_res_num][0], rotate_sidechain_atom[mc.sel_res_num][0]);   
  NewDeltaContacts(rotate_sidechain_natoms[mc.sel_res_num][0], rotate_sidechain_atom[mc.sel_res_num][0], sidechain_not_rotated[mc.sel_res_num][0]);
  
  return;
}


