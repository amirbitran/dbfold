void SetupMatrixStuff(void);
void Fold(void);

/* 
 * the_matrix is a 3D array of struct cell
 */
void SetupMatrixStuff(void) {
  int i, j, k;
  struct cell *the_cell;

  for (i=0; i<MATRIX_SIZE; i++)
    for (j=0; j<MATRIX_SIZE; j++)
      for (k=0; k<MATRIX_SIZE; k++) 
	the_matrix[i][j][k].natoms = 0;

  for (i=0; i<natoms; i++) { 
    FindLatticeCoordinates(&native[i]);

    the_cell = &the_matrix[native[i].X][native[i].Y][native[i].Z];
    the_cell->atom_list[(the_cell->natoms)++] = i;

  }
}

void Fold(void) {
  int i;

  int irep;
  int itmp;
  int ntmp;
  float etmp;
  int partner;

  char rmsd_filename[100];
  char temp_filename[100];
  float E_pot_now = 0.0;
  float E_tor_now = 0.0;
  float E_sct_now = 0.0;
  float E_aro_now = 0.0;
  float E_hbond_now = 0.0;
  float E_constraint_now = 0.0; //AB

  /* rmsd by jyang */

  /* the following two moved to backbone.h */
  /* struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE]; */

  /* 
   * struct alignment contains int NFRAG as well as two arrays of
   * struct fragments of size MAXFRAG (=30).  
   * struct fragments contains two ints: x1, x2
   */
  struct alignment align;

  NFRAG = 1;
  align.seqptr[1].x1 = 1;
  align.structptr[1].x1 = 1;
  align.seqptr[1].x2 = nresidues;
  align.structptr[1].x2 = nresidues;

  mc_flags.init = !NO_NEW_CLASHES;

  /* Save protein data into orig_native */
  for (i=0; i<natoms; i++)
    CopyAtom(native[i], &orig_native[i]);

  CenterProtein(&native, natoms);  
  SetupMatrixStuff();

  Contacts();  
  fprintf(STATUS, "INITIAL CLASHES\t%d\n", nclashes); 
  if (!mc_flags.init) {
    fprintf(STATUS, "\nturning off clashes:\n");
    TurnOffNativeClashes(1);
  }
  fprintf(STATUS, "\n");
  
  
    /* wmj */
  { /* Initialize contactstring */
    orig_contactstring = malloc (sizeof *orig_contactstring * nresidues * nresidues);
    fill_calpha_contact_string (native_residue, struct_native, orig_contactstring);  //AB made this change so that this reads structure_file instead of native_file...still ok to pass native_residue structure so long as structure file has same # of atoms as native file
    //fill_calpha_contact_string (native_residue, orig_native, orig_contactstring);  //can potentially pass struct_native instead of orig_native, while still passing native_residue, since the function just uses native_residue to figure out where each atom is located (at what index)...this should not change so long as I don't change the protein size
  }
  

  /*AB*/
   //Read_substructures(substructure_path, substructures);
   /*End AB*/
    

  /* setup initial energies */

  //align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &native_rms);

  /* 
   * native (struct_f2) is the structure under MC simulation
   * struct_native (struct_f1) is the reference structure for RMSD
   */
  for (i=0; i<nalign; ++i) {
    struct_f1[i+1].CA.x = struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;  //struct_native is our array for the STRUCTURE FILE: that is the "goal" structure that the protein aims to look like
    struct_f1[i+1].CA.y = struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
    struct_f1[i+1].CA.z = struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
    struct_f2[i+1].CA.x = native[native_residue[map_to_seq[i]].CA].xyz.x;
    struct_f2[i+1].CA.y = native[native_residue[map_to_seq[i]].CA].xyz.y;
    struct_f2[i+1].CA.z = native[native_residue[map_to_seq[i]].CA].xyz.z;
  }
  native_rms = getrms(struct_f1, struct_f2, align);

  ResetEnergies(0); //Here we actually fill in the values for prev_E (which are, ironically, the current energies, but these will indeed become previous once a move is made)
  rms_RMSDmin = 100.;
  mcstep_RMSDmin = 0;
  mcstep_Emin = 0;
  Emin = prev_E;
  Emin_pot = prev_E_pot;
  Emin_hbond = prev_E_hbond;
  Emin_tor = prev_E_tor;
  Emin_sct = prev_E_sct;
  Emin_aro = prev_E_aro;
  Emin_constraint = prev_E_constraint;
  fprintf(STATUS, "ENERGY = %.2f(clashes) + %.2f(rmsd) + %.2f(potential) + %.2f(Aro) + %.2f(hbond) + %.2f(tor) + %.2f(sct)\n\n",
    weight_clash, weight_rms, weight_potential, ARO_WEIGHT, weight_hbond, TOR_WEIGHT, SCT_WEIGHT);
  for (i=0; i<natoms; i++) {
    prev_native[i] = native[i];
    native_Emin[i] = prev_native[i];
  }
  for (i=0; i < nresidues; ++i)
    cur_rotamers[i] = 0;

  naccepted = 0;
  nrejected = 0;
  nothers = 0;
  n_sidechain_accepted = 0;

  natives=number_of_calpha_native_contacts (native_residue, native, orig_contactstring);  
  fprintf(STATUS, "         step #    energy contact rmsd   natives potnl    sctor      hbond       Aro   torsion accept reject others   temp setpoint E_constraint weighted_mean_constraint_distance\n---------------------------------------------------------------------------------------------\n");


  /* main folding loop */
  for (mcstep=0; mcstep<MC_STEPS; mcstep++) {
  
  	/* AB: Turn off knowledge based moves at step MAX_CLUSTERSTEP */
  	if (mcstep ==MAX_CLUSTERSTEP){
  		USE_CLUSTER = 0;
  		fprintf(STATUS, "At MC step %10ld, USE_CLUSTER was set to %5.3f \n", mcstep, USE_CLUSTER);
  	}
  	else if (mcstep == 0){
  		fprintf(STATUS, "At MC step %10ld, USE_CLUSTER has value %5.3f \n", mcstep, USE_CLUSTER);
  	}
  	
  	/*Print test by AB */
  	//if (mcstep ==5) {
  	//	fprintf(STATUS, "At MC step %10ld, USE_CLUSTER is set to %5.3f \n", mcstep, USE_CLUSTER);
  	//}
  	//else if (mcstep == 520000) {
  	//	fprintf(STATUS, "At MC step %10ld, USE_CLUSTER is set to %5.3f \n", mcstep, USE_CLUSTER);
  	//}

    /******** replica exchange ****************/
    if (mcstep % MC_REPLICA_STEPS == 0) {
		//ierr=MPI_Barrier(mpi_world_comm);  //Added by AB to ensure synchronization--seems like things were getting out of sync for some reason?
        for (i=0; i<nprocs; i++)
          replica_index[i] = i;
        ierr = MPI_Allgather(&E, 1, MPI_FLOAT, Enode, 1, MPI_FLOAT, mpi_world_comm);
        ierr = MPI_Allgather(&natives, 1, MPI_INT, Nnode, 1, MPI_INT, mpi_world_comm);

        if (myrank == 0) {  //node 0 will mediate the exchange
          for (irep=0; irep<MAX_EXCHANGE; irep++) {
            sel_num = (int) (threefryrand()*(nprocs-1)); /*Who initiates exchange?*/
            /*AB: This was previously nprocs-2, but I saw no reason why second to last node shouldn't be able to initiate exchange
            When we say (int)*random*(n_procs-1), we keep in mind that (int) applies a floor function, so this will draw valeus between 0 and nprocs-2 
            Moreover nprocs-2 corresponds to the second ot last node due to zero indexing
            
            Previously we were only drawing up to third to last node */
            
            //Now, we need to choose the "partner" with whom sel_num exchanges
            if (sel_num%NODES_PER_TEMP == NODES_PER_TEMP -1 ){ //the last setpoint for a given temperature
            /*In this case, the exchange is initiated by the  lowest setpoint node in a given temp
            It wouldn't't make sense for this node to exchange with the node above it in index, since this would 
            involve exchanging with  a partner who has a very different setpoint than the initiator--the acceptance ratio woudl be very low
            
            Thus, we always exchange with the node that is one temperature step above it (i.e. whose rank is NODES_PER_TEMP higher)
            */
            	partner	=sel_num+NODES_PER_TEMP;
            }
            else if (sel_num>=nprocs-NODES_PER_TEMP){
            /*This is saying that the exchange is initiated by somebody with the highest temperature, so obviously can't go up in temperature
            */
            	partner=sel_num+1;
            }
            else {  //we either exchnage with the node that is one temperature above or one contact setpoint below, with 50% probability each
            	if (threefryrand() <0.5) {
            		partner=sel_num+1;
            	}
            	else {
            		partner=sel_num+NODES_PER_TEMP;
            	}
            }

            //fprintf(STATUS, "irep : %d, sel_num : %d\n", irep, sel_num);
            //fflush(STATUS);

            delta_E = Enode[partner]-Enode[sel_num];
            delta_T = 1.0/Tnode[partner]-1.0/Tnode[sel_num];
            delta_N = k_bias*((Nnode[partner]-Cnode[sel_num])*(Nnode[partner]-Cnode[sel_num]) - (Nnode[sel_num]-Cnode[sel_num])*(Nnode[sel_num]-Cnode[sel_num]) ) + k_bias *(  (Nnode[sel_num]-Cnode[partner])*(Nnode[sel_num]-Cnode[partner]) - (Nnode[partner]-Cnode[partner])*(Nnode[partner]-Cnode[partner]) );
            delta_all = delta_E * delta_T - delta_N;
            
            if (delta_all >= 0 || threefryrand() < expf(delta_all)) {
              itmp = replica_index[sel_num];
              replica_index[sel_num] = replica_index[partner];
              replica_index[partner] = itmp;  //for instance, say 5 initiates exchange with 6...then replica_index[5] now gets a value of 6, and vice versa
              etmp = Enode[sel_num];
              Enode[sel_num] = Enode[partner];
              Enode[partner] = etmp;
              ntmp = Nnode[sel_num];
              Nnode[sel_num] = Nnode[partner];
              Nnode[partner] = ntmp;  
              //accepted_replica[sel_num][partner]++;            
              accepted_replica[sel_num]++;
              //fprintf(STATUS, "Node %d successfully exchanged with Node %d, with a delta_all of %8.3f \n", sel_num, partner, delta_all);
            }
            else {
              //rejected_replica[sel_num][partner]++; 
              rejected_replica[sel_num]++;
              //fprintf(STATUS, "Node %d UNSUCCESSFULLY exchanged with Node %d, with a delta_all of %8.3f \n", sel_num, partner, delta_all);
            }
          }
        }

        ierr = MPI_Bcast(replica_index, nprocs, MPI_INT, 0, mpi_world_comm);
        ierr = MPI_Bcast(Enode, nprocs, MPI_FLOAT, 0, mpi_world_comm);
        ierr = MPI_Bcast(Nnode, nprocs, MPI_INT, 0, mpi_world_comm);
        ierr = MPI_Bcast(accepted_replica, nprocs, MPI_INT, 0, mpi_world_comm);
        ierr = MPI_Bcast(rejected_replica, nprocs, MPI_INT, 0, mpi_world_comm);

        //fprintf(STATUS, "Replica Index\n");
        //for (i=0; i<nprocs; i++) fprintf(STATUS, "%5d %5d %8.3f\n", i, replica_index[i], Enode[i]);

		if (mcstep % MC_PRINT_STEPS == 0) {
			fprintf(STATUS, "RPLC %10ld E : %8.3f Natives : %d FROM %2d(T=%5.3f,setpoint=%d ) E : %8.3f, Natives  :  %d, accepted : %5d, rejected : %5d  \n", mcstep, E, natives, replica_index[myrank], Tnode[replica_index[myrank]], Cnode[replica_index[myrank]],Enode[myrank], Nnode[myrank], accepted_replica[myrank], rejected_replica[myrank]);
        	//fprintf(STATUS, "RPLC %10ld E : %8.3f Natives : %d FROM %2d(T=%5.3f,setpoint=%d ) E : %8.3f, Natives  :  %d, accepted : %5d, rejected : %5d  \n", mcstep, E, natives, replica_index[myrank], Tnode[replica_index[myrank]], Cnode[replica_index[myrank]],Enode[myrank], Nnode[myrank], accepted_replica[myrank][replica_index[myrank]], rejected_replica[myrank][replica_index[myrank]]);
        }
        fflush(STATUS);

        for (i=0; i<natoms; i++) {   //Temporary arrays to store all atom coordinates as they stood before any exchanges happened...these will be transferred later
          buf_out[3*i] = native[i].xyz.x;
          buf_out[3*i+1] = native[i].xyz.y;
          buf_out[3*i+2] = native[i].xyz.z;
        }

        for (i=0; i<nprocs; i++){ //normally, replica_index[i] should equal i, unless an exchange occurred
          if (replica_index[i] != i) {  //there is a discrepancy, indicating an exchange occurred...for instance, replica_index[5]=6 and replica_index[6]=5
            if (myrank == i) {  //For instance, if my rank is 6 and replica_index[6]=5, then I received an exchange from node 5
              ierr = MPI_Recv(buf_in, 3*natoms, MPI_FLOAT, replica_index[i], (i+2), mpi_world_comm, &mpi_status);  //now, receive the atomic coordinates from my partner
            } else if (myrank == replica_index[i]) {   //For instance, if my rank is 5 and replica_index[6]=5, then this means I gave an exchange to node 6
              ierr = MPI_Send(buf_out, 3*natoms, MPI_FLOAT, i, (i+2), mpi_world_comm);  //give my coordinates to my partner 
            }
          }
        }

        if (replica_index[myrank] != myrank) { //if exchange occurred, transfer temporary data to real-deal atomic data structure, and recompute energies with the updated atomic info
          for (i=0; i<natoms; i++) {
            native[i].xyz.x = buf_in[3*i];
            native[i].xyz.y = buf_in[3*i+1];
            native[i].xyz.z = buf_in[3*i+2];
          }
        }
      
      /* Re-center */
//      if (!USE_ROTAMERS) {
	CenterProtein(&native, natoms);
	SetupMatrixStuff();
	for (i=0; i<natoms; i++)
	  prev_native[i] = native[i];
//      }
      /* Re-set values that are susceptible to round-off error */ 
      Contacts();
      if (!mc_flags.init)
	TurnOffNativeClashes(0);
      ResetEnergies(0); //Recompute energies
      GetChi();
    
    natives= number_of_calpha_native_contacts (native_residue, native, orig_contactstring);
    } /************ END: replica exchange ********************/

    /* if backbone move was made and accepted we update native_rms */
    if (backbone_accepted == 1) {
      for (i=0; i<nalign; ++i) {
        /* commented out because struct_native doesn't get altered */
        struct_f1[i+1].CA.x = struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;
        struct_f1[i+1].CA.y = struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
        struct_f1[i+1].CA.z = struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
        struct_f2[i+1].CA.x = native[native_residue[map_to_seq[i]].CA].xyz.x;
        struct_f2[i+1].CA.y = native[native_residue[map_to_seq[i]].CA].xyz.y;
        struct_f2[i+1].CA.z = native[native_residue[map_to_seq[i]].CA].xyz.z;
      }
      native_rms = getrms(struct_f1, struct_f2, align);
    }
      
    /* Print Output */
    if (mcstep % MC_PRINT_STEPS == 0) {
      //align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &native_rms);
      /* rmsd by jyang */
      for (i=0; i<nalign; ++i) {
        /* commented out because struct_native doesn't get altered */
        struct_f1[i+1].CA.x = struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;
        struct_f1[i+1].CA.y = struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
        struct_f1[i+1].CA.z = struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
        struct_f2[i+1].CA.x = native[native_residue[map_to_seq[i]].CA].xyz.x;
        struct_f2[i+1].CA.y = native[native_residue[map_to_seq[i]].CA].xyz.y;
        struct_f2[i+1].CA.z = native[native_residue[map_to_seq[i]].CA].xyz.z;
      }
      native_rms = getrms(struct_f1, struct_f2, align);
      TypeContacts();
      E_pot_now = FullAtomEnergy();
      E_sct_now = sctenergy();
      E_tor_now = torsionenergy();
      E_aro_now = aromaticenergy();
      E_hbond_now = HydrogenBonds();
      
        	 if (strcmp(constraint_file, "None")!=0 ){ //AB
 	 E_constraint_now =  Compute_constraint_energy (native_residue, native); //AB
  		}	else { //AB
 	 E_constraint_now =0; //AB
   } //AB

      
      
      
     fprintf(STATUS,"STEP %10ld  %8.2f %6d %5.2f %d %9.2f %9.2f %9.2f    %6.2f %8.2f %6.2f %6.2f %6.2f %6.3f, %d %6.3f %6.3f \n", 
	      mcstep, E, ncontacts, native_rms,
	      natives,
	      E_pot_now, E_sct_now, E_hbond_now, E_aro_now, E_tor_now,
	      100*(Float)naccepted/(Float)MC_PRINT_STEPS,
	      100*(Float)nrejected/(Float)MC_PRINT_STEPS,
	      100*(Float)nothers/(Float)MC_PRINT_STEPS,
	      MC_TEMP, number_of_contacts_setpoint, E_constraint_now, mean_constraint_distance); 
      fflush(STATUS);
      n_sidechain_accepted = 0;
      naccepted = 0;
      nrejected = 0;
      nothers = 0;
    } 

    if (mcstep % 10000 == 0) {
      /* Re-center */
      CenterProtein(&native, natoms);
      SetupMatrixStuff();
      for (i=0; i<natoms; i++)
        prev_native[i] = native[i];

      /* Re-set values that are susceptible to round-off error */
      Contacts();
      if (!mc_flags.init)
	TurnOffNativeClashes(0);
      ResetEnergies(mcstep);
      GetChi();
    }
    
    /* Output Structure */
    if (PRINT_PDB) {
      if (mcstep % MC_PDB_PRINT_STEPS == 0) {
	sprintf(temp_filename, "%s.%ld", pdb_out_file, mcstep);  
    //    if ((MC_TEMP < 0.18) || (mcstep > MC_STEPS - 100000))
	  PrintPDB(temp_filename);
      }
    }

    if (E < Emin) {
      mcstep_Emin = mcstep;
      Emin = E;
      Emin_pot = E_pot;
      Emin_hbond = E_hbond;
      Emin_tor = E_tor;
      Emin_aro = E_aro;
      Emin_constraint = E_constraint; //AB
      rms_Emin = native_rms;
      for (i=0; i<natoms; i++)
	native_Emin[i] = native[i];
    }

    if (native_rms < rms_RMSDmin) {
      mcstep_RMSDmin = mcstep;
      rms_RMSDmin = native_rms;
      E_RMSDmin = E;
      E_RMSDmin_pot = E_pot;
      E_RMSDmin_hbond = E_hbond;
      E_RMSDmin_tor = E_tor;
      E_RMSDmin_sct = E_sct;
      E_RMSDmin_aro = E_aro;
      E_RMSDmin_constraint = E_constraint; //AB
      for (i=0; i<natoms; i++)
	native_RMSDmin[i] = native[i];
    }
    /* Make a move */
    MakeMove(STEP_SIZE, USE_GLOBAL_BB_MOVES);
  } /* end main folding loop */

  sprintf(temp_filename, "%s_Emin.pdb", pdb_out_file);
  PrintPDB_Emin(temp_filename);
  fprintf(STATUS, "\nMC step at Emin: %10ld\n", mcstep_Emin);
  fprintf(STATUS, "Emin:%8.2f  Emin_pot:%8.2f  E_hb:%7.2f  E_tor:%7.2f  E_sct:%7.2f E_aro:%7.2f E_constraint:%7.2f \n",
          Emin, Emin_pot,  Emin_hbond, Emin_tor, Emin_sct, Emin_aro, Emin_constraint);
  fprintf(STATUS, "rmsd at Emin:%8.2f  ", rms_Emin);
  fprintf(STATUS, "Pdb file at Emin: %s\n", temp_filename);

  sprintf(rmsd_filename, "%s_RMSDmin.pdb", pdb_out_file);
  PrintPDB_RMSDmin(rmsd_filename);
  fprintf(STATUS, "\nMC step at RMSDmin: %10ld\n", mcstep_RMSDmin);
  fprintf(STATUS, "rms_Rmin:%8.2f  E_RMSDmin:%8.2f E_Rmin_pot:%8.2f E_Rhb:%8.2f E_Rtor:%8.2f E_Rsct:%8.2f E_Raro:%8.2f E_Rconstraint:%8.2f \n",
         rms_RMSDmin, E_RMSDmin, E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro, E_RMSDmin_constraint);
  fprintf(STATUS, "Pdb file at RMSDmin: %s\n", rmsd_filename);
  
  
  

  //The following is added by AB
  //if (myrank<nprocs-1) {
  	//fprintf(STATUS, "Final exchange statistics with node %i (T = %5.3f, setpoint=%d):  accepted = %i, rejected = %i \n", myrank+1, Tnode[myrank+1], Cnode[myrank+1], accepted_replica[myrank][myrank + 1], rejected_replica[myrank][myrank + 1]);
  	//if (myrank < nprocs- NODES_PER_TEMP) {
  	//	fprintf(STATUS, "Final exchange statistics with node %i (T = %5.3f, setpoint=%d):  accepted = %i, rejected = %i \n", myrank+NODES_PER_TEMP, Tnode[myrank+NODES_PER_TEMP], Cnode[myrank+NODES_PER_TEMP], accepted_replica[myrank][myrank + NODES_PER_TEMP], rejected_replica[myrank][myrank + NODES_PER_TEMP]);
  		
  	//}
  //}
  

  return;
}
