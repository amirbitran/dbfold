#if DEBUG
void Debug();
void DebugContacts();
void CheckForDebugContacts(int, int);

void Debug() {
  int i, j, k;
  short temp;

  DebugContacts();
  for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++)
      if (data[i][j].clashes!=debug_clashes[i][j])
	fprintf(STATUS, "CLASH_MISMATCH\t(%d, %d)\t%d\t%d\n", i, j, data[i][j].clashes, debug_clashes[i][j]);
  for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++)
      if (data[i][j].contacts !=debug_contacts[i][j]) {
	fprintf(STATUS, "CONTACT_MISMATCH\t(%d, %d)\t%d\t%d", i, j, data[i][j].contacts, debug_contacts[i][j]);
	for (k=0; k<total_pairs; k++)
	  if (ab[k].a == i && ab[k].b == j) {
	    fprintf(STATUS, "%d\n", k);
	    break;
	  }
	if (k==total_pairs)
	  fprintf(STATUS, "ERROR\n");
      }
  temp=0;
   for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++)
      temp+=data[i][j].contacts;
  fprintf(STATUS, "REPORTED\t%d\t%d\t%d\n", nclashes, ncontacts, temp); 
  temp=0;
  for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++)
      temp+=debug_contacts[i][j];
  fprintf(STATUS, "DEBUG\t\t%d\t%d\t%d\n", debug_nclashes, debug_ncontacts, temp); 
    
  return;
}
  
void DebugContacts() {
  int i, j;
  debug_nclashes=0;
  debug_ncontacts=0;
  for (i=0; i<natoms; i++)
    for (j=i+1; j<natoms; j++) {
      debug_clashes[i][j]=0;
      debug_contacts[i][j]=0;
      if (data[i][j].check_contacts || data[i][j].check_clashes)
	CheckForDebugContacts(i, j);
    }

  return;
}


void CheckForDebugContacts(int a, int b) {
  
  distance = (native[a].xyz_int.x-native[b].xyz_int.x)*(native[a].xyz_int.x-native[b].xyz_int.x) + (native[a].xyz_int.y-native[b].xyz_int.y)*(native[a].xyz_int.y-native[b].xyz_int.y) + (native[a].xyz_int.z-native[b].xyz_int.z)*(native[a].xyz_int.z-native[b].xyz_int.z);

  if (data[a][b].check_clashes && distance < hard_core[native[a].smogtype][native[b].smogtype]){
    debug_clashes[a][b]=debug_clashes[b][a]=1;
    debug_nclashes++;
  }
  if (data[a][b].check_contacts && (distance <= contact_distance[native[a].smogtype][native[b].smogtype].b) && (distance >= contact_distance[native[a].smogtype][native[b].smogtype].a)) {
    debug_contacts[a][b]=debug_contacts[b][a]=1;
    debug_ncontacts++;
  }

  if (debug_contacts[a][b]!= data[a][b].contacts)
    fprintf(STATUS, "CONTACTS %d, %d --> %f\t%f\t%f\n", a, b, sqrt(distance)/100, sqrt(hard_core[native[a].smogtype][native[b].smogtype])/100, sqrt(contact_distance[native[a].smogtype][native[b].smogtype].b)/100);
  if (debug_clashes[a][b] != data[a][b].clashes)
    fprintf(STATUS, "CLASHES %d, %d --> %ld\t%ld, %ld\n", a, b, distance, hard_core[native[a].smogtype][native[b].smogtype], contact_distance[native[a].smogtype][native[b].smogtype].b);
  
  return;
}


#endif
