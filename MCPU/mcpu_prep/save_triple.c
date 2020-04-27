/* Program: save_triple.c
 * Purpose: To convert triple.energy data from residue specific to general one
 * */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int seq[1000];
char line[250], line_tar[250];
char aa_singlet[] = "ARNDCQEGHILKMFPSTWYV";

int main(int argc, char *argv[]) {

	int i, j, n;
	int energy;
	int ang_PCA, ang_bCA, phim, psim;
	int obs, tot_obs;
	int r1, r2, r3;
	int res1, res2, res3;
	int ang[4];
  
  	char fn[100];
  	
	FILE *f_org, *f_tar;
	
	if(argc != 2) {
		printf("Usage: save_triple <PDB_ID>\n");
		exit(1);
	}
	
	
	sprintf(fn, "%s.fasta", argv[1]);
	f_tar = fopen(fn,"r");
	
	n=0;
	while(fgets(line_tar,200,f_tar)!=NULL) {
		assert(line_tar);
		if(line_tar[0] == '>') continue;

		printf("%s", line_tar);
		i=0;
		while(line_tar[i] != '\n') {
			for(j=0;j<20;j++) {
				if(aa_singlet[j] == line_tar[i]) break;
			}
			printf("%d %c %d\n", n, line_tar[i], j);
			seq[n++] = j;
			i++;
		}
		
	}
	
	fprintf(stdout, "Total %d aa reads.\n", n);
	
	fclose(f_tar);
	
	sprintf(fn, "%s.triple", argv[1]);
	f_tar = fopen(fn, "w");
	
	for(i=0;i<n-2;i++) {
		
		
		r1 = seq[i];
		r2 = seq[i+1];
		r3 = seq[i+2];
	
    	f_org = fopen("triple.energy", "r");
    	assert(f_org);
    	

    	while(fgets(line,200,f_org)!=NULL) {
    		sscanf(line,"%d%d%d %d%d%d%d %d %d %d %*s", 
    			&res1, &res2, &res3, &ang_PCA, &ang_bCA, &phim, &psim, &energy, &obs, &tot_obs);
      		//printf("%s", line);
      		//printf("res1 = %d, res2 = %d, res3 = %d\n",res1, res2, res3);
      		if((r1==res1)&&(r2==res2)&&(r3==res3))
        		fprintf(f_tar, "%4d %2d %2d %2d %3d %3d %3d %3d %5d %7d %7d\n",
               		i, res1, res2, res3, ang_PCA, ang_bCA, phim, psim, energy, obs, tot_obs);
     	}
		
		fclose(f_org);
	}
	
	fclose(f_tar);
	
	sprintf(fn, "%s.sctorsion", argv[1]);
	f_tar = fopen(fn, "w");
	
	for(i=0;i<n-3;i++) {
			
		r1 = seq[i];
		r2 = seq[i+1];
		r3 = seq[i+2];

		f_org = fopen("sct.energy","r");
    	assert(f_org);

		while(fgets(line,200,f_org)!=NULL)	{
		
			sscanf(line,"%d%d%d %d%d%d%d %d %d %d %*s",
	     			&res1, &res2, &res3 , &ang[0], &ang[1], &ang[2], &ang[3], &energy, &obs, &tot_obs);
			if((r1==res1)&&(r2==res2)&&(r3==res3))
        		fprintf(f_tar, "%4d %2d %2d %2d %3d %3d %3d %3d %5d %7d %7d\n",
               		i, res1, res2, res3, ang[0], ang[1], ang[2], ang[3], energy, obs, tot_obs);

		}
		fclose(f_org);
	}
     
	fclose(f_tar);
     
     

	return 0;
}

