#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau); 

struct fragments {
  int x1; 
  int x2; 
}; 

struct alignment {
  struct fragments seqptr[MAXFRAG];
  struct fragments structptr[MAXFRAG];
  int NFRAG;
};

 
struct backbone {
  struct vector N;
  struct vector CA;
  struct vector C;
  struct vector CB;
  struct vector O;
};

float getrms(struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn);
int kearsley(float a[5][5], struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn);
void jacobi(float a[5][5], int n, float d[5], float v[5][5], int *nrot);

int NFRAG; 

float getrms(struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn){
  // given the alignment, program RMS aligns the two protein backbone chains
  // for the C-alphas
  // (1) construct the kearsley matrix
  // (2) run the jacobi rotation with the kearsley matrix
  // (3) output the lowest eigenvalue

  int i;
  float a[5][5], d[5], v[5][5];
  int nrot = 0;
  int algnlength = 0;
  float lowest = 10000000000.;
  int nnn = 4;

  kearsley(a, struct1, struct2, algn);  // return the kearsley 'a' matrix
  jacobi(a, nnn, d, v, &nrot);
  
  for (i=1; i<=NFRAG; i++){
    algnlength += algn.seqptr[i].x2 - algn.seqptr[i].x1 + 1; 
  }
  
  // find lowest eigenvalue in d vector
  for (i=1; i<=4; i++) {
    if (d[i] < lowest) lowest = d[i];
  }
  lowest = sqrt(lowest/algnlength);
  return lowest;
}

int kearsley(float a[5][5], struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn){
  int i, j, k, algnlength = 0;
  struct vector centroid1, centroid2;
//  float xmxm = 0, ymym = 0, zmzm = 0, ypzm = 0, ymzp = 0, xmzp = 0, xpzm = 0;
//  float ypyp = 0, zpzp = 0, xmym = 0, xpyp = 0, xpxp = 0;
//  float xpym = 0, xmyp = 0, xmzm = 0, xpzp = 0, ymzm = 0, ypzp = 0;
  float xm, xp, ym, yp, zm, zp;
  
  // first translate the two molecules in the alignment region
  // struct1 is the real structure of the SEQUENCE
  // struct2 is the structure of teh TEMPLATE
 
  centroid1.x = 0; centroid1.y = 0; centroid1.z = 0;

  for (i=0; i<=4; i++)
    for (j=0; j<=4; j++)
      a[i][j] = 0;

  for (i=1; i<=NFRAG; i++)
    for (j=algn.seqptr[i].x1; j<=algn.seqptr[i].x2; j++)
      algnlength++;
 
  for (i=1; i<=NFRAG; i++){
    for (j=algn.seqptr[i].x1; j<=algn.seqptr[i].x2; j++){
      if(struct1[j].CA.x > -9990){
        centroid1.x += struct1[j].CA.x/algnlength;
        centroid1.y += struct1[j].CA.y/algnlength;
        centroid1.z += struct1[j].CA.z/algnlength;
      }
    }
  }
 
  for (i=1; i<=NFRAG; i++){
    for (j=algn.seqptr[i].x1; j<=algn.seqptr[i].x2; j++){
      struct1[j].CA.x -= centroid1.x;  
      struct1[j].CA.y -= centroid1.y;  
      struct1[j].CA.z -= centroid1.z;  
    }
  }

  centroid2.x = 0; centroid2.y = 0; centroid2.z = 0;

  for (i=1; i<=NFRAG; i++){
    for (j=algn.structptr[i].x1; j<=algn.structptr[i].x2; j++){
      if(struct2[j].CA.x > -9990){
        centroid2.x += struct2[j].CA.x/algnlength;
        centroid2.y += struct2[j].CA.y/algnlength;
        centroid2.z += struct2[j].CA.z/algnlength;
      }
    }
  }

  for (i=1; i<=NFRAG; i++){
    for (j=algn.structptr[i].x1; j<=algn.structptr[i].x2; j++){
      struct2[j].CA.x -= centroid2.x;
      struct2[j].CA.y -= centroid2.y;
      struct2[j].CA.z -= centroid2.z;
    }
  }

  // ok, translated, now construct the matrix into 'a'

  for (i=1; i<=NFRAG; i++){
    for (j=algn.seqptr[i].x1; j<=algn.seqptr[i].x2; j++){
      k = j - algn.seqptr[i].x1 + algn.structptr[i].x1;
      if(struct1[j].CA.x > -9990 && struct2[k].CA.x > -9990){
        xm = struct1[j].CA.x - struct2[k].CA.x; xp = struct1[j].CA.x + struct2[k].CA.x;
        ym = struct1[j].CA.y - struct2[k].CA.y; yp = struct1[j].CA.y + struct2[k].CA.y;
        zm = struct1[j].CA.z - struct2[k].CA.z; zp = struct1[j].CA.z + struct2[k].CA.z;
        a[1][1] += xm*xm + ym*ym + zm*zm;
        a[1][2] += yp*zm - ym*zp;
        a[1][3] += xm*zp - xp*zm;
        a[1][4] += xp*ym - xm*yp;
        a[2][2] += yp*yp + zp*zp + xm*xm;
        a[2][3] += xm*ym - xp*yp;
        a[2][4] += xm*zm - xp*zp;
        a[3][3] += xp*xp + zp*zp + ym*ym;
        a[3][4] += ym*zm - yp*zp;
        a[4][4] += xp*xp + yp*yp + zm*zm;
      }
    }
  }

  // symmetrize

  for (i=2; i<=4; i++)
    for (j=1; j<i; j++){
	a[i][j] = a[j][i];
    }
//  fprintf(STATUS, "a matrix:\n"); 
//  for (i=1; i<=4; i++)
//    fprintf(STATUS, "%12.2f\t%12.2f\t%12.2f\t%12.2f\n", a[i][1], a[i][2], a[i][3], a[i][4]);

  for (i=1; i<=NFRAG; i++){
    for (j=algn.seqptr[i].x1; j<=algn.seqptr[i].x2; j++){
      struct1[j].CA.x += centroid1.x;
      struct1[j].CA.y += centroid1.y;
      struct1[j].CA.z += centroid1.z;
    }
  }

  for (i=1; i<=NFRAG; i++){
    for (j=algn.structptr[i].x1; j<=algn.structptr[i].x2; j++){
      struct2[j].CA.x += centroid2.x;
      struct2[j].CA.y += centroid2.y;
      struct2[j].CA.z += centroid2.z;
    }
  }
  
  return 0;
}

void jacobi(float a[5][5], int n, float d[5], float v[5][5], int *nrot){
  // pass in address of nrot to modify

  int j, iq, ip, i;
  float tresh, theta, tau, t, sm, s, h, g, c, b[5], z[5];

  // jacobi should clean up a bunch of functions
  for (ip=1; ip<=4; ip++){
    for (iq=1; iq<=4; iq++) v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }
  for (ip=1; ip<=4; ip++){
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  *nrot = 0;

  for (i=1; i<=50; i++) {  // 50 jacobi iterations
    sm = 0.0;
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0)  {
      return;   // reached threshhold, no need to free b, z, whatever since we hardcode.
    }
    if (i < 4) 
      tresh = 0.2*sm/(n*n);
    else
      tresh = 0.0;
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++) {
        g = 100.0*fabs(a[ip][iq]);
        if (i>4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
           && (float)(fabs(d[iq])+g) == (float) fabs(d[iq]))
           a[ip][iq] = 0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq]-d[ip];
          if ((float)(fabs(h)+g) == (float) fabs(h))
            t = (a[ip][iq])/h;
          else {
            theta = 0.5*h/(a[ip][iq]);
            t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0)
              t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0.0;
          for (j=1; j<=ip-1; j++) {
            ROTATE(a, j, ip, j, iq)
          }
          for (j=ip+1; j<=iq-1; j++) {
            ROTATE(a, ip, j, j, iq)
          }
          for (j=iq+1; j<=n; j++) {
            ROTATE(a, ip, j, iq, j)
          }
          for (j=1; j<=n; j++) {
            ROTATE(v, j, ip, j, iq)
          }
          ++(*nrot);
        }
      }
    }
    for (ip=1; ip<=n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  fprintf(stderr, "Too many iterations.\n");
  fprintf(stderr, "ABORT.\n");
  exit(1);
}
