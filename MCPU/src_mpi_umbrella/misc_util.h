Float Tiny(Float x);
void squeeze(char *, int);
Float GaussianNum();

Float Tiny(Float x) {

  return ((Float) ((int) (x*PRECISION)))/PRECISION;
}

void squeeze(char s[], int c) {

   int i, j;
   for (i = j= 0; s[i] != '\0'; i++)
     if (s[i] != c)
       s[j++] = s[i];
   s[j] = '\0';

   return;
}

Float GaussianNum() {

  static int iset=0;
  static Float gset;
  
  Float fac, rsq, v1, v2;
  
  if (iset==0) {
    
    do {
      
      v1 = 2.0*threefryrand()-1.0; 
      v2 = 2.0*threefryrand()-1.0;
      
      rsq = v1*v1 + v2*v2;

    } while (rsq >= 1.0 || rsq == 0.0);
    
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    
    return v2*fac;

  }
  else {
    iset=0;
    return gset;
  }

}
