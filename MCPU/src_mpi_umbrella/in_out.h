void read_pdb_backbone(char pdb_name[200], char res_name[5][4], double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res, int end_res);
void write_pdb_backbone(char pdb_name[200], char res_name[5][4], double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res, int end_res);

FILE *fpdb;

//!----------------------------------------------------------------------
//! Copyright (C) 2003 
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.  
//!----------------------------------------------------------------------
//!-----------------------------------------------------------------------
//MODULE in_out
//!----------------------------------------------------------------------------
//  integer, parameter :: dp_io = kind(1.0d0)
//!----------------------------------------------------------------------------
//CONTAINS
//!-----------------------------------------------------------------------
//subroutine read_pdb_backbone(pdb_name, res_name, r_n, r_a, r_c, stt_res, end_res)
//!-----------------------------------------------------------------------
//! read backbone atom coord from a pdb file
//!-----------------------------------------------------------------------
void read_pdb_backbone(char pdb_name[], char res_name[5][4], double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res, int end_res)
 {
//  implicit none
//  character(len=100), intent(in) :: pdb_name
//  integer, intent(in) ::  stt_res, end_res
//  character(len=4), intent(out) :: res_name(:)
//  real(dp_io), intent(out) :: r_n(:,:), r_a(:,:), r_c(:,:)
//  character(len=10) :: string, atmname, resname
  char atmname[10], resname[10];
//  integer :: res_no, unit=10, ioerror, ir
  int res_no, ir;
//  real(dp_io) :: x, y, z
  double x, y, z;
//!-----------------------------------------------------------------------

//  open(unit, file = trim(pdb_name))
  if((fpdb=fopen(pdb_name, "r"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", pdb_name);
    exit(1);
   }
//  do    ! exit at ATOM
//     read(unit, 50) string
//     if (string == 'ATOM') exit
//  end do
//  backspace(unit)

//  do    ! exit if not ATOM
//     read(unit, 50) string
//     if (string /= 'ATOM') exit
//     backspace(unit)
//     read(unit, 56, iostat=ioerror) atmname, resname, res_no, x, y, z
//     if (ioerror < 0) exit 
//     if (res_no > end_res) exit
//     if (res_no >= stt_res) then
//        ir = res_no - stt_res + 1
//        res_name(ir) = resname
//!        print*, ir, atmname
//        if (atmname == 'N') r_n(:, ir) = (/ x, y, z /)
//        if (atmname == 'CA') r_a(:, ir) = (/ x, y, z /)
//        if (atmname == 'C')  r_c(:, ir) = (/ x, y, z /)
//     end if
//  end do
 // close(unit)
  while(1)
   {
    fscanf(fpdb, "%*s%*s %s%s %*s  %d%lf%lf%lf", atmname, resname, &res_no, &x, &y, &z);
    if (res_no>end_res)
     {
      fclose(fpdb);
      return;
     }
    if (res_no>=stt_res)
     {
      ir = res_no - stt_res;
      strcpy(res_name[ir], resname);
      if(strcmp(atmname, "N")==0)
       {
	r_n[ir][0]=x;
	r_n[ir][1]=y;
	r_n[ir][2]=z;
       }
      else if(strcmp(atmname, "CA")==0)
       {
	r_a[ir][0]=x;
	r_a[ir][1]=y;
	r_a[ir][2]=z;
       }
      else if(strcmp(atmname, "C")==0)
       {
	r_c[ir][0]=x;
	r_c[ir][1]=y;
	r_c[ir][2]=z;
       }
      else if(strcmp(atmname, "O")==0)
       {
	r_o[ir][0]=x;
	r_o[ir][1]=y;
	r_o[ir][2]=z;
       }
      else if(strcmp(atmname, "CB")==0)
       {
	r_s[ir][0][0]=x;
	r_s[ir][0][1]=y;
	r_s[ir][0][2]=z;
       }
     }
   }
  
//50 format (a6)
//56 format (13x, a2, 2x, a3, 2x, i4, 4x, 3f8.3)

//end subroutine read_pdb_backbone
 }
//!-----------------------------------------------------------------------
//subroutine write_pdb_backbone(pdb_name, res_name, r_n, r_a, r_c, stt_res, end_res)
void write_pdb_backbone(char pdb_name[], char res_name[5][4], double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3], double r_s[5][MSAR][3], int stt_res, int end_res)
//!-----------------------------------------------------------------------
//! write backbone atom coord in pdb format
//!-----------------------------------------------------------------------
 {
//  implicit none
//  character(len=100), intent(in) :: pdb_name
//  integer, intent(in) :: stt_res, end_res
//  character(len=4), intent(in) :: res_name(:)
//  real(dp_io), intent(in) :: r_n(:, :), r_a(:, :), r_c(:, :)
//  integer :: i, j, k, unit = 20, res_no, ir
  int res_no, ir, k=0;

//  open(unit, file = trim(pdb_name))
  if((fpdb=fopen(pdb_name, "w"))==NULL)
   {
    fprintf(STATUS, "ERROR: Can't open the file: %s!\n", pdb_name);
    exit(1);
   }
	
//  k = 0
  k = 0;
//  do res_no = stt_res, end_res
//     ir = res_no - stt_res + 1
//     k = k+1
//     write(unit, 72) 'ATOM', k, 'N  ', res_name(ir), res_no, r_n(1:3, ir)
//     k = k+1
//     write(unit, 72) 'ATOM', k, 'CA ', res_name(ir), res_no, r_a(1:3, ir)
//     k = k+1
//     write(unit, 72) 'ATOM', k, 'C  ', res_name(ir), res_no, r_c(1:3, ir)
//  end do
//  write(unit, "(a3)") 'END'
//  close (unit)
  for (res_no=stt_res;res_no<=end_res;res_no++)
   {
    ir = res_no - stt_res;
    k++;
    fprintf(fpdb, "ATOM   %4d  N   %s A%4d    %8.3lf%8.3lf%8.3lf\n", 
	    k, res_name[ir], res_no, r_n[ir][0], r_n[ir][1], r_n[ir][2]);
    k++;
    fprintf(fpdb, "ATOM   %4d  CA  %s A%4d    %8.3lf%8.3lf%8.3lf\n", 
	    k, res_name[ir], res_no, r_a[ir][0], r_a[ir][1], r_a[ir][2]);
    k++;
    fprintf(fpdb, "ATOM   %4d  C   %s A%4d    %8.3lf%8.3lf%8.3lf\n", 
	    k, res_name[ir], res_no, r_c[ir][0], r_c[ir][1], r_c[ir][2]);
    fprintf(fpdb, "ATOM   %4d  O   %s A%4d    %8.3lf%8.3lf%8.3lf\n", 
	    k, res_name[ir], res_no, r_o[ir][0], r_o[ir][1], r_o[ir][2]);
    fprintf(fpdb, "ATOM   %4d  CB  %s A%4d    %8.3lf%8.3lf%8.3lf\n", 
	    k, res_name[ir], res_no, r_s[ir][0][0], r_s[ir][0][1], r_s[ir][0][2]);
   }
  fclose(fpdb);
  return;
//72 format(a4, 3x, i4, 2x, a3, 1x, a3, 2x, i4, 4x, 3f8.3)

//end subroutine write_pdb_backbone
//!-----------------------------------------------------------------------
//END MODULE in_out
//!-----------------------------------------------------------------------
 }
