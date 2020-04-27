void DoRotation(int, int, int, int, Float, short, short *);

void DoRotation(int a, int b, int c, int d, Float delta_angle, short rotate_natoms, short *rotate_atom) {
  int j;
  Float uu, ww, vv, f;
  Float t, norm;
  struct vector temp_xyz, u, v;

  /* make rotation matrix */
  temp_xyz.x = native[b].xyz.x;
  temp_xyz.y = native[b].xyz.y;
  temp_xyz.z = native[b].xyz.z;

  v.x = native[c].xyz.x-temp_xyz.x;
  v.y = native[c].xyz.y-temp_xyz.y;
  v.z = native[c].xyz.z-temp_xyz.z;

  norm = 1/sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  v.x *= norm;
  v.y *= norm;
  v.z *= norm;  /*We have computed vector between b and c, and now we normalize it!*/

  t = tan(delta_angle*0.5);
  
  uu = t*v.x; 
  vv = t*v.y; 
  ww = t*v.z; 
  f = 2 / (1 + t*t);

  rot_mat_00 = 0.5 * f * (1 + uu*uu - vv*vv - ww*ww);
  rot_mat_01 = f * (uu*vv - ww);
  rot_mat_02 = f * (uu*ww + vv);
  rot_mat_10 = f * (uu*vv + ww);
  rot_mat_11 = 0.5 * f * (1 - uu*uu + vv*vv - ww*ww);
  rot_mat_12 = f * (vv*ww - uu);
  rot_mat_20 = f * (uu*ww - vv);
  rot_mat_21 = f * (vv*ww + uu);
  rot_mat_22 = 0.5 * f * (1 - uu*uu - vv*vv + ww*ww);

  /* rotate each vector */

  for (j=0; j<rotate_natoms; j++) {
    temp_atom = &native[rotate_atom[j]];
    u.x = temp_atom->xyz.x - temp_xyz.x;
    u.y = temp_atom->xyz.y - temp_xyz.y;
    u.z = temp_atom->xyz.z - temp_xyz.z;
    temp_atom->xyz.x = rot_mat_00*u.x + rot_mat_01*u.y + rot_mat_02*u.z + temp_xyz.x;
    temp_atom->xyz.y = rot_mat_10*u.x + rot_mat_11*u.y + rot_mat_12*u.z + temp_xyz.y;
    temp_atom->xyz.z = rot_mat_20*u.x + rot_mat_21*u.y + rot_mat_22*u.z + temp_xyz.z;
  }

  return;

}

 
