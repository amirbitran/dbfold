struct vector {
  Float x, y, z;
};

struct int_vector {
  long int x, y, z;
};

Float D2(struct vector, struct vector);    /* calculates the square of the distance */
void Zero(struct vector *);
void MakeVector(struct vector, struct vector, struct vector *);
Float Dot(struct vector, struct vector);
Float Norm(struct vector);
void Copy(struct vector, struct vector *);
void Inverse(struct vector *);
void Scale(Float, struct vector *);
void Add(struct vector, struct vector *);
void Normalize(struct vector *);
void CrossProduct(struct vector, struct vector, struct vector *);
void bisect(struct vector, struct vector, struct vector *);
void RotateX(Float, struct vector *);
void RotateY(Float, struct vector *);
void RotateZ(Float, struct vector *);
void MakeRotationMatrix(Float, struct vector);
void Rotate(Float[3][3], struct vector *);
Float Angle(struct vector, struct vector);

Float D2(struct vector xyz1, struct vector xyz2) {

  Float a = xyz1.x-xyz2.x;
  Float b = xyz1.y-xyz2.y;
  Float c = xyz1.z-xyz2.z;

  return a*a+b*b+c*c;

}

void Zero(struct vector *A) {

  (*A).x =0;
  (*A).y =0;
  (*A).z =0;

  return;
}

void MakeVector(struct vector A, struct vector B, struct vector *C) {
  
  (*C).x = B.x - A.x;
  (*C).y = B.y - A.y;
  (*C).z = B.z - A.z;

  return;

}

Float Angle(struct vector A, struct vector B) {
  Scale(1/Norm(A), &A);
  Scale(1/Norm(B), &B);
  return (acos(Tiny(Dot(A, B))));
}

Float Dot(struct vector A, struct vector B) {

  return A.x*B.x+A.y*B.y+A.z*B.z;

}

Float Norm(struct vector A) {

  return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);

}

void Copy(struct vector A, struct vector *B) {

  (*B).x = A.x;
  (*B).y = A.y;
  (*B).z = A.z;

  return;

}

void Scale(Float factor, struct vector *A) {

  (*A).x *= factor;
  (*A).y *= factor;
  (*A).z *= factor;

  return;

}

void Inverse(struct vector *A) {

  (*A).x *= -1;
  (*A).y *= -1;
  (*A).z *= -1;

  return;

}

void Add(struct vector A, struct vector *B) {

  (*B).x += A.x;
  (*B).y += A.y;
  (*B).z += A.z;

  return;

}

void Normalize(struct vector *A) {

  Float norm;

  norm = Norm(*A);
  (*A).x /=norm;
  (*A).y /=norm;
  (*A).z /=norm;
  
  return;
}

void CrossProduct(struct vector A, struct vector B, struct vector *C) {

  (*C).x = A.y*B.z - A.z*B.y;
  (*C).y = A.z*B.x - A.x*B.z;
  (*C).z = A.x*B.y - A.y*B.x;

  return;

}
void bisect(struct vector A, struct vector B, struct vector *C) {

  (*C).x = (A.x + B.x)/2.;
  (*C).y = (A.y + B.y)/2.;
  (*C).z = (A.z + B.z)/2.;

  return;

}
  
void RotateZ(Float angle, struct vector *A) {
  
  Float x_new, y_new;

  x_new = (*A).x*cos(angle)+(*A).y*sin(angle);
  y_new = -(*A).x*sin(angle)+(*A).y*cos(angle);

  (*A).x = x_new;
  (*A).y = y_new;

  return;
}

void RotateX(Float angle, struct vector *A) {
  
  Float y_new, z_new;
  
  y_new = (*A).y*cos(angle)-(*A).z*sin(angle);
  z_new = (*A).y*sin(angle)+(*A).z*cos(angle);

  (*A).y = y_new;
  (*A).z = z_new;

  return;
}

void RotateY(Float angle, struct vector *A) {
  
  Float x_new, z_new;
  
  x_new = (*A).x*cos(angle)-(*A).z*sin(angle);
  z_new = (*A).x*sin(angle)+(*A).z*cos(angle);

  (*A).x = x_new;
  (*A).z = z_new;

  return;
}

void MakeRotationMatrix(Float angle, struct vector C) {

  Float u, w, v, f;
  Float t;
  
  t = tan(angle/2.0);
  
  u = t*C.x; 
  v = t*C.y; 
  w = t*C.z; 
  f = 1 / (1 + t*t);

  rot_mat[0][0] = f * (1 + u*u - v*v - w*w);
  rot_mat[0][1] = f * 2 * (u*v - w);
  rot_mat[0][2] = f * 2 * (u*w + v);
  rot_mat[1][0] = f * 2 * (u*v + w);
  rot_mat[1][1] = f * (1 - u*u + v*v - w*w);
  rot_mat[1][2] = f * 2 * (v*w - u);
  rot_mat[2][0] = f * 2 * (u*w - v);
  rot_mat[2][1] = f * 2 * (v*w + u);
  rot_mat[2][2] = f * (1 - u*u - v*v + w*w);

  return;

}
  

void Rotate(Float rot_mat[3][3], struct vector *A) {

  Float x_new, y_new, z_new;
  
  x_new = rot_mat[0][0]*(*A).x + rot_mat[0][1]*(*A).y + rot_mat[0][2]*(*A).z;
  y_new = rot_mat[1][0]*(*A).x + rot_mat[1][1]*(*A).y + rot_mat[1][2]*(*A).z;
  z_new = rot_mat[2][0]*(*A).x + rot_mat[2][1]*(*A).y + rot_mat[2][2]*(*A).z;

  (*A).x = x_new;
  (*A).y = y_new;
  (*A).z = z_new;

  return;

}

