#define XYZ_ApTBmA(dst,a,t,b) \
do { \
  dst.x = a.x + t*(b.x-a.x); \
  dst.y = a.y + t*(b.y-a.y); \
  dst.z = a.z + t*(b.z-a.z); \
 }  while(0)

void trilinear()
{
  
  /* 
     tri-linear interpolation 
  */
  
  /* z-dir */
  
  XYZ_ApTBmA(u0, grad[0], iver[2], grad[4]); /* u0 = grad[0] + iver[2]*(grad[4]-grad[0]); */
  XYZ_ApTBmA(u1, grad[1], iver[2], grad[5]); /* u1 = grad[1] + iver[2]*(grad[5]-grad[1]); */
  XYZ_ApTBmA(u2, grad[2], iver[2], grad[6]); /* u2 = grad[2] + iver[2]*(grad[6]-grad[2]); */
  XYZ_ApTBmA(u3, grad[3], iver[2], grad[7]); /* u3 = grad[3] + iver[2]*(grad[7]-grad[3]); */
  
  /* y-dir */
  XYZ_ApTBmA(v0, u0, iver[1], u1); /* v0 = u0 + iver[1]*(u1-u0); */
  XYZ_ApTBmA(v1, u3, iver[1], u2); /* v1 = u3 + iver[1]*(u2-u3); */
  
  /* x-dir */
  XYZ_ApTBmA(h, v0, iver[0], v1); /*  h = v0 + iver[0]*(v1-v0); */
}
  
