BEGIN {
  Bohr2Angs = 0.529177;
}

/0ENDPOINTS/ {
  getline; getline;
  o[0] = $1 * Bohr2Angs;
  o[1] = $2 * Bohr2Angs;
  o[2] = $3 * Bohr2Angs;
  
  getline;
  for(j=0; j<2; j++) {
    for(i=0; i<3; i++)
      v[j,i] = $(i+1) * Bohr2Angs - o[i];
    getline;
  }
}

END {
#
# WARNING: WIEN is writing grid as (den(i,j),j=1,ny),i=1,nx, BUT
#          XCrySDen is reading it as (den(i,j),i=1,nx),j=1,ny
#  FIXING: interchange vec-X with vec-Y and interchange nx with ny
# 
  printf "%f %f %f\n", o[0], o[1], o[2];
  printf "%f %f %f\n", v[1,0], v[1,1], v[1,2];       
  printf "%f %f %f\n", v[0,0], v[0,1], v[0,2];       
}
