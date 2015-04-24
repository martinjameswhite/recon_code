#include	<cmath>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"



void	shift_obj(std::vector<struct particle>& P,
                  const std::vector<float>& phi) {
// Computes the displacement field from phi using second-order accurate
// finite difference and shifts the data and randoms.
// The displacements are pulled from the grid onto the positions of the
// particles using CIC.
  const float fac = 0.5*Ng;	// The 1/(2h) factor for the derivatives.
  float sum=0;
#pragma omp parallel for shared(P,phi) reduction(+:sum)
  for (int nn=0; nn<P.size(); ++nn) {
    // This is written out in gory detail both to make it easier to
    // see what's going on and to encourage the compiler to optimize
    // and vectorize the code as much as possible.
    int  ix0 = Ng*P[nn].pos[0]; float dx = Ng*P[nn].pos[0]-ix0;
    int  iy0 = Ng*P[nn].pos[1]; float dy = Ng*P[nn].pos[1]-iy0;
    int  iz0 = Ng*P[nn].pos[2]; float dz = Ng*P[nn].pos[2]-iz0;
    int  ixp = (ix0+1   )%Ng;
    int  ixP = (ix0+2   )%Ng;
    int  ixm = (ix0-1+Ng)%Ng;
    int  iyp = (iy0+1   )%Ng;
    int  iyP = (iy0+2   )%Ng;
    int  iym = (iy0-1+Ng)%Ng;
    int  izp = (iz0+1   )%Ng;
    int  izP = (iz0+2   )%Ng;
    int  izm = (iz0-1+Ng)%Ng;
    float px,py,pz,wt;
    wt = (1-dx)*(1-dy)*(1-dz);
    px = (phi[Ng*Ng*ixp+Ng*iy0+iz0]-phi[Ng*Ng*ixm+Ng*iy0+iz0])*fac*wt;
    py = (phi[Ng*Ng*ix0+Ng*iyp+iz0]-phi[Ng*Ng*ix0+Ng*iym+iz0])*fac*wt;
    pz = (phi[Ng*Ng*ix0+Ng*iy0+izp]-phi[Ng*Ng*ix0+Ng*iy0+izm])*fac*wt;
    wt = ( dx )*(1-dy)*(1-dz);
    px+= (phi[Ng*Ng*ixP+Ng*iy0+iz0]-phi[Ng*Ng*ix0+Ng*iy0+iz0])*fac*wt;
    py+= (phi[Ng*Ng*ixp+Ng*iyp+iz0]-phi[Ng*Ng*ixp+Ng*iym+iz0])*fac*wt;
    pz+= (phi[Ng*Ng*ixp+Ng*iy0+izp]-phi[Ng*Ng*ixp+Ng*iy0+izm])*fac*wt;
    wt = (1-dx)*( dy )*(1-dz);
    px+= (phi[Ng*Ng*ixp+Ng*iyp+iz0]-phi[Ng*Ng*ixm+Ng*iyp+iz0])*fac*wt;
    py+= (phi[Ng*Ng*ix0+Ng*iyP+iz0]-phi[Ng*Ng*ix0+Ng*iy0+iz0])*fac*wt;
    pz+= (phi[Ng*Ng*ix0+Ng*iyp+izp]-phi[Ng*Ng*ix0+Ng*iyp+izm])*fac*wt;
    wt = (1-dx)*(1-dy)*( dz );
    px+= (phi[Ng*Ng*ixp+Ng*iy0+izp]-phi[Ng*Ng*ixm+Ng*iy0+izp])*fac*wt;
    py+= (phi[Ng*Ng*ix0+Ng*iyp+izp]-phi[Ng*Ng*ix0+Ng*iym+izp])*fac*wt;
    pz+= (phi[Ng*Ng*ix0+Ng*iy0+izP]-phi[Ng*Ng*ix0+Ng*iy0+iz0])*fac*wt;
    wt = ( dx )*( dy )*(1-dz);
    px+= (phi[Ng*Ng*ixP+Ng*iyp+iz0]-phi[Ng*Ng*ix0+Ng*iyp+iz0])*fac*wt;
    py+= (phi[Ng*Ng*ixp+Ng*iyP+iz0]-phi[Ng*Ng*ixp+Ng*iy0+iz0])*fac*wt;
    pz+= (phi[Ng*Ng*ixp+Ng*iyp+izp]-phi[Ng*Ng*ixp+Ng*iyp+izm])*fac*wt;
    wt = ( dx )*(1-dy)*( dz );
    px+= (phi[Ng*Ng*ixP+Ng*iy0+izp]-phi[Ng*Ng*ix0+Ng*iy0+izp])*fac*wt;
    py+= (phi[Ng*Ng*ixp+Ng*iyp+izp]-phi[Ng*Ng*ixp+Ng*iym+izp])*fac*wt;
    pz+= (phi[Ng*Ng*ixp+Ng*iy0+izP]-phi[Ng*Ng*ixp+Ng*iy0+iz0])*fac*wt;
    wt = (1-dx)*( dy )*( dz );
    px+= (phi[Ng*Ng*ixp+Ng*iyp+izp]-phi[Ng*Ng*ixm+Ng*iyp+izp])*fac*wt;
    py+= (phi[Ng*Ng*ix0+Ng*iyP+izp]-phi[Ng*Ng*ix0+Ng*iy0+izp])*fac*wt;
    pz+= (phi[Ng*Ng*ix0+Ng*iyp+izP]-phi[Ng*Ng*ix0+Ng*iyp+iz0])*fac*wt;
    wt = ( dx )*( dy )*( dz );
    px+= (phi[Ng*Ng*ixP+Ng*iyp+izp]-phi[Ng*Ng*ix0+Ng*iyp+izp])*fac*wt;
    py+= (phi[Ng*Ng*ixp+Ng*iyP+izp]-phi[Ng*Ng*ixp+Ng*iy0+izp])*fac*wt;
    pz+= (phi[Ng*Ng*ixp+Ng*iyp+izP]-phi[Ng*Ng*ixp+Ng*iyp+iz0])*fac*wt;
    float rx = box.ctr[0]+box.L*(P[nn].pos[0]-0.5);
    float ry = box.ctr[1]+box.L*(P[nn].pos[1]-0.5);
    float rz = box.ctr[2]+box.L*(P[nn].pos[2]-0.5);
    float r2 = rx*rx+ry*ry+rz*rz;
    float pr = beta*(px*rx+py*ry+pz*rz);
    P[nn].pos[0] -= px + pr*rx/r2;
    P[nn].pos[1] -= py + pr*ry/r2;
    P[nn].pos[2] -= pz + pz*rz/r2;
    sum += px*px+py*py+pz*pz;
  }
  sum = sqrt(sum/3/P.size());
  std::cout<<"# RMS 1D shift is "<<sum*box.L<<" Mpc/h."<<std::endl;
}
