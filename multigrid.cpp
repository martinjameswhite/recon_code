#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<sstream>
#include	<algorithm>
#include	<vector>
#include	<exception>

#include	"global.h"
#include	"multigrid.h"

// The multigrid code for solving our modified Poisson-like equation.
// See the notes for details.
// The only place that the equations appear explicitly is in
//	gauss_seidel, jacobi and residual
//
//
// The finite difference equation we are solving is written schematically
// as A.v=f where A is the matrix operator formed from the discretized
// partial derivatives and f is the source term (density in our case).
// The matrix, A, is never explicitly constructed.
// The solution, phi, is in v.
//
//
// Author:	Martin White	(UCB/LBNL)
// Written:	20-Apr-2015
// Modified:	20-Apr-2015
//





void MultiGrid::gauss_seidel(std::vector<float>& v,
                             const std::vector<float>& f,
                             const int N) {
// NOTE: For the cross terms the red-black scheme doesn't quite work ...
// for this reason this code isn't being used currently.
// Does an update using red-black Gauss-Seidel (RBGS).  This, and in
// residual below, is where the explicit equation we are solving appears.
// See notes for more details.
  const float h  = 1.0/N;
  const float h2 = h*h;
  const float len= box.L*h;
  const int   Nit= 5;
  for (int iter=0; iter<Nit; ++iter) {
    for (int rb=0; rb<2; ++rb) {
#pragma omp parallel for shared(v,f,rb)
      for (int ix=0; ix<N; ++ix) {
        float rx=box.ctr[0]+len*(ix-N/2);
        int ixp =(ix+1  )%N;
        int ixm =(ix-1+N)%N;
        for (int iy=0; iy<N; ++iy) {
          float ry=box.ctr[1]+len*(iy-N/2);
          int iyp =(iy+1  )%N;
          int iym =(iy-1+N)%N;
          int st  =((ix+iy)%2==rb)?0:1;
          for (int iz=st; iz<N; iz+=2) {
            float rz=box.ctr[2]+len*(iz-N/2);
            float g =beta/(rx*rx+ry*ry+rz*rz+1e-30);
            int izp=(iz+1  )%N;
            int izm=(iz-1+N)%N;
            int ii = N*N*ix+N*iy+iz;
            v[ii] = h2*f[ii]+
                    (1+g*rx*rx)*(v[N*N*ixp+N*iy +iz ]+v[N*N*ixm+N*iy +iz ])+
                    (1+g*ry*ry)*(v[N*N*ix +N*iyp+iz ]+v[N*N*ix +N*iym+iz ])+
                    (1+g*rz*rz)*(v[N*N*ix +N*iy +izp]+v[N*N*ix +N*iy +izm])+
                    (g*rx*ry/2)*(v[N*N*ixp+N*iyp+iz ]+v[N*N*ixm+N*iym+iz ]
                                -v[N*N*ixm+N*iyp+iz ]-v[N*N*ixp+N*iym+iz ])+
                    (g*rx*rz/2)*(v[N*N*ixp+N*iy +izp]+v[N*N*ixm+N*iy +izm]
                                -v[N*N*ixm+N*iy +izp]-v[N*N*ixp+N*iy +izm])+
                    (g*ry*rz/2)*(v[N*N*ix +N*iyp+izp]+v[N*N*ix +N*iym+izm]
                                -v[N*N*ix +N*iym+izp]-v[N*N*ix +N*iyp+izm]);
            v[ii]/= 6+2*beta;
          }
        }
      }
    }
  }
}





void MultiGrid::jacobi(std::vector<float>& v,
                       const std::vector<float>& f, const int N) {
// Does an update using damped Jacobi.  This, and in residual below,
// is where the explicit equation we are solving appears.
// See notes for more details.
  std::vector<float> jac;
  try{jac.resize(N*N*N);}catch(std::exception& e) {myexception(e);}
  const float w  = 0.4;	// The "damping" factor.
  const float h  = 1.0/N;
  const float h2 = h*h;
  const float len= box.L*h;
  const int   Nit= 5;
  for (int iter=0; iter<Nit; ++iter) {
#pragma omp parallel for shared(v,f,jac)
    for (int ix=0; ix<N; ++ix) {
      float rx= box.ctr[0]+len*(ix-N/2);
      int ixp = (ix+1  )%N;
      int ixm = (ix-1+N)%N;
      for (int iy=0; iy<N; ++iy) {
        float ry= box.ctr[1]+len*(iy-N/2);
        int iyp = (iy+1  )%N;
        int iym = (iy-1+N)%N;
        for (int iz=0; iz<N; ++iz) {
          float rz= box.ctr[2]+len*(iz-N/2);
          float g = beta/(rx*rx+ry*ry+rz*rz+1e-30);
          int izp = (iz+1  )%N;
          int izm = (iz-1+N)%N;
          int ii  = N*N*ix+N*iy+iz;
          jac[ii] = h2*f[ii]+
                    (1+g*rx*rx)*(v[N*N*ixp+N*iy +iz ]+v[N*N*ixm+N*iy +iz ])+
                    (1+g*ry*ry)*(v[N*N*ix +N*iyp+iz ]+v[N*N*ix +N*iym+iz ])+
                    (1+g*rz*rz)*(v[N*N*ix +N*iy +izp]+v[N*N*ix +N*iy +izm])+
                    (g*rx*ry/2)*(v[N*N*ixp+N*iyp+iz ]+v[N*N*ixm+N*iym+iz ]
                                -v[N*N*ixm+N*iyp+iz ]-v[N*N*ixp+N*iym+iz ])+
                    (g*rx*rz/2)*(v[N*N*ixp+N*iy +izp]+v[N*N*ixm+N*iy +izm]
                                -v[N*N*ixm+N*iy +izp]-v[N*N*ixp+N*iy +izm])+
                    (g*ry*rz/2)*(v[N*N*ix +N*iyp+izp]+v[N*N*ix +N*iym+izm]
                                -v[N*N*ix +N*iym+izp]-v[N*N*ix +N*iyp+izm])+
                    (g*rx*len) *(v[N*N*ixp+N*iy +iz ]-v[N*N*ixm+N*iy +iz ])+
                    (g*ry*len) *(v[N*N*ix +N*iyp+iz ]-v[N*N*ix +N*iym+iz ])+
                    (g*rz*len) *(v[N*N*ix +N*iy +izp]-v[N*N*ix +N*iy +izm]);
          jac[ii]/= 6+2*beta;
        }
      }
    }
#pragma omp parallel for shared(v,jac)
    for (int nn=0; nn<v.size(); ++nn)
      v[nn] = (1-w)*v[nn] + w*jac[nn];
  }
}





std::vector<float> MultiGrid::residual(const std::vector<float>& v,
                                       const std::vector<float>& f,
                                       const int N) {
// Returns the residual, r=f-Av, keeping track of factors of h=1/N.
  const float h  = 1.0/N;
  const float len= box.L*h;
  std::vector<float> r;
  try {r.resize(N*N*N);}catch(std::exception& e) {myexception(e);}
  // First compute the operator A on v, keeping track of periodic
  // boundary conditions and ignoring the 1/h terms.
  // Note the relative signs here and in jacobi (or gauss_seidel).
#pragma omp parallel for shared(r,v,f)
  for (int ix=0; ix<N; ++ix) {
    float rx=box.ctr[0]+len*(ix-N/2);
    int ixp =(ix+1  )%N;
    int ixm =(ix-1+N)%N;
    for (int iy=0; iy<N; ++iy) {
      float ry=box.ctr[1]+len*(iy-N/2);
      int iyp =(iy+1  )%N;
      int iym =(iy-1+N)%N;
      for (int iz=0; iz<N; ++iz) {
        float rz=box.ctr[2]+len*(iz-N/2);
        float g =beta/(rx*rx+ry*ry+rz*rz+1e-30);
        int izp =(iz+1  )%N;
        int izm =(iz-1+N)%N;
        int ii  = N*N*ix+N*iy+iz;
        r[ii] = (6+ 2*beta)* v[ii] -
                (1+g*rx*rx)*(v[N*N*ixp+N*iy +iz ]+v[N*N*ixm+N*iy +iz ])-
                (1+g*ry*ry)*(v[N*N*ix +N*iyp+iz ]+v[N*N*ix +N*iym+iz ])-
                (1+g*rz*rz)*(v[N*N*ix +N*iy +izp]+v[N*N*ix +N*iy +izm])-
                (g*rx*ry/2)*(v[N*N*ixp+N*iyp+iz ]+v[N*N*ixm+N*iym+iz ]
                            -v[N*N*ixm+N*iyp+iz ]-v[N*N*ixp+N*iym+iz ])-
                (g*rx*rz/2)*(v[N*N*ixp+N*iy +izp]+v[N*N*ixm+N*iy +izm]
                            -v[N*N*ixm+N*iy +izp]-v[N*N*ixp+N*iy +izm])-
                (g*ry*rz/2)*(v[N*N*ix +N*iyp+izp]+v[N*N*ix +N*iym+izm]
                            -v[N*N*ix +N*iym+izp]-v[N*N*ix +N*iyp+izm])-
                (g*rx*len) *(v[N*N*ixp+N*iy +iz ]-v[N*N*ixm+N*iy +iz ])-
                (g*ry*len) *(v[N*N*ix +N*iyp+iz ]-v[N*N*ix +N*iym+iz ])-
                (g*rz*len) *(v[N*N*ix +N*iy +izp]-v[N*N*ix +N*iy +izm]);
      }
    }
  }
  // Now subtract it from f -- remember the 1/h^2 from the derivatives.
#pragma omp parallel for shared(r,f)
  for (int nn=0; nn<r.size(); ++nn)
    r[nn] = f[nn] - r[nn]*N*N;
  return(r);
}






std::vector<float> MultiGrid::prolong(const std::vector<float>& v2h,
                                      const int N) {
// Transfer a vector, v2h, from the coarse grid with spacing 2h to a
// fine grid with spacing 1h using linear interpolation and periodic BC.
// The length, N, is of the coarse-grid vector, v2h.
// This is simple, linear interpolation in a cube.
  const int N2=2*N;
  std::vector<float> v1h;
  try {v1h.resize(N2*N2*N2);} catch(std::exception& e) {myexception(e);}
#pragma omp parallel for shared(v2h,v1h)
  for (int ix=0; ix<N; ++ix) {
    int ixp=(ix+1)%N;
    for (int iy=0; iy<N; ++iy) {
      int iyp=(iy+1)%N;
      for (int iz=0; iz<N; ++iz) {
        int izp=(iz+1)%N;
        v1h[N2*N2*(2*ix+0)+N2*(2*iy+0)+(2*iz+0)]= v2h[N*N*ix +N*iy +iz ];
        v1h[N2*N2*(2*ix+1)+N2*(2*iy+0)+(2*iz+0)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ixp+N*iy +iz ])/2;
        v1h[N2*N2*(2*ix+0)+N2*(2*iy+1)+(2*iz+0)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ix +N*iyp+iz ])/2;
        v1h[N2*N2*(2*ix+0)+N2*(2*iy+0)+(2*iz+1)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ix +N*iy +izp])/2;
        v1h[N2*N2*(2*ix+1)+N2*(2*iy+1)+(2*iz+0)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ixp+N*iy +iz ]+
                                                  v2h[N*N*ix +N*iyp+iz ]+
                                                  v2h[N*N*ixp+N*iyp+iz ])/4;
        v1h[N2*N2*(2*ix+0)+N2*(2*iy+1)+(2*iz+1)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ix +N*iyp+iz ]+
                                                  v2h[N*N*ix +N*iy +izp]+
                                                  v2h[N*N*ix +N*iyp+izp])/4;
        v1h[N2*N2*(2*ix+1)+N2*(2*iy+0)+(2*iz+1)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ixp+N*iy +iz ]+
                                                  v2h[N*N*ix +N*iy +izp]+
                                                  v2h[N*N*ixp+N*iy +izp])/4;
        v1h[N2*N2*(2*ix+1)+N2*(2*iy+1)+(2*iz+1)]=(v2h[N*N*ix +N*iy +iz ]+
                                                  v2h[N*N*ixp+N*iy +iz ]+
                                                  v2h[N*N*ix +N*iyp+iz ]+
                                                  v2h[N*N*ix +N*iy +izp]+
                                                  v2h[N*N*ixp+N*iyp+iz ]+
                                                  v2h[N*N*ixp+N*iy +izp]+
                                                  v2h[N*N*ix +N*iyp+izp]+
                                                  v2h[N*N*ixp+N*iyp+izp])/8;
      }
    }
  }
  return(v1h);
}




std::vector<float> MultiGrid::restrict(const std::vector<float>& v1h,
                                       const int N) {
// Transfer a vector, v1h, from the fine grid with spacing 1h to a coarse
// grid with spacing 2h using full weighting and periodic BC.
// The length, N, is of the fine-grid vector (v1h) and is assumed even,
// the code doesn't check.
  const int N2=N/2;
  std::vector<float> v2h;
  try {v2h.resize(N2*N2*N2);} catch(std::exception& e) {myexception(e);}
#pragma omp parallel for shared(v2h,v1h)
  for (int ix=0; ix<N2; ++ix) {
    int ix0= 2*ix;
    int ixp=(2*ix+1  )%N;
    int ixm=(2*ix-1+N)%N;
    for (int iy=0; iy<N2; ++iy) {
      int iy0= 2*iy;
      int iyp=(2*iy+1  )%N;
      int iym=(2*iy-1+N)%N;
      for (int iz=0; iz<N2; ++iz) {
        int iz0= 2*iz;
        int izp=(2*iz+1  )%N;
        int izm=(2*iz-1+N)%N;
        v2h[N2*N2*ix+N2*iy+iz] = (8*v1h[N*N*ix0+N*iy0+iz0]+
                                  4*v1h[N*N*ixp+N*iy0+iz0]+
                                  4*v1h[N*N*ixm+N*iy0+iz0]+
                                  4*v1h[N*N*ix0+N*iyp+iz0]+
                                  4*v1h[N*N*ix0+N*iym+iz0]+
                                  4*v1h[N*N*ix0+N*iy0+izp]+
                                  4*v1h[N*N*ix0+N*iy0+izm]+
                                  2*v1h[N*N*ixp+N*iyp+iz0]+
                                  2*v1h[N*N*ixm+N*iyp+iz0]+
                                  2*v1h[N*N*ixp+N*iym+iz0]+
                                  2*v1h[N*N*ixm+N*iym+iz0]+
                                  2*v1h[N*N*ixp+N*iy0+izp]+
                                  2*v1h[N*N*ixm+N*iy0+izp]+
                                  2*v1h[N*N*ixp+N*iy0+izm]+
                                  2*v1h[N*N*ixm+N*iy0+izm]+
                                  2*v1h[N*N*ix0+N*iyp+izp]+
                                  2*v1h[N*N*ix0+N*iym+izp]+
                                  2*v1h[N*N*ix0+N*iyp+izm]+
                                  2*v1h[N*N*ix0+N*iym+izm]+
                                    v1h[N*N*ixp+N*iyp+izp]+
                                    v1h[N*N*ixm+N*iyp+izp]+
                                    v1h[N*N*ixp+N*iym+izp]+
                                    v1h[N*N*ixm+N*iym+izp]+
                                    v1h[N*N*ixp+N*iyp+izm]+
                                    v1h[N*N*ixm+N*iyp+izm]+
                                    v1h[N*N*ixp+N*iym+izm]+
                                    v1h[N*N*ixm+N*iym+izm])/64.0;
      }
    }
  }
  return(v2h);
}






void MultiGrid::vcycle(std::vector<float>& v,
                       const std::vector<float>& f, const int N) {
// Does one V-cycle, with a recursive strategy, replacing v in the process.
  jacobi(v,f,N);
  if (N>4 && N%2==0) {
    // Not at coarsest level -- recurse coarser.
    const int N2=N/2;
    std::vector<float> f2h=restrict(residual(v,f,N),N);
    // Make a vector of zeros as our first guess.
    std::vector<float> v2h;
    try {v2h.resize(N2*N2*N2);} catch(std::exception& e) {myexception(e);}
    // and recursively call ourself
    vcycle(v2h,f2h,N2);
    // take the residual and prolong it back to the finer grid
    std::vector<float> v1h=prolong(v2h,N2);
    // and correct our earlier guess.
    for (int nn=0; nn<v.size(); ++nn) v[nn] += v1h[nn];
  }
  jacobi(v,f,N);
}



void MultiGrid::print_error(const std::vector<float>& v,
                            const std::vector<float>& f, const int N) {
// For debugging purposes, prints an estimate of the residual.
  std::vector<float> r=MultiGrid::residual(v,f,N);
  float res1=0,res2=0,src1=0,src2=0;
#pragma omp parallel for shared(v,f) reduction(+:src1,src2,res1,res2)
  for (int nn=0; nn<r.size(); ++nn) {
    src1 += fabs(f[nn]);
    src2 += f[nn]*f[nn];
    res1 += fabs(r[nn]);
    res2 += r[nn]*r[nn];
  }
  src2 = sqrt(src2);
  res2 = sqrt(res2);
  std::cout<<"# Source   L1 norm is "<<std::fixed
           <<std::setprecision(7)<<src1<<std::endl;
  std::cout<<"# Residual L1 norm is "<<std::fixed
           <<std::setprecision(7)<<res1<<std::endl;
  std::cout<<"# Source   L2 norm is "<<std::fixed
           <<std::setprecision(7)<<src2<<std::endl;
  std::cout<<"# Residual L2 norm is "<<std::fixed
           <<std::setprecision(7)<<res2<<std::endl;
}




std::vector<float> MultiGrid::fmg(const std::vector<float>& f1h,
                                  const int N) {
// The full multigrid cycle, also done recursively.
  std::vector<float> v1h;
  if (N>4 && N%2==0) {
    // Recurse to a coarser grid.
    std::vector<float> f2h,v2h;
    f2h = restrict(f1h,N);
    v2h = fmg(f2h,N/2);
    v1h = prolong(v2h,N/2);
  }
  else {
    // Start with a guess of zero -- should be no memory issues at this
    // coarsest level.
    v1h.resize(N*N*N);
    std::fill(v1h.begin(),v1h.end(),0.0);
  }
  const int Niter=6;
  for (int iter=0; iter<Niter; ++iter) {
    vcycle(v1h,f1h,N);
  }
  if (N==Ng) print_error(v1h,f1h,N);
  return(v1h);
}
