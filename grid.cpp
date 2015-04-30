#include	<cmath>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"

#ifndef	BRUTEFORCESMOOTH
extern void fft_smooth(std::vector<float>& delta, const int N, const float Rf);
#endif


void	enclosing_box(const std::vector<struct particle>& R) {
// Uses the randoms to figure out an enclosing box -- stored as a
// global structure "box".
  // First work out the maximum and minimum coordinates.
  for (int idim=0; idim<3; ++idim) {
    box.min[idim]= 1e30;
    box.max[idim]=-1e30;
  }
  for (int nn=0; nn<R.size(); ++nn) {
    for (int idim=0; idim<3; ++idim) {
      if (R[nn].pos[idim]>box.max[idim])
        box.max[idim] = R[nn].pos[idim];
      if (R[nn].pos[idim]<box.min[idim])
        box.min[idim] = R[nn].pos[idim];
    }
  }
  // Now use this to define a center and a length, allowing room
  // for padding.
  const float padfactor=1.5;
  box.L=-1;
  for (int idim=0; idim<3; ++idim) {
    box.ctr[idim] = (box.min[idim]+box.max[idim])/2;
    if (box.L<box.max[idim]-box.min[idim])
      box.L = padfactor*(box.max[idim]-box.min[idim]);
  }
}





void	remap_pos(std::vector<struct particle>& D,
                  std::vector<struct particle>& R1,
                  std::vector<struct particle>& R2) {
// Remaps the positions into the box.
  enclosing_box(R1);
#pragma omp parallel for shared(D)
  for (int nn=0; nn<D.size(); ++nn)
    for (int idim=0; idim<3; ++idim)
      D[nn].pos[idim] = (D[nn].pos[idim]-box.ctr[idim])/box.L+0.5;
#pragma omp parallel for shared(R1)
  for (int nn=0; nn<R1.size(); ++nn)
    for (int idim=0; idim<3; ++idim)
      R1[nn].pos[idim] = (R1[nn].pos[idim]-box.ctr[idim])/box.L+0.5;
#pragma omp parallel for shared(R2)
  for (int nn=0; nn<R2.size(); ++nn)
    for (int idim=0; idim<3; ++idim)
      R2[nn].pos[idim] = (R2[nn].pos[idim]-box.ctr[idim])/box.L+0.5;
}




std::vector<float> cic_assign(const std::vector<struct particle>& P) {
// Generates a density by CIC assignment.
  std::vector<float> dg;
  try {dg.resize(Ng*Ng*Ng);} catch(std::exception& e) {myexception(e);}
  std::fill(dg.begin(),dg.end(),0.0); // Shouldn't need this.
  for (int nn=0; nn<P.size(); ++nn) {
    int         ii[3],ix,iy,iz;
    float       dx,dy,dz,wt;
    wt=P[nn].wt;
    ix=Ng*P[nn].pos[0]; ii[0]=(ix+1)%Ng; dx=Ng*P[nn].pos[0]-ix;
    iy=Ng*P[nn].pos[1]; ii[1]=(iy+1)%Ng; dy=Ng*P[nn].pos[1]-iy;
    iz=Ng*P[nn].pos[2]; ii[2]=(iz+1)%Ng; dz=Ng*P[nn].pos[2]-iz;
    if (ix<0 || ix>=Ng || iy<0 || iy>=Ng || iz<0 || iz>=Ng) {
      std::cout << "Index out of range:"<<std::endl
                << "ix="<<ix<<", iy="<<iy<<", iz="<<iz<<std::endl
                << "dx="<<dx<<", dy="<<dy<<", dz="<<dz<<std::endl
                << "pos["<<nn<<"]=("
                << P[nn].pos[0]<<","<<P[nn].pos[1]<<","<<P[nn].pos[2]<<")"
                << std::endl;
      std::cout.flush();
      myexit(1);
    }
    dg[(Ng*Ng* ix  +Ng* iy  + iz  )] += (1.0-dx)*(1.0-dy)*(1.0-dz)*wt;
    dg[(Ng*Ng*ii[0]+Ng* iy  + iz  )] +=      dx *(1.0-dy)*(1.0-dz)*wt;
    dg[(Ng*Ng* ix  +Ng*ii[1]+ iz  )] += (1.0-dx)*     dy *(1.0-dz)*wt;
    dg[(Ng*Ng* ix  +Ng* iy  +ii[2])] += (1.0-dx)*(1.0-dy)*     dz *wt;
    dg[(Ng*Ng*ii[0]+Ng*ii[1]+ iz  )] +=      dx *     dy *(1.0-dz)*wt;
    dg[(Ng*Ng*ii[0]+Ng* iy  +ii[2])] +=      dx *(1.0-dy)*     dz *wt;
    dg[(Ng*Ng* ix  +Ng*ii[1]+ii[2])] += (1.0-dx)*     dy *     dz *wt;
    dg[(Ng*Ng*ii[0]+Ng*ii[1]+ii[2])] +=      dx *     dy *     dz *wt;
  }
  return(dg);
}






void	smooth(std::vector<float>& delta, const int N, const float Rf) {
// Performs a Gaussian smoothing using brute-force convolution.
  // Take a copy of the data we're smoothing.
  std::vector<float> ss;
  try{
    ss.assign(delta.begin(),delta.end());
  } catch(std::exception& e) {myexception(e);}
  // Now set up the smoothing stencil.
  // The number of grid points to search: >=rmax*Rf.
  const float rmax= 2.5;
  const int   rad = (int)(rmax*Rf*N+1.0);
  const int   Ns  = 2*rad+1;
  const float fact= 0.5/(N*Rf)/(N*Rf);
  std::vector<float> kern;
  try{kern.resize(Ns*Ns*Ns);} catch(std::exception& e) {myexception(e);}
#pragma omp parallel for shared(kern)
  for (int dx=-rad; dx<=rad; ++dx) {
    for (int dy=-rad; dy<=rad; ++dy) {
      for (int dz=-rad; dz<=rad; ++dz) {
        int ii  = Ns*Ns*(rad+dx)+Ns*(dy+rad)+(dz+rad);
        float r2= fact*(dx*dx+dy*dy+dz*dz);
        if (r2<rmax*rmax/2.0)
          kern[ii]= exp(-r2);
        else
          kern[ii]= 0;
      }
    }
  }
#pragma omp parallel for shared(delta,ss,kern)
  for (int ix=0; ix<N; ++ix)
    for (int iy=0; iy<N; ++iy)
      for (int iz=0; iz<N; ++iz) {
        float sumd=0,sumw=0;
        for (int dx=-rad; dx<=rad; ++dx) {
          int iix = (ix+dx+N)%N;
          int i1  = N*N*iix;
          int j1  = Ns*Ns*(rad+dx);
          for (int dy=-rad; dy<=rad; ++dy) {
            int iiy = (iy+dy+N)%N;
            int i2  = N*iiy;
            int j2  = Ns*(rad+dy);
            for (int dz=-rad; dz<=rad; ++dz) {
              int iiz = (iz+dz+N)%N;
              int ii  = i1+i2+iiz;
              int jj  = j1+j2+(dz+rad);
              sumd += kern[jj]*ss[ii];
              sumw += kern[jj];
            }
          }
        }
        delta[N*N*ix+N*iy+iz] = sumd/(sumw+1e-30);
      }
}




void	print_stats(const char pre[], const std::vector<float>& delta) {
// Print some statistics about the density field.
  double avg=0,rms=0;
#pragma omp parallel for reduction(+:avg,rms)
  for(int nn=0; nn<delta.size();++nn) {
    avg+=delta[nn];
    rms+=delta[nn]*delta[nn];
  }
  avg /= delta.size();
  rms /= delta.size();
  double dmax=-1e30;
  for(int nn=0; nn<delta.size();++nn)
    if (delta[nn]>dmax) dmax=delta[nn];
  std::cout<<pre<<"delta is "
           <<std::fixed<<std::setw(8)<<std::setprecision(4)<<avg
           <<" +/- "<<sqrt(rms-avg*avg)
           <<", max="<<dmax<<std::endl;
}




std::vector<float> make_grid(const std::vector<struct particle>& D,
                             const std::vector<struct particle>& R,
                             const float Rf) {
// Bins the data and randoms onto a grid to produce the density
// contrast, delta = rho/rhobar - 1.
  // Assign the data and randoms to grids using CIC.
  std::vector<float> dg=cic_assign(D);
  std::vector<float> rg=cic_assign(R);
  // We remove any points which have too few randoms for a decent
  // density estimate -- this is "fishy", but it tames some of the
  // worst swings due to 1/eps factors.  Better would be an interpolation
  // or a pre-smoothing (or many more randoms).
  const float RanMin=0.75;
  int count=0;
  for (int nn=0; nn<rg.size(); ++nn)
    if (rg[nn]>0 && rg[nn]<RanMin) {
      dg[nn]=rg[nn]=0;
      count++;
    }
  std::cout<<"# Removed "<<count<<" low nr grid points from survey."<<std::endl;
  // We will need the total weights--here we want double precision.
  double dwt=0; for (int nn=0; nn<dg.size(); ++nn) dwt += dg[nn];
  double rwt=0; for (int nn=0; nn<rg.size(); ++nn) rwt += rg[nn];
  const float norm=rwt/(dwt+1e-30);
  // Now compute delta--use mean density where we have no randoms.
  std::vector<float> delta;
  try {delta.resize(Ng*Ng*Ng);} catch(std::exception& e) {myexception(e);}
#pragma omp parallel for shared(delta,dg,rg,bias)
  for (int ix=0; ix<Ng; ++ix)
    for (int iy=0; iy<Ng; ++iy)
      for (int iz=0; iz<Ng; ++iz) {
        int ii = Ng*Ng*ix+Ng*iy+iz;
        if (rg[ii]>0)
          delta[ii] = (dg[ii]/rg[ii]*norm - 1.0)/bias;
        else
          delta[ii] = 0;
      }
  // At this stage also remove the mean, so the source passed to the
  // multigrid routine is genuinely mean 0.  So as to not disturb the
  // padding regions, we only compute and subtract the mean for the
  // regions with delta!=0.
  double avg=0,cnt=0;
  for (int nn=0; nn<delta.size(); ++nn)
    if (delta[nn]!=0) {
      avg += delta[nn];
      cnt += 1;
    }
  avg /= cnt+1e-30;
  for (int nn=0; nn<delta.size(); ++nn)
    if (delta[nn]!=0)
      delta[nn] -= avg;
  print_stats("# Pre-smoothing:  ",delta);
  // Smooth this grid--pass smoothing length in box units.
  // This is the most time consuming routine by far (apart from ascii I/O)
  // so we also provide an FFT-based option which is much faster.
#ifdef	BRUTEFORCESMOOTH
  smooth(delta,Ng,Rf/box.L);
#else
  fft_smooth(delta,Ng,Rf/box.L);
#endif
  print_stats("# Post-smoothing: ",delta);
  return(delta);
}
