#ifndef	BRUTEFORCESMOOTH
#include	<cmath>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"fftw3.h"
#include	"global.h"





void	fft_smooth(std::vector<float>& delta, const int N, const float Rf) {
// Performs a Gaussian smoothing using the FFTW library (v3), assumed to be
// installed already.  Rf is assumed to be in box units.
  // Make temporary vectors.  FFTW uses double precision.
  std::vector<double> vi,vout;
  try{
    vi.resize(N*N*N);
    vout.resize(2*N*N*(N/2+1));
  } catch(std::exception& e) {myexception(e);}
  fftw_complex *vo = (fftw_complex *)&vout[0];
#pragma omp parallel for shared(delta,vi)
  for (int nn=0; nn<delta.size(); ++nn) vi[nn]=delta[nn];
  // Generate the FFTW plan files.
  fftw_plan fplan = fftw_plan_dft_r2c_3d(N,N,N,&vi[0],&vo[0],FFTW_ESTIMATE);
  fftw_plan iplan = fftw_plan_dft_c2r_3d(N,N,N,&vo[0],&vi[0],FFTW_ESTIMATE);
  fftw_execute(fplan);
  // Now multiply by the smoothing filter.
  const double fact=0.5*Rf*Rf*(2*M_PI)*(2*M_PI);
#pragma omp parallel for shared(vo)
  for (int ix=0; ix<N; ++ix) {
    int iix = (ix<=N/2)?ix:ix-N;
    for (int iy=0; iy<N; ++iy) {
      int iiy = (iy<=N/2)?iy:iy-N;
      for (int iz=0; iz<N/2+1; ++iz) {
        int iiz = (iz<=N/2)?iz:iz-N;
        int ip  = N*(N/2+1)*ix+(N/2+1)*iy+iz;
        double smth = exp(-fact*(iix*iix+iiy*iiy+iiz*iiz));
        vo[ip][0] *= smth;
        vo[ip][1] *= smth;
      }
    }
  }
  vo[0][0]=vo[0][1]=0;	// Set the mean to zero.
  fftw_execute(iplan);
#pragma omp parallel for shared(delta,vi)
  for (int nn=0; nn<delta.size(); ++nn) delta[nn]=vi[nn]/N/N/N;
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(iplan);
}
#endif
