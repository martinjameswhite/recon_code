#include	<cmath>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"
#include	"multigrid.h"

// A code to perform reconstruction, using the lowest order algorithm.
//
// Author:	Martin White	(UCB/LBNL)
// Written:	20-Apr-2015
// Modified:	20-Apr-2015
//
//
// This code uses a standard multigrid technique to compute the
// displacements given an observed density field and then move the
// objects and "randoms" back along the displacement vector.
//
// Want to get rid of continual allocation and deallocation of arrays.
// I think we just need two arrays at each level, could call them A & B,
// or V1 and V2.
//


// Global variables.
float	bias=1,beta=0;
struct	Box box;



void    myexit(const int flag)
{
  std::cout.flush();
  std::cerr.flush();
  exit(flag);
}


void    myexception(const std::exception& e)
{
  std::cout<<": Exception: "<<e.what()<<std::endl;
  std::cout.flush();
  std::cerr.flush();
  exit(1);
}









int	main(int argc, char **argv)
{
  float Rf;
  if (argc!=7) {
    std::cout<<"Usage: recon <data-file> <random-file> <random-file>"
             <<" <bias> <f-growth> <R-filter>"<<std::endl;
    myexit(1);
  }
  bias = atof(argv[4]);		// Sets global variable.
  beta = atof(argv[5])/bias;	// Sets global variable.
  Rf   = atof(argv[6]);
  LCDM lcdm(0.30);

#ifdef	TESTMG
  // Make a cosine wave (only in x-direction) and solve for it with
  // beta=0.
  box.L=1.25;  box.ctr[0]=box.ctr[1]=box.ctr[2]=0.6;
  const float dt=2*M_PI/Ng;
  std::vector<float> src(Ng*Ng*Ng);
  for (int ix=0; ix<Ng; ++ix)
    for (int iy=0; iy<Ng; ++iy)
      for (int iz=0; iz<Ng; ++iz)
        src[Ng*Ng*ix+Ng*iy+iz] = cos(ix*dt);
  bias = 1; beta = 0;
  std::vector<float> ans = MultiGrid::fmg(src,Ng);
  for (int ix=0; ix<Ng; ++ix)
    std::cout<<std::fixed<<std::setprecision(6)
             <<ix<<" "<<cos(ix*dt)/(2*M_PI*2*M_PI)
             <<" "<<ans[Ng*Ng*ix+Ng*0   +Ng/4]
             <<" "<<ans[Ng*Ng*ix+Ng*Ng/2+0]
             <<" "<<ans[Ng*Ng*ix+Ng*0   +Ng/2]<<std::endl;
  return(0);
#endif

  // Read the data and figure out the 3D positions within
  // and enclosing box.
  std::vector<struct particle> D = read_data(argv[1],lcdm);
  std::vector<struct particle> R1= read_data(argv[2],lcdm);
  std::vector<struct particle> R2= read_data(argv[3],lcdm);
  std::cout<<"# Read "<<std::setw(10)<<D.size()
           <<" objects from "<<argv[1]<<std::endl;
  std::cout<<"# Read "<<std::setw(10)<<R1.size()
           <<" randoms from "<<argv[2]<<std::endl;
  std::cout<<"# Read "<<std::setw(10)<<R2.size()
           <<" randoms from "<<argv[3]<<std::endl;
  remap_pos(D,R1,R2);
  std::cout<<"# Enclosing survey in a box of side "<<box.L<<" Mpc/h."
           <<std::endl;
  std::cout<<"# Grid/mesh size is "<<box.L/Ng<<" Mpc/h"
           <<" and filter scale is "<<Rf<<" Mpc/h."
           <<std::endl;

  write_data(D ,"data_raw.xyzw");
  write_data(R2,"rand_raw.xyzw");
  
  // Make the density (contrast) grid and solve for the
  // displacement potential, phi.
  std::vector<float> delta = make_grid(D,R1,Rf);
  std::vector<float> phi   = MultiGrid::fmg(delta,Ng);

  // Shift the particles and randoms back -- if you want to not enhance
  // the line-of-sight shift for the randoms you need to change beta before
  // calling shift_obj.
  shift_obj(D ,phi);
  shift_obj(R2,phi);

  write_data(D ,"data_rec.xyzw");
  write_data(R2,"rand_rec.xyzw");

  return(0);
}
