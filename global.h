#ifndef	_GLOBAL_H_
#define	_GLOBAL_H_
#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>


// The size of the grid is Ng^3.
const	int Ng=1024;


struct particle {
  float pos[3],wt;
};

struct Box {
  float min[3],max[3];
  float ctr[3],L;
};


// Some helpful global variables.
extern	float	bias,beta;
extern	struct	Box box;

// and the method prototypes.

extern void    myexit(const int flag);
extern void    myexception(const std::exception& e);

#include "lcdm.h"

extern struct particle fill_particle(const double ra, const double dec,
                                     const double  z, const double  wt,
                                     const LCDM& lcdm);

extern std::vector<struct particle> read_data(const char fname[],
                                              const LCDM& lcdm);

extern void write_data(const std::vector<struct particle>& P,
                       const char fname[]);

extern void remap_pos(std::vector<struct particle>& D,
                      std::vector<struct particle>& R1,
                      std::vector<struct particle>& R2);


extern std::vector<float> make_grid(const std::vector<struct particle>& D,
                                    const std::vector<struct particle>& R,
                                    const float Rf);


extern void shift_obj(std::vector<struct particle>& P,
                      const std::vector<float>& phi);

#endif
