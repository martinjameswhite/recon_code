#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"









struct particle fill_particle(const double ra, const double dec,
                              const double  z, const double  wt,
                              const LCDM& lcdm) {
  // Converts from ra, dec and z to 3D Cartesian position.
  struct particle curp;
  if (z<0.0 || z>2.0) {
    std::cout << "z="<<z<<" out of range."<<std::endl;
    myexit(1);
  }
  if (dec<-90.0 || dec>90.0) {
    std::cout << "DEC="<<dec<<" out of range."<<std::endl;
    myexit(1);
  }
  double r    =  lcdm.chiz(z);
  double theta= (90.0 - dec)*M_PI/180.;
  double phi  = (ra        )*M_PI/180.;
  curp.pos[0] = r*sin(theta)*cos(phi);
  curp.pos[1] = r*sin(theta)*sin(phi);
  curp.pos[2] = r*cos(theta);
  curp.wt     = wt;
  return(curp);
}







std::vector<struct particle> read_data(const char fname[], const LCDM& lcdm) {
  // Load the data from an ascii file.
  std::ifstream fs(fname);
  if (!fs) {
    std::cerr<<"Unable to open "<<fname<<" for reading."<<std::endl;
    myexit(1);
  }
  std::string buf;
  do {        // Skip any preamble comment lines.
    getline(fs,buf);
  } while(!fs.eof() && buf[0]=='#');
  std::vector<struct particle> P;
  try {P.reserve(10000000);} catch(std::exception& e) {myexception(e);}
  while (!fs.eof()) {
    double ra,dec,z,wt;
#ifdef  READWEIGHT
    std::istringstream(buf) >> ra >> dec >> z >> wt;
#else
    std::istringstream(buf) >> ra >> dec >> z;
    wt = 1.0;
#endif
    struct particle curp = fill_particle(ra,dec,z,wt,lcdm);
    try {
      P.push_back(curp);
    } catch(std::exception& e) {myexception(e);}
    getline(fs,buf);
  }
  fs.close();
  return(P);
}



void	write_data(const std::vector<struct particle>& P, const char fname[]) {
// Writes the normalized positions to a file.
  std::ofstream ofs(fname,std::ios::trunc);
  if (!ofs) {
    std::cerr<<"Unable to open "<<fname<<" for writing."<<std::endl;
    myexit(1);
  }
  for (int nn=0; nn<P.size(); ++nn)
    ofs<<std::fixed<<std::setw(15)<<std::setprecision(4)
       <<box.ctr[0]+box.L*(P[nn].pos[0]-0.5)
       <<std::fixed<<std::setw(15)<<std::setprecision(4)
       <<box.ctr[1]+box.L*(P[nn].pos[1]-0.5)
       <<std::fixed<<std::setw(15)<<std::setprecision(4)
       <<box.ctr[2]+box.L*(P[nn].pos[2]-0.5)
       <<std::fixed<<std::setw(15)<<std::setprecision(4)
       <<P[nn].wt<<std::endl;
  ofs.close();
}
