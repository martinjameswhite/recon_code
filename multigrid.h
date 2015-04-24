#include	<cmath>
#include	<vector>

class MultiGrid {
// A class to hold the main multigrid methods -- this is specialized to
// the modified Poisson equation described in the notes.  All of the methods
// are declared static in multigrid.h, since the class as currently
// configured holds no state information.
private:
static std::vector<float> residual(const std::vector<float>& v,
                                   const std::vector<float>& f,  const int N);
static std::vector<float> prolong(const std::vector<float>& v2h, const int N);
static std::vector<float> restrict(const std::vector<float>& v1h,const int N);
static void vcycle(std::vector<float>& v,
                   const std::vector<float>& f, const int N);
static void gauss_seidel(std::vector<float>& v,
                         const std::vector<float>& f, const int N);
static void jacobi(std::vector<float>& v,
                   const std::vector<float>& f, const int N);
static void print_error(const std::vector<float>& v,
                        const std::vector<float>& f,  const int N);
public:
  static std::vector<float> fmg(const std::vector<float>& f1h, const int N);
};
