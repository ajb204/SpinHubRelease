/*Legal_notice:

Copyright 1996-2001 by Soren M. Kristensen, Department of Chemistry,
University of Copenhagen, DK-2100 Copenhagen O, Denmark.  All rights
reserved.

This software and its related software and documentation is provided 'as
is' without express or implied warrenty. The Departmet of Chemistry, The
University of Copenhagen and the author make no warrenties as to any matter
whatsoever with respect to the program and the related software and
documentation. In particular, any and all warranties of merchantability and
fitness for any particular purpose are expressly excluded.

Legal_notice*/

#include "fuda_utils_perlman.H"

namespace FUDA
{
  // Some useful constants.
  const double pi = 3.14159265358979323844;
  const double log2 = log(2.0);

  // Some global utility functions.
  double calc_enorm(int n, double[]);
  double get_machine_double_eps();
  double calc_sd(int nfree, double enorm);
  int rand();
  void srand(unsigned int seed);
  unsigned int rand_seed();
  double rand_uniform();
  double rand_gauss(double mean, double sd);
  double gamma_func(double x);
  double log_gamma_func(double x);
  double chi2_distrib_dens(double x, unsigned int n);
  double chi2_distrib_dens2(double x, unsigned int n);
  double chi2_distrib_crit_old(double alpha, unsigned int n);
  

  // Function proto type.
  typedef int (*FUNC)(void *fs,
		      double p[],
		      int dp_flg[],
		      double dp[],
		      double *value);

  // Product function declaration.
  int nd_func(void *fs, double p[], int dp_flg[], double dp[], double *value);

  // nd_func product function structure (fs).
  class Nd_func_fs
    {
      // Number of functions (Dimension of N-dimensional function).
      unsigned int nfunc;
      // Number of parameters per dimension (array of dim nfunc).
      unsigned int *func_nparam;
      // Pointer to pointers to normalized 1D functions (array of dim nfunc).
      FUNC *func_ref;
      // pointer to array of function structures to pass to fncs (array
      // of dim nfunc).
      void **func_struct;
    public:

      Nd_func_fs(unsigned int n)
	{
	  nfunc = n;

	  // Allocate.
	  if (nfunc>0)
	    {
	      func_nparam = new unsigned int[nfunc];
	      func_ref = new FUNC[nfunc];
	      func_struct = new void*[nfunc];
	    }
	  else nfunc = 0;
	}

      ~Nd_func_fs()
	{
	  // Eventually deallocate.
	  if (nfunc>0)
	    {
	      delete[] func_nparam;
	      delete[] func_ref;
	      delete[] func_struct;
	    }
	}

      // Set data for function number ifunc.
      void set_func(unsigned int ifunc, FUNC func,
		    unsigned int np, void *fs)
	{
	  func_nparam[ifunc] = np;
	  func_ref[ifunc] = func;
	  func_struct[ifunc] = fs;
	}

      // Get data.
      unsigned int get_nfunc(){ return (nfunc); }
      unsigned int get_func_nparam(int i) { return (func_nparam[i]); }
      FUNC get_func_ref(unsigned int i) { return (func_ref[i]); }
      void *get_func_struct(int i) { return (func_struct[i]); }
  };

}
