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
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>


/************************************************************************/
/* MinpackLM. Class implementing MINPACK lmstr Levenberg-Marquardt
   minimization. */
class MinpackLM
{
  /* The stack is used for saving a reference to the MinpackLM object
     before calling lmstr_. The function MinpackLM::get_object_ref()
     is used within lmstr_eval_ to get a pointer to the MinpackLM
     object. */
  static std::list<MinpackLM*> stack;

  Fuda *fuda; // fuda object to which this MinpackLM object refers.
  bool initialized; // True if initialize() has been called.
  bool minimized; /* True if minimize() has been called since the last
		     call to initialize(). */
  unsigned int eval_mod_count; /* MinpackLM structure is in sync with
                                  fuda->eval struc if init is true
                                  and eval_mod_count is equal to
                                  fuda->eval_get_mod_count(). the
                                  get_sync() function returns true if
                                  the MinpackLM structure is in
                                  sync. */

  // Numerical derivative flag and epsilon.
  bool numderiv;  // Flag specifying numerical differntiation.
  double numderiv_eps;  // Numerical differentiation epsilon.
  
  // Parameters controlling lmstr.
  double tol;     // lmstr1 tol (see lm_set_tol(tol)).
  double ftol;    // lmstr exit tolerance.
  double xtol;    // lmstr exit tolerance.
  double gtol;    // lmstr exit tolerance.
  double ctol;    // lmstr covar tolerance.
  int maxfev;     // lmstr maximum function evaluation with iflag=1.
  double factor;  // lmstr initial step bound factor.
  int nprint;     // Print status every nprint iteration.
  bool scale_covar; // Scale covariance matrix with sd*sd.

  // Minimization output values.
  double enorm;   // lmstr euclidian norm after fit.
  double sd;      // lmstr overall standard deviation after fit.
  
  // Storage allocation related variables.
  bool allocated; // Arrays have previously been allocated.
  unsigned int m_max; // Current allocation of number of points.    
  unsigned int n_max; /* Current allocation of number of dependent params.
			 m_max and n_max determine the memory allocation 
			 related values: ldfjac=n_max; lwa=5*n_max+m_max; plus
			 allocation size of x, fvec, fjac, ipvt and wa. */ 
  unsigned int lwa_max; // Currently allocated size of work space. */
  
  // Structure holding all lmstr arguments.
  struct Arg
  {
    int m;        // Number of data points.
    int n;        // Number of dependent parameters.
    double *x;    // Initial parameters on entry, final on exit.
    double *fvec; // Final calculated data point values.
    double *fjac; // n*n twodimensional array with encoded Jacobian.
    int ldfjac;   // Leading dimension of fjac, must be >= n.
    double ftol;  // Tolerance in sum of squares.
    double xtol;  // Tolerance in x vector.
    double gtol;  /* Tolerance in orthogonality between function vector and 
		     columns of the jacobian */
    int maxfev;   // Maximum number of function evaluations with iflag=1.
    double *diag; // Parameter scaling factors if mode=2.
    int mode;     // Parameter scaling mode: 1=internal scaling, 2=use diag.
    double factor;// Initial step bound factor. 100 is recommended.
    int nprint;   // Print status for every nprint iterations.
    int info;     // Output status.
    int nfev;     // Output number of function evaluations.
    int njev;     // Output number of calls to calculate derrivatives.
    int *ipvt;    // int[n].
    double *qtf;  // double[n].
    double *wa1;  // double[n].
    double *wa2;  // double[n].
    double *wa3;  // double[n].
    double *wa4;  // double[m].
  } arg;

  // Private functions.
  void push_object_ref(MinpackLM *lm);
  void pop_object_ref();
  void set_initialized();
  void clear_initialized();
  void set_minimized();
  void clear_minimized();

 public:

  MinpackLM(Fuda *fuda);

  static MinpackLM *top_object_ref();
  

  Fuda *get_fuda();
  bool is_sync();
  void initialize();
  bool is_initialized();
  void minimize();
  bool is_minimized();
  bool fit_is_ok();
  bool is_converged();
  double get_enorm();
  double get_sd();
  double get_value(unsigned int i);
  double get_esd(unsigned int i);
  double get_covar(unsigned int i, unsigned int j);

  void set_tol(double tol);
  double get_tol();
  void set_ftol(double tol);
  double get_ftol();
  void set_xtol(double tol);
  double get_xtol();
  void set_gtol(double tol);
  double get_gtol();
  void set_ctol(double tol);
  double get_ctol();
  void set_maxfev(int maxfev);
  int get_maxfev();
  void set_factor(double factor);
  double get_factor();
  void set_scale_covar(bool scale_covar);
  bool get_scale_covar();
  void set_numderiv(bool numderiv);
  bool get_numderiv();
  void set_numderiv_eps(double eps);
  double get_numderiv_eps();
  void set_nprint(int nprint);
  int get_nprint();
  int get_info();
  int get_nfev();
  int get_njev();

  void print();
};


/************************************************************************/
/* Uferr. Uf error signalling namespace holding all exceptions created by
   uf routines. */
namespace Uferr
{
  struct LmTolInvalid {
    double tol;
    std::string name;
    LmTolInvalid(double t, std::string n) : tol(t), name(n) {}
  };
  struct LmTolsInvalid {};
  struct LmMaxfevInvalid {
    int maxfev;
    LmMaxfevInvalid(int maxf) : maxfev(maxf) {}
  };
  struct LmFactorInvalid {
    double factor;
    LmFactorInvalid(double f) : factor(f) {}
  };
  struct LmEpsfcnInvalid {
    double value;
    LmEpsfcnInvalid(double eps) : value(eps) {}
  };
  struct LmNoFit {};
  struct LmFitNotConverged {};
  struct LmParamIndexInvalid {
    unsigned int i;
    LmParamIndexInvalid(unsigned int a_i) : i(a_i) {}
  };

}


extern "C" void lmstr_eval_(int *m,
			    int *n,
			    double x[],
			    double fvec[],
			    double fjrow[],
			    int *iflag);

