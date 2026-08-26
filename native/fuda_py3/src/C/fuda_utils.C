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

#include <cstdlib>
#include <ctime>
#include "fuda_classes.H"
#include "fuda_utils.H"
#include "minpack.H"

///////////////////////////////////////////////////////////////
// UF namespace utility functions.
///////////////////////////////////////////////////////////////

// Calculate euclidian norm of double array.
double FUDA::calc_enorm(int n, double x[])
{
  // Wrapper for minpack enorm_ function.
  return (enorm_(&n,x));
}

// Return machine (double) precision.
double FUDA::get_machine_double_eps()
{
  // Wrapper for minpack dpmpar_ function.
  const int i=1;
  return (dpmpar_(&i));
}

// Calculate overall standard deviation.
double FUDA::calc_sd(int nfree, double enorm)
{
  return (enorm/sqrt(nfree));
}

// Product function. Creates a new function of one or more other functions.
int FUDA::nd_func(void *fs, double p[], int dp_flg[], double dp[], double *value) 
{
  // Type-cast reference to get nd_fs structure.
  FUDA::Nd_func_fs& ndfs = *(FUDA::Nd_func_fs*)fs;
  
  // References to Scale factor (intens) and return function value.
  double& I = p[0];      // intensity  is always the first parameter.
  double& y = *value;     // return value.
  
  // Keepers for local storage.
  static unsigned int alloc_nfunc = 0;
  static double *norm = NULL;
  
  // If nfunc is zero or less we deallocate storage and return immediately.
  unsigned int nfunc = ndfs.get_nfunc();
  if (nfunc<=0)
    {
      if (norm!=NULL) delete norm;
      return (0);
    }
  
  // Make sure we have space allocated to normalisation factors.
  if (nfunc>alloc_nfunc)
    {
      // Maybe deallocate previous storage.
      if (norm!=NULL) delete norm;
      
      // Allocate normalization factor array.
      alloc_nfunc = nfunc;
      norm = new double[alloc_nfunc];
    }
  
  /* Calculate the normalised function contribution for each
     function. If dp_flg[] is set we also calculate corresponding
     partial derrivatives.  First, we initialize the point value
     and the normalisation factors for each function. Second, we
     calculate the point value and the normalisation constants.
     The calculated value of norm[ifunc] is identical to the point
     value divided by the ifunc'th function contribution. The
     calculation is robust in cased where function contributions
     in one or more dimensions are zero. */
  
  /* Initialize normalized return value. Later we multiply with I
     to get final return value */
  double y_norm = 1.0;
  
  // Initialize normalization factors.
  for (unsigned int ifunc=0; ifunc<nfunc; ifunc++)
    norm[ifunc] = y_norm;
  
  // Loop over functions.
  unsigned int offset = 1;  /* this offsets for the intensity param
			       which is the first. */
  for (unsigned int ifunc=0; ifunc<nfunc; ifunc++)
    {
      // Calculate function values for each dimension.
      double dvalue;
      int rtn = (ndfs.get_func_ref(ifunc))
	(ndfs.get_func_struct(ifunc),
	 &p[offset],
	 &dp_flg[offset],
	 &dp[offset],
	 &dvalue);
      
      // Check return value and abort on exception.
      if (rtn!=0) return(rtn);
      
      // Update y_norm.
      y_norm *= dvalue;
      
      // Update normalization factors.
      for (unsigned int jfunc=0; jfunc<nfunc; jfunc++)
	if (jfunc!=ifunc) norm[jfunc] *= dvalue; 
      
      // Count up offset.
      offset += ndfs.get_func_nparam(ifunc);
    }
  
  // We now turn to calculating the derrivatives.
  
  // dy/dI is simply y_norm. 
  if (dp_flg[0]) dp[0] = y_norm;
  
  /* The rest have already been calculated by the normalized
     function, but needs scaling with the normalization constants
     for the respective dimensions. */
  offset = 1;  /* This offsets for the intensity param which is
		  the first. */
  for (unsigned int ifunc=0; ifunc<nfunc; ifunc++)
    {
      for (unsigned int iparam=0;
	   iparam<ndfs.get_func_nparam(ifunc);
	   iparam++)
	if (dp_flg[offset+iparam])
	  dp[offset+iparam] *= I*norm[ifunc];
      
      // Count up offset.
      offset += ndfs.get_func_nparam(ifunc);
    }
  
  // Finally, we calculate the return value.
  y = I*y_norm;
  
  return (0);
}


/* Random number */
int FUDA::rand()
{
  return (::rand());
}


/* Set random number generator seed */
void FUDA::srand(unsigned int seed)
{
  ::srand(seed);
}


/* Get a random number from the real-time clock to use as a seed number. */
unsigned int FUDA::rand_seed()
{
  unsigned int t = time(0);
  return(t);
}


/* Uniform deviate [0;1[ */
double FUDA::rand_uniform()
{
  return (double(::rand())/(1.0+double(RAND_MAX)));
}


/* Normal distributed noise with mean and sd. This is a hack from
   slatec routine rgauss.f */
double FUDA::rand_gauss(double mean, double sd)
{
  int repeat = 6;
  double val = -repeat;
  for(int i=0; i<2*repeat; i++)
    val += double(::rand())/double(RAND_MAX);
  return (mean+sd*val);
}


double FUDA::gamma_func(double x)
{
  /* Returns the gamma function of x according to Stirlings formula
     including 3 terms of Stirlings asymtothic series. */

  double gamma =
    sqrt(2.0*pi*(x-1.0))
    *pow(x-1.0,x-1.0)
    *exp(-(x-1.0))
    *(1.0+1.0/12.0/(x-1.0)
      +1.0/228.0/pow(x-1.0,2)
      -139.0/51840.0/pow(x-1.0,3));
    return (gamma);  
}


double FUDA::log_gamma_func(double x)
{
  /* Returns the natural lograrithm of the gamma function of x
     according to Stirlings formula including 3 terms of Stirlings
     asymtothic series. */

  double y = x-1.0;
  double log_gamma =
    0.5*log(2.0*pi*y) + y*log(y) - y 
    + log(1.0 + 1.0/(12.0*y) + 1.0/(288.0*y*y) - 139.0/(51840.0*y*y*y));

  return (log_gamma);  
}


double FUDA::chi2_distrib_dens(double x, unsigned int n)
{
  /* Returns the density of the chi-square distribution for a given
     chi2 value x and the number of degrees of freedom. */

  static unsigned int nprev=0;
  static double c_prev=0;
  double c;
  
  // Calculate normalization constant or reuse previous value.
  if (n==nprev) c = c_prev;
  else 
    {
      c = 1.0/(pow(2.0,n/2.0) * FUDA::gamma_func(n/2.0));
      nprev = n;
      c_prev = c;
    }

  // Calculate and return density.
  double dens = c * pow(x,(n-2)/2.0) * exp(-x/2.0);
  
  return (dens);
}

double FUDA::chi2_distrib_dens2(double x, unsigned int n)
{
  /* Returns the density of the chi-square distribution for a given
     chi2 value x and the number of degrees of freedom. */
  
  static unsigned int nprev=0;
  static double log_gamma=0;
  
  // Calculate log_gamma or reuse previous value.
  if (n!=nprev)
    {
      log_gamma = FUDA::log_gamma_func(n/2.0);
      nprev = n;
    }

  // Calculate log to density.
  double log_dens = (n-2.0)/2.0*log(x) - x/2.0
    - n/2.0*log(2.0) - log_gamma;
  
  // Calculate and return density.
  return (exp(log_dens));
}

double FUDA::chi2_distrib_crit_old (double alpha, unsigned int n)
{
  /* Returns the alpha critical value of the chi2 distribution for n
     degrees of freedom. */

  double target = 1.0-alpha;
  double delta_chi2 = 0.01;

  double sum = 0.0;
  double chi2 = -delta_chi2/2.0;
  while (sum<target)
    {
      chi2 += delta_chi2;
      sum += delta_chi2*FUDA::chi2_distrib_dens(chi2,n);
    }
  chi2 += delta_chi2/2.0;

  return (chi2);
}
