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
#include <stack>
#include <vector>
#include <ctime>
#include <cmath>
#include "minpack.H"
#include "fuda_classes.H"
#include "fuda_utils.H"
#include "fuda_MinpackLM.H"

// Static class members.
std::list<MinpackLM*> MinpackLM::stack;

// Functions returning the largest/smallest value of any type.
template <class T> const T& max(const T& a, const T& b)
{
  return (a<b) ? b : a;
}

// Functions returning the largest/smallest value of any type.
template <class T> const T& min(const T& a, const T& b)
{
  return (a<b) ? a : b;
}

// Constructor for MinpackLM data structure.
MinpackLM::MinpackLM(Fuda *fuda_ref)
{
  // Save fuda referenc.
  fuda = fuda_ref;

  // Initialize flags and control variables.
  clear_initialized();
  clear_minimized();
  eval_mod_count = 0;
  allocated = 0;

  // Some default settings.
  numderiv = 0; 
  numderiv_eps = FUDA::get_machine_double_eps(); 
  maxfev = 1000; 
  nprint = 1;
  factor = 100.0;
  ctol = 0.0;
  scale_covar= 1;
}


void MinpackLM::push_object_ref(MinpackLM *lm)
{
  stack.push_back(lm);
}


void MinpackLM::pop_object_ref()
{
  stack.pop_back();
}


MinpackLM *MinpackLM::top_object_ref()
{
  return (stack.back());
}


Fuda *MinpackLM::get_fuda()
{
  return (fuda);
}


void MinpackLM::set_initialized()
{
  initialized = 1;
}


void MinpackLM::clear_initialized()
{
  initialized = 0;
}


void MinpackLM::set_minimized()
{
  minimized = 1;
}


void MinpackLM::clear_minimized()
{
  minimized = 0;
}


bool MinpackLM::is_initialized()
{
  return(initialized);
}


bool MinpackLM::is_minimized()
{
  return(initialized && minimized);
}


bool MinpackLM::is_sync()
{
  /* The MinpackLM object is in sync with fuda->eval structure if
     initialize() has been called and fuda->eval_get_mod_count() is
     unchanged */
  return(initialized && (eval_mod_count==fuda->eval_get_mod_count()));
}


// Set tolerance for lm minimization routine.
void MinpackLM::set_tol(double tolerance)
{ 
  // Is the value valid.
  if (tolerance<=0.0 && gtol<=0.0) throw Uferr::LmTolsInvalid();

  // This corresponds to the minpack lmstr1 behaviour.
  tol = tolerance;
  ftol = tolerance;
  xtol = tolerance;
  gtol = 0.0;
}


double MinpackLM::get_tol() 
{
  return(tol);
}


void MinpackLM::set_ftol(double tolerance)
{
  ftol = tolerance;
}


double MinpackLM::get_ftol()
{
  return(ftol);
}


void MinpackLM::set_xtol(double tolerance)
{
  xtol = tolerance;
}


double MinpackLM::get_xtol()
{
  return(xtol);
}


void MinpackLM::set_gtol(double tolerance) 
{
  gtol = tolerance;
}


double MinpackLM::get_gtol()
{
  return(gtol);
}


void MinpackLM::set_ctol(double tolerance) 
{
  ctol = tolerance;
}


double MinpackLM::get_ctol() 
{
  return(ctol);
}


void MinpackLM::set_maxfev(int max_func_eval)
{
  maxfev = max_func_eval; 
}


int MinpackLM::get_maxfev()
{
  return(maxfev);
}


void MinpackLM::set_factor(double a_factor)
{
  factor = a_factor; 
}


double MinpackLM::get_factor()
{
  return(factor);
}


void MinpackLM::set_numderiv_eps(double eps)
{
  // Get machine precision. Don't set it to less than the machine eps.
  numderiv_eps = max(fabs(eps),FUDA::get_machine_double_eps());
}


double MinpackLM::get_numderiv_eps() 
{
  return(numderiv_eps);
}


void MinpackLM::set_numderiv(bool a_numderiv)
{
  numderiv = a_numderiv; 
}


bool MinpackLM::get_numderiv()
{
  return(numderiv);
}


void MinpackLM::set_scale_covar(bool a_scale_covar)
{
  scale_covar = a_scale_covar; 
}


bool MinpackLM::get_scale_covar()
{
  return(scale_covar);
}


void MinpackLM::set_nprint(int a_nprint)
{
  nprint = a_nprint;
}


int MinpackLM::get_nprint() 
{
  return(nprint);
}


int MinpackLM::get_info()
{
  return(arg.info);
}


int MinpackLM::get_nfev()
{
  return(arg.nfev);
}


int MinpackLM::get_njev()
{
  return(arg.njev);
}


double MinpackLM::get_sd()
{
  return(sd);
}


double MinpackLM::get_enorm()
{
  return(enorm);
}


// Return true if the fit went OK but not necessarily converged.
bool MinpackLM::fit_is_ok()
{
  if (is_sync() && is_minimized() && arg.info>=1)
    return (1);
  else
    return (0);
}


// Return true if the fit converged.
bool MinpackLM::is_converged()
{
  if (fit_is_ok() && arg.info<=4)
    return (1);
  else
    return (0);
}


// Initialize/setup lm struc. which is needed for an lmstr minimization.
void MinpackLM::initialize()
{
  // Function for initializing the MinpackLM object for a minimization.

  // If eval and MinpackLM object is in sync we do nothing.
  if (is_sync() && fuda->eval_is_sync()) return;

  /* The fuda->eval structure must be in rsync and the Func and Parm
     must be sync. */
  if (!fuda->eval_is_rsync() ||
      !fuda->func_is_sync() ||
      !fuda->param_is_sync()) fuda->eval_init();

  // Get and calculate allocation requirements.
  unsigned int m_new = fuda->eval_get_ndata(); // Number of data points.
  unsigned int n_new = fuda->eval_get_nfree(); // Number of free parameters.

  // Allocate or reallocate arrays.

  // Number of data points (number of functions in minpack terminology).
  if (m_max<m_new || !allocated)
    {
      // Eventually deallocate.
      if (allocated)
	{
	  delete(arg.fvec);
	  delete(arg.wa4);
	}
      
      // Allocate.
      arg.fvec = new double[m_new];
      arg.wa4 = new double[m_new];
      
      // Set new allocation size.
      m_max = m_new;
    }
  
  // Number of parameters.
  if (n_max<n_new || !allocated)
    {
      // Eventually deallocate.
      if (allocated)
	{
	  delete(arg.x);
	  delete(arg.diag);
	  delete(arg.qtf);
	  delete(arg.wa1);
	  delete(arg.wa2);
	  delete(arg.wa3);
	  delete(arg.fjac);
	  delete(arg.ipvt);
	}
      
      // Allocate.
      arg.x    = new double[n_new];
      arg.diag = new double[n_new];
      arg.qtf  = new double[n_new];
      arg.wa1  = new double[n_new];
      arg.wa2  = new double[n_new];
      arg.wa3  = new double[n_new];
      arg.fjac = new double[n_new*n_new];
      arg.ipvt = new int[n_new];
      
      // Set new allocation size.
      n_max = n_new;
    }

  // Set allocated flag.
  allocated = 1;

  // Then we setup the arg structure.

  // Dimensions.
  if (n_new<1) throw Uferr::NoFreeParam();
  if (n_new>m_new) throw Uferr::NparamGtNdata(n_new, m_new);
  arg.m = m_new;
  arg.n = n_new;

  // Initial paramters.
  for (unsigned int i=0; i< fuda->eval_get_nfree(); i++)
    arg.x[i] = fuda->eval_get_free(i)->get_value();

  // Leading dimension size of fjac.
  arg.ldfjac = arg.n;

  // Tolerances.
  if (ftol < 0.0) throw Uferr::LmTolInvalid(ftol,"ftol");
  if (xtol < 0.0) throw Uferr::LmTolInvalid(xtol,"xtol");
  if (gtol < 0.0) throw Uferr::LmTolInvalid(gtol,"gtol"); 
  if (ftol==0.0 && ftol==0.0 && ftol==0.0)
    throw Uferr::LmTolsInvalid();
  arg.ftol = ftol;
  arg.xtol = xtol;
  arg.gtol = gtol;

  // Check ctol used for calculation of the covariance matrix.
  if (ctol < 0.0) throw Uferr::LmTolInvalid(ctol,"ctol"); 

  // Max func. eval.
  if (maxfev<1) throw Uferr::LmMaxfevInvalid(maxfev);
  arg.maxfev = maxfev;

  // Scaling mode is hard wired.
  arg.mode = 1;

  // Initial step bound factor.
  if (factor<=0.0) throw Uferr::LmFactorInvalid(factor);
  arg.factor = factor;

  // Initial nprint reporting divider.
  arg.nprint = nprint;

  // Get and save fuda.eval mod_count.
  eval_mod_count = fuda->eval_get_mod_count();

  // Finally, set initialized flag and clear minimized flag.
  initialized = 1;
  minimized = 0;
}


// Perform an Levenberg-Maquardt minimization.
void MinpackLM::minimize()
{
  /* Function for peforming a lmstr (minpack Levenberg-Marquart)
     minimization. */

  // The lm structure must be in sync.
  initialize();

  // push MinpackLM object reference to stack.
  push_object_ref(this);

  // Call the minpack minimizer.
  lmstr_(lmstr_eval_,
	 &arg.m,
	 &arg.n,
	 arg.x,
	 arg.fvec,
	 arg.fjac,
	 &arg.ldfjac,
	 &arg.ftol,
	 &arg.xtol,
	 &arg.gtol,
	 &arg.maxfev,
	 arg.diag,
	 &arg.mode,
	 &arg.factor,
	 &arg.nprint,
	 &arg.info,
	 &arg.nfev,
	 &arg.njev,
	 arg.ipvt,
	 arg.qtf,
	 arg.wa1,
	 arg.wa2,
	 arg.wa3,
	 arg.wa4);

  // Pop MinpackLM object reference from stack.
  pop_object_ref();

  // If the fit succeeded, we proceed.
  if (arg.info>=1)
    {
      // Calculate final euclidian norm and standard deviation.
      enorm = FUDA::calc_enorm(arg.m, arg.fvec);
      // NB! nfree=m-n, but if nfree=0, we use nfree=1.
      sd = FUDA::calc_sd(max(arg.m-arg.n,1), enorm);


      // If the fit converged, we calculate the covariance matrix.
      if (arg.info<=4)
	{	  
	  // Calculate the covariance matrix - it replaces arg.fjac.
	  covar_(&arg.n,
		 arg.fjac,
		 &arg.ldfjac,
		 arg.ipvt,
		 &ctol,
		 arg.wa1);
	  
	  /* If the covariance matrix is scaled by sd^2, the
	     uncertaities are no longer governed by the experimental
	     uncertainties, but by the size of sd. Thus, with scaling
	     turne on, a perfect fit (e.g. of simulated data without
	     noise) will result in very small estimated uncertainties
	     on the parameters regardless of the size of the errors
	     specified for the fitted data points. With scaling turned
	     off, even a perfect fit will have large uncertainties on
	     the estimated parameters if the uncertainties on the data
	     are large. */
	  if (scale_covar)
	    {
	      // Scale covariance matrix with the square of sd.
	      double sqr_sd = pow(sd, 2);
	      int ioffset = 0;
	      for (int i=0; i<arg.n; i++)
		{
		  for (int j=ioffset; j<ioffset+arg.n; j++)
		    arg.fjac[j] *=sqr_sd;
		  ioffset += arg.ldfjac;
		}
	    }
	}
    }

  // Set minimized flag.
  minimized = 1;
}


// Return value for i'th free paramter.
double MinpackLM::get_value(unsigned int i)
{
  if (!fit_is_ok()) throw Uferr::LmNoFit();
  if ((int) i >= arg.n) throw Uferr::LmParamIndexInvalid(i);
  return (arg.x[i]);
}


// Return esd for i'th free paramter.
double MinpackLM::get_esd(unsigned int i)
{
  if (!fit_is_ok()) throw Uferr::LmNoFit();
  if (!is_converged()) throw Uferr::LmFitNotConverged();
  if ((int) i >= arg.n) throw Uferr::LmParamIndexInvalid(i);
  return ( sqrt(arg.fjac[i*arg.n+i]) );
}


// Return element from covariance matrix.
double MinpackLM::get_covar(unsigned int i, unsigned int j)
{
  if (!is_converged()) throw Uferr::LmFitNotConverged();
  if ((int) i>=arg.n) throw Uferr::LmParamIndexInvalid(i);
  if ((int) j>=arg.n) throw Uferr::LmParamIndexInvalid(j);
  return (arg.fjac[i*arg.n+j]);
}


/* lmstr_eval_ is the function to supply to lmstr (minpack). This
   procedure then calls the generic Fuda::eval_obs_calc function for
   each data point a calculation has to be performed for and
   eventually calculates derrivatives of the dependent parameters. */
void lmstr_eval_(int *m,
		 int *n,
		 double x[],
		 double fvec[],
		 double fjrow[],
		 int *iflag)
{
  // Get reference to MinpackLM structure and fuda structure.
  MinpackLM& lm = *MinpackLM::top_object_ref();
  Fuda& uf = *lm.get_fuda();

  // Branch out according to iflag.
  if (*iflag==1)
    {
      // We calculate function values for all data points and place them
      // in fvec[].

      // Loop over data.
      for (unsigned int d_index=0; d_index<uf.eval_get_ndata(); d_index++)
	{
          // Evaluate function value.
	  try {	    
	    fvec[d_index] = uf.eval_obs_calc(d_index, x, fjrow, 0);
	  }
	  catch (Uferr::FtypeCallError& e) {
	    std::cout << "lmstr_eval_ - FtypeCallError:\n"
		      << e.msg << "\n";
	    *iflag = -1;
	    break;
	  }	      
	  catch (...) {
	    std::cout << "lmstr_eval_ - unexpected exception caught "
		      << "from function value call to fuda method "
	              << "eval_obs_calc\n";
	    *iflag = -1;
	    break;
	  }	  
	}
    }
  else if (*iflag>1)
    {
      /* We calculate derivatives for the iflag-1'th data point (which is
      the datapoint with index iflag-2, when indexing from zero as in c++)
      and place them in fjrow[]. */

      // Calculate index of data point.
      int d_index = *iflag-2;

      /* Calculate derivatives of function with respect to dependent
         paramters */

      // Branch out according to numderiv flag.
      if (lm.get_numderiv())
	{
	  /* Numerical diferentiation by simple forward difference
             approximation adapted from minpack fdjac2 routine. */
	  double eps = sqrt(lm.get_numderiv_eps());
	  double fval0 = 0.0;
	  try {
	    fval0 = uf.eval_obs_calc(d_index, x, fjrow, 0);
	  }
	  catch (Uferr::FtypeCallError& e) {
	    std::cout << "lmstr_eval_ - FtypeCallError:\n"
		      << e.msg << "\n";
	    *iflag = -1;
	  }	      
	  catch (...) {
	    std::cout << "lmstr_eval_ - unexpected exception caught "
		      << "from function value call to fuda method "
	              << "eval_obs_calc\n";
	    *iflag = -1;
	  }	  
	  

	  // Loop over free parameters.
	  for (int j=0; j<*n; j++)
	    {
	      double temp = x[j];
	      double h = eps*fabs(temp);
	      if (h==0.0) h = eps;
	      x[j] = temp+h;
	      double fval1;
	      try {
		fval1 = uf.eval_obs_calc(d_index, x, fjrow, 0);
	      }	      
	      catch (Uferr::FtypeCallError& e) {
		std::cout << "lmstr_eval_ - FtypeCallError:\n"
			  << e.msg << "\n";
		*iflag = -1;
		break;
	      }	      
	      catch (...) {
		std::cout << "lmstr_eval_ - unexpected exception caught "
			  << "from function value call to fuda method "
			  << "eval_obs_calc\n";
		*iflag = -1;
		break;
	      }	  

	      x[j] = temp;
	      fjrow[j] = (fval1-fval0)/h;
	    }
	}
      else
	{
	  // Use analytical expressions for derivatives.  
	  try {
	    uf.eval_obs_calc(d_index, x, fjrow, 1);
	  }
	  catch (Uferr::FtypeCallError& e) {
	    std::cout << "lmstr_eval_ - FtypeCallError:\n"
		      << e.msg << "\n";
	    *iflag = -1;
	  }	      
	  catch (...) {
	    std::cout << "lmstr_eval_ - unexpected exception caught "
		      << "from derivative call to fuda method eval_obs_calc";
	    *iflag = -1;
	  }	  
	}
    }
  else if (*iflag==0)
    {
      /* Print status */

      // Calc. number of degrees of freedom.
      int nfree = *m-*n;
      // NB! nfree=m-n, but if nfree=0, we set nfree=1.
      if (nfree==0) nfree++;
      
      // Calc. enorm.
      double enorm = FUDA::calc_enorm(*m, fvec);

      // Calc. sd.
      double sd = FUDA::calc_sd(nfree, enorm);

      // Report sd.
      std::cout << "Iter: sd =  " << sd << "   "
		<< "enorm =  " << enorm <<"\n";
    }
}

