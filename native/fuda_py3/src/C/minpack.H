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
/* Declarations for FORTRAN MINPACK routines */ 

// Function type expected by lmstr_.
typedef void (*LMSTR_FUNC)(int *m,
			   int *n,
			   double x[],
			   double fvec[],
			   double fjrow[],
			   int *iflag);

// Function type expected by lmdif_.
typedef void (*LMDIF_FUNC)(int *m,
			   int *n,
			   double x[],
			   double fvec[],
			   int *iflag);

// Levenberg-Marquart minimization with minimized storage use.
extern "C" void lmstr_(LMSTR_FUNC fcn,
		       int *m,
		       int *n,
		       double x[],
		       double fvec[],
		       double fjac[],
		       int *ldfjac,
		       double *ftol,
		       double *xtol,
		       double *gtol,
		       int *maxfev,
		       double diag[],
		       int *mode,
		       double *factor,
		       int *nprint,
		       int *info,
		       int *nfev,
		       int *njev,
		       int ipvt[],
		       double qtf[],
		       double wa1[],
		       double wa2[],
		       double wa3[],
		       double wa4[]);

// Calculate the covariance matrix from fjac matrix returned by lmstr_.
extern "C" void covar_(int *n,
		       double r[],
		       int *ldr,
		       int ipvt[],
		       double *tol,
		       double wa[]);

// Calculate euclidian norm of a vector.
extern "C" double enorm_(int *n,
			 double x[]);

// Return machine specific parameters.
extern "C" double dpmpar_(const int *i);











