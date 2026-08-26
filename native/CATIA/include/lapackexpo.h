#ifndef LAPACKEXP_H
#define LAPACKEXP_H

/*
*  Turn on HAVE_LAPACK_CONFIG_H to redefine C-LAPACK datatypes
*/
#ifdef HAVE_LAPACK_CONFIG_H
#include "lapacke_config.h"
#endif

#include "lapacke_mangling.h"

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
#ifndef lapack_int
#define lapack_int     int
#endif

#ifndef lapack_logical
#define lapack_logical lapack_int
#endif

/* f2c, hence clapack and MacOS Accelerate, returns double instead of float
 * for sdot, slange, clange, etc. */
#if defined(LAPACK_F2C)
    typedef double lapack_float_return;
#else
    typedef float lapack_float_return;
#endif

/* Complex types are structures equivalent to the
* Fortran complex types COMPLEX(4) and COMPLEX(8).
*
* One can also redefine the types with his own types
* for example by including in the code definitions like
*
* #define lapack_complex_float std::complex<float>
* #define lapack_complex_double std::complex<double>
*
* or define these types in the command line:
*
* -Dlapack_complex_float="std::complex<float>"
* -Dlapack_complex_double="std::complex<double>"
*/

#ifndef LAPACK_COMPLEX_CUSTOM

/* Complex type (single precision) */
#ifndef lapack_complex_float
#include <complex.h>
#define lapack_complex_float    float _Complex
#endif

#ifndef lapack_complex_float_real
#define lapack_complex_float_real(z)       (creal(z))
#endif

#ifndef lapack_complex_float_imag
#define lapack_complex_float_imag(z)       (cimag(z))
#endif

/* Complex type (double precision) */
#ifndef lapack_complex_double
#include <complex.h>
#define lapack_complex_double   double _Complex
#endif

#ifndef lapack_complex_double_real
#define lapack_complex_double_real(z)      (creal(z))
#endif

#ifndef lapack_complex_double_imag
#define lapack_complex_double_imag(z)       (cimag(z))
#endif

#endif /* LAPACK_COMPLEX_CUSTOM */

/* Callback logical functions of one, two, or three arguments are used
*  to select eigenvalues to sort to the top left of the Schur form.
*  The value is selected if function returns TRUE (non-zero). */

 
typedef lapack_logical (*LAPACK_S_SELECT2) ( const float*, const float* );
typedef lapack_logical (*LAPACK_S_SELECT3)
    ( const float*, const float*, const float* );
typedef lapack_logical (*LAPACK_D_SELECT2) ( const double*, const double* );
typedef lapack_logical (*LAPACK_D_SELECT3)
    ( const double*, const double*, const double* );

typedef lapack_logical (*LAPACK_C_SELECT1) ( const lapack_complex_float* );
typedef lapack_logical (*LAPACK_C_SELECT2)
    ( const lapack_complex_float*, const lapack_complex_float* );
typedef lapack_logical (*LAPACK_Z_SELECT1) ( const lapack_complex_double* );
typedef lapack_logical (*LAPACK_Z_SELECT2)
    ( const lapack_complex_double*, const lapack_complex_double* );

#define LAPACK_lsame LAPACK_GLOBAL(lsame,LSAME)
lapack_logical LAPACK_lsame( char* ca,  char* cb,
                              lapack_int lca, lapack_int lcb );


/*----------------------------------------------------------------------------*/
/* This is in alphabetical order (ignoring leading precision). */


#define LAPACK_dgpadm LAPACK_GLOBAL(dgpadm,DGPADM)
void LAPACK_dgpadm(
		   lapack_int *ideg,
		   lapack_int *m,
		   double *t,
		   double *H,
		   lapack_int const *ldh,
		   double *wsp,
		   int *lwsp,
		   double *ipiv,
		   int *iexph,
		   int *ns,
		   lapack_int *iflag
		   );
  
#define LAPACK_zgpadm LAPACK_GLOBAL(zgpadm,ZGPADM)
void LAPACK_zgpadm(
		   lapack_int *ideg,
		   lapack_int *m,
		   double *t,
		   lapack_complex_double *H,
		   lapack_int const *ldh,
		   lapack_complex_double *wsp,
		   int *lwsp,
		   lapack_complex_double *ipiv,
		   int *iexph,
		   int *ns,
		   lapack_int *iflag
		   );




#if defined(BIND_FORTRAN_LOWERCASE_UNDERSCORE) || defined(BIND_FORTRAN_LOWERCASE)
// Allow manual override of the defaults, e.g. if you want to use a fortran
// lib compiled with gcc from MSVC
#else

// First we need to know what the conventions for linking
// C with Fortran is on this platform/toolset
#if defined(__GNUC__) || defined(__ICC) || defined(__sgi) || defined(__COMO__) || defined(__KCC)
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#elif defined(__IBMCPP__) || defined(_MSC_VER)
#define BIND_FORTRAN_LOWERCASE
#else
#error do not know how to link with fortran for the given platform
#endif

#endif


//from boost numerical bindings.

  
// Next we define macro's to convert our symbols to 
// the current convention
#if defined(BIND_FORTRAN_LOWERCASE_UNDERSCORE)
#define FORTRAN_ID( id ) id##_
#elif defined(BIND_FORTRAN_LOWERCASE)
#define FORTRAN_ID( id ) id
#else
#error do not know how to bind to fortran calling convention
#endif

  

#define BLAS_DGEMM FORTRAN_ID( dgemm )
#define BLAS_ZGEMM FORTRAN_ID( zgemm )
#define BLAS_DGEMV FORTRAN_ID( dgemv )
#define BLAS_DGEMV2 FORTRAN_ID( dgemv2 )
#define BLAS_ZGEMV FORTRAN_ID( zgemv )

  void   BLAS_DGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double     *alpha, const double     *a, const int *lda, const double     *b, const int *ldb, const double     *beta, double     *c, const int *ldc);

  void   BLAS_DGEMV(const char *trans, const int *m, const int *n, const double   *alpha, const double   *a, const int *lda, const double   *x, const int *incx, const double   *beta, double   *y, const int *incy) ;

  void   BLAS_DGEMV2(const int *n, const double   *a, double   *y) ; //propagate vector y by matrix a.

  void   BLAS_ZGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k, lapack_complex_double *alpha, lapack_complex_double *a, const int *lda, lapack_complex_double *b, const int *ldb, lapack_complex_double *beta, lapack_complex_double *c, const int *ldc);

void   BLAS_ZGEMV(const char *trans, const int *m, const int *n,const lapack_complex_double *alpha, const lapack_complex_double *a, const int *lda, const lapack_complex_double *x, const int *incx, const lapack_complex_double *beta, lapack_complex_double *y, const int *incy) ;

  
  
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LAPACKEXPO_H */
