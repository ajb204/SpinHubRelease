*> \brief \b DGEMV2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMV(N,A,X,Y)
*
*       .. Scalar Arguments ..
*       INTEGER N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(N,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMV2  performs one of the matrix-vector operations
*>
*>     y := alpha*A**T*x,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> n by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of rows/columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N )
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMV2(N,A,X)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*  -- Simplified version, alpha=1, beta=0, matrices aligned and A is square. For use in CATIA
*     
*     .. Scalar Arguments ..
      INTEGER N
*     .. Array Arguments ..
      DOUBLE PRECISION A(N,*),X(*)
*  =====================================================================
*     .. Parameters ..
      DOUBLE PRECISION Y(N)
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,J,JY
*     First form  y := beta*y.
*                  DO 10 I = 1,N
*                      Y(I) = ZERO
*   10             CONTINUE
*      Y=0
*     Form  y := alpha*A**T*x + y.
*      JY = 1
      DO 100 J = 1,N
         TEMP = 0
         DO 90 I = 1,N
            TEMP = TEMP + A(I,J)*X(I)
   90    CONTINUE
         Y(J) = TEMP
  100 CONTINUE
* copy up array
      DO 200 I = 1,N
         X(I) = Y(I)
  200 CONTINUE
      
*     R      ETURN
*
*     End of DGEMV2 .
*
      END
