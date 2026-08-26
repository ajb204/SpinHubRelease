/***/
/* vutil: simple vector utilities.
/***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prec.h"
#include "rand.h"
#include "vutil.h"

#define FABS(X) (((X) < 0.0) ? -(X) : (X))

static int rSortProc();
static int rASortProc();

static int rSortProc( rI, rJ )
 
   float *rI, *rJ;
{
    if (*rI > *rJ) return( 1 );
    if (*rI < *rJ) return( -1 );

    return (0);
}

static int rASortProc( rI, rJ )

   float *rI, *rJ;
{
    float aI, aJ;

    aI = FABS( *rI );
    aJ = FABS( *rJ );

    if (aI > aJ) return( 1 );
    if (aI < aJ) return( -1 );

    return (0);
}

int vSort64( src, size )

   float   *src;
   NMR_INT size;
{
    (void) qsort( (void *)src, (size_t)size, sizeof(float), rSortProc );
    return( 0 );
}

int vASort64( src, size )

   float   *src;
   NMR_INT size;
{
    (void) qsort( (void *)src, (size_t)size, sizeof(float), rASortProc );
    return( 0 );
}

/***/
/* iSet: set elements of vec to val.
/***/

int iSet64( vec, n, val )

   int     *vec, val;
   NMR_INT n;
{
    NMR_INT i;

    if (!vec || n < 1) return( 0 );
 
    for( i = 0; i < n; i++ ) *vec++ = val;

    return( 0 );
}

int vvProd64( dest, srcX, srcY, nX, nY )

   float   *dest, *srcX, *srcY;
   NMR_INT nX, nY;
{
    float   *xPtr, y;
    NMR_INT ix, iy;

    if (!dest || !srcX || !srcY) return( 0 );
    if (nX < 1 || nY < 1) return( 0 );

    for( iy = 0; iy < nY; iy++ )
       {
        y    = *srcY++;
        xPtr = srcX;

        for( ix = 0; ix < nX; ix++ ) *dest++ += (*xPtr++)*y;
       }

    return( 0 );
}

int vvvProd64( dest, srcX, srcY, srcZ, nX, nY, nZ )

   float   *dest, *srcX, *srcY, *srcZ;
   NMR_INT nX, nY, nZ;
{
    NMR_INT ix, iy, iz;

    if (!dest || !srcX || !srcY || !srcZ) return( 0 );
    if (nX < 1 || nY < 1 || nZ < 1) return( 0 );

    for( iz = 0; iz < nZ; iz++ )
       {
        for( iy = 0; iy < nY; iy++ )
           {
            for( ix = 0; ix < nX; ix++ ) *dest++ += srcX[ix]*srcY[iy]*srcZ[iz];
           }
       }

    return( 0 );
}

int vvvvProd64( dest, srcX, srcY, srcZ, srcA, nX, nY, nZ, nA )

   float   *dest, *srcX, *srcY, *srcZ, *srcA;
   NMR_INT nX, nY, nZ, nA;
{
    NMR_INT ix, iy, iz, ia;

    if (!dest || !srcX || !srcY || !srcZ || !srcA) return( 0 );
    if (nX < 1 || nY < 1 || nZ < 1 || nA < 1) return( 0 );

    for( ia = 0; ia < nA; ia++ )
       {
        for( iz = 0; iz < nZ; iz++ )
           {
            for( iy = 0; iy < nY; iy++ )
               {
                for( ix = 0; ix < nX; ix++ ) *dest++ += srcX[ix]*srcY[iy]*srcZ[iz]*srcA[ia];
               }
           }
       }

    return( 0 );
}

int vvCopy64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    while( size-- ) *dest++ = *src++;

    return( 0 );
}

int vvCopyOff64( dest, src, size, destOff, srcOff )

   float   *dest, *src;
   NMR_INT size, destOff, srcOff;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    dest += destOff;
    src  += srcOff;

    while( size-- ) *dest++ = *src++;

    return( 0 );
}

int vvSetSign64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
   NMR_INT i;
   float   v;

   if (!dest || !src || size < 1) return( 0 );

   for( i = 0; i < size; i++ )
      {
       v     = *dest < 0.0 ? -(*dest) : *dest;
       *dest =  *src < 0.0 ? -v : v;

       src++; 
       dest++;
      }

   return( 0 );
}

/***/
/* vSign: replace values in vec with -1 0 or 1 according to their signs.
/***/

int vSign64( vec, length )

   float   *vec;
   NMR_INT length;
{
    if (!vec || length < 1) return( 0 );

    while( length-- )
       {
        if (*vec < 0.0)
           *vec = -1.0;
        else if (*vec > 0.0)
           *vec = 1.0;
        else
           *vec = 0.0;

        vec++;
       }

   return( 0 );
}

/***/
/* vFixNaN: replace NaN values in vec with 0.
/***/

int vFixNaN64( vec, length )

   float   *vec;
   NMR_INT length;
{
    int n, status;

    if (!vec || length < 1) return( 0 );

    n = 0;

    while( length-- )
       {
#ifdef NMR64
        status = fpclassify( *vec );

        if (status != FP_NORMAL && status != FP_ZERO)
           {
            *vec = 0.0;
            n++;
           }
#else
        if (isnan( *vec ) || isinf( *vec ) || (*vec < 0.0 && *vec > -1.0e-24) || (*vec > 0.0 && *vec < 1.0e-24))
           {
            *vec = 0.0;
            n++;
           }
#endif
        vec++;
       }

   return( n );
}


int vNeg64( v, size )

   float   *v;
   NMR_INT size;
{
   NMR_INT i;

   if (!v || size < 1) return( 0 );

   for( i = 0; i < size; i++ )
      {
       *v = -(*v);
       v++;
      }

   return( 0 );
}

int vAbs64( v, size )

   float   *v;
   NMR_INT size;
{
   NMR_INT i;

   if (!v || size < 1) return( 0 );

   for( i = 0; i < size; i++ )
      {
       *v = *v < 0.0 ? -(*v) : *v;
       v++;
      }

   return( 0 );
}

/***/
/* vAlt: sign-alternate data in vec (+ -).
/***/

int vAlt64( vec, length )

   float   *vec;
   NMR_INT length;
{
    if (!vec || length < 1) return( 0 );

    length /= 2;
    vec++;

    while( length-- )
       {
        *vec = -(*vec);
        vec += 2;
       }

    return( 0 );
}

int vRev64( vec, length )

   float   *vec;
   NMR_INT length;
{
    float *rPtr, rtemp;

    if (!vec || length < 1) return( 0 );

    rPtr    = vec + length - 1;
    length /= 2;

    while( length-- )
       {
        rtemp   = *vec;
        *vec++  = *rPtr;
        *rPtr-- = rtemp;
       }

    return( 0 );
}

int vvMult64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ *= *src++;

    return( 0 );
}

int vvDiv64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ )
       {
        if (*src == 0.0)
           *dest++ = 0.0;
        else
           *dest++ /= *src;
       }

    return( 0 );
}

int vvSwap64( v1, v2, size )

   float   *v1, *v2;
   NMR_INT size;
{
    float rtemp;

    if (!v1 || !v2 || size < 1) return( 0 );

    while( size-- )
       {
        rtemp = *v1;
        *v1++ = *v2;
        *v2++ = rtemp;
       }

    return( 0 );
}


int vvAdd64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ += *src++;

    return( 0 );
}

int vvSub64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ -= *src++;

    return( 0 );
}

int vvMask64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ )
       {
        if (*src++ == 0.0) *dest = 0.0;
        dest++;
       }

    return( 0 );
}

int vvIMask64( dest, src, size )

   float   *dest, *src;
   NMR_INT size;
{  
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) 
       { 
        if (*src++ != 0.0) *dest = 0.0;
        dest++; 
       } 
 
    return( 0 );
}

int vsSet64( dest, size, val )

   float   *dest, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ = val;

    return( 0 );
}

int vsAdd64( dest, size, val )

   float   *dest, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ += val;

    return( 0 );
}

int vsAddGRand64( dest, size, val )

   float   *dest, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ += val*gaussRand();

    return( 0 );
}

int vsMult64( dest, size, val )

   float   *dest, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ *= val;

    return( 0 );
}

int vsClipBelow64( dest, size, thresh )

   float   *dest, thresh;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    if (thresh < 0.0) thresh = -thresh;

    for( i = 0; i < size; i++ )
       {
        if (*dest > -thresh && *dest < thresh) *dest = 0.0;;
        dest++;
       }

    return( 0 );
}

int vsClipAbove64( dest, size, thresh )

   float   *dest, thresh;
   NMR_INT size; 
{
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    if (thresh < 0.0) thresh = -thresh;

    for( i = 0; i < size; i++ ) 
       { 
        if (*dest < -thresh)
           *dest = -thresh;
        else if (*dest > thresh)
           *dest = thresh;

        dest++;
       }

    return( 0 );
}

int vsMultThresh64( vec, length, val, thresh )

   float   *vec, val, thresh;
   NMR_INT length;
{
    float   pthresh, mthresh;
    NMR_INT i;

    pthresh = FABS( thresh );
    mthresh = -pthresh;

    if (thresh < 0.0)
       {
        for( i  = 0; i < length; i++ )
           {
            if (*vec > pthresh)
               {
                *vec -= pthresh + (*vec - pthresh)*val;
               }
            else if (*vec < mthresh)
               {
                *vec -= mthresh + (*vec - mthresh)*val;
               }
            else
               {
                *vec = 0.0;
               }

            vec++;
           }
       }
    else
       {
        for( i  = 0; i < length; i++ )
           {
            if (*vec > pthresh)
               {
                *vec = pthresh + (*vec - pthresh)*val;
               }
            else if (*vec < mthresh)
               {
                *vec = mthresh + (*vec - mthresh)*val;
               }

            vec++;
           }
       }

    return( 0 );
}

int vqMult64( dest, size, valR, valI )

   float   *dest, valR, valI;
   NMR_INT size;
{
    float   *bR, *bI, tmpR, tmpI;
    NMR_INT i;

    if (!dest || size < 1) return( 0 );

    bR = dest; bI = dest + size;

    for( i = 0; i < size; i++ )
       {
        tmpR = valR*(*bR) - valI*(*bI);
        tmpI = valI*(*bR) + valR*(*bI);

        *bR++ = tmpR;
        *bI++ = tmpI;
       }

    return( 0 );
}

int vvsMult64( dest, src, size, val )

   float  *dest, *src, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ = val*(*src++);

    return( 0 );
}

int vsvAdd64( dest, src, size, val )

   float   *dest, *src, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ += val*(*src++);

    return( 0 );
}

int vsvMult64( dest, src, size, val )

   float   *dest, *src, val;
   NMR_INT size;
{
    NMR_INT i;

    if (!dest || !src || size < 1) return( 0 );

    for( i = 0; i < size; i++ ) *dest++ *= val*(*src++);

    return( 0 );
}

float vSum64( src, size )

   float   *src;
   NMR_INT size;
{
    float   sum;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    sum = 0.0;

    for( i = 0; i < size; i++ ) sum += *src++;

    return( sum );
}

float vSumPos64( src, size )

   float   *src;
   NMR_INT size;
{
    float   sum;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    sum = 0.0;

    for( i = 0; i < size; i++ )
       {
        if (*src > 0.0) sum += *src;
        src++;
       }

    return( sum );
}

float vSumNeg64( src, size )

   float   *src;
   NMR_INT size;
{
    float   sum;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    sum = 0.0;

    for( i = 0; i < size; i++ ) 
       { 
        if (*src < 0.0) sum += *src;
        src++;
       }

    return( sum );
}

float vMin64( src, size )

   float   *src;
   NMR_INT size;
{
    float   v;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    v = *src;

    for( i = 0; i < size; i++ )
       {
        if (*src < v) v = *src;
        src++;
       }

    return( v );
}

float vMax64( src, size )

   float   *src;
   NMR_INT size;
{
    float   v;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    v = *src;

    for( i = 0; i < size; i++ )
       {
        if (*src > v) v = *src;
        src++;
       }

    return( v );
}

float vMinAbs64( src, size )

   float   *src;
   NMR_INT size;
{
    float   vA, qA;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    vA = *src < 0.0 ? -(*src) : *src;

    for( i = 0; i < size; i++ )
       {
        qA = *src < 0.0 ? -(*src) : *src;
        if (qA < vA) vA = qA;
        src++;
       }

    return( vA );
}

float vMaxAbs64( src, size )

   float   *src;
   NMR_INT size;
{
    float   vA, qA;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    vA = *src < 0.0 ? -(*src) : *src;

    for( i = 0; i < size; i++ )
       {
        qA = *src < 0.0 ? -(*src) : *src;
        if (qA > vA) vA = qA;
        src++;
       }

    return( vA );
}

float vAMax64( src, size )

   float   *src;
   NMR_INT size;
{
    float   vA, qA;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    vA = *src;

    for( i = 0; i < size; i++ )
       {
        if (FABS(*src) > FABS(vA)) vA = *src;
        src++;
       }

    return( vA );
}

float vAMin64( src, size )

   float   *src;
   NMR_INT size;
{
    float   vA, qA;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    vA = *src;

    for( i = 0; i < size; i++ )
       {
        if (FABS(*src) < FABS(vA)) vA = *src;
        src++;
       }

    return( vA );
}

float vSumSq64( src, size )

   float   *src;
   NMR_INT size;
{
    float   sum;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    sum = 0.0;

    for( i = 0; i < size; i++ )
       {
        sum += (*src)*(*src);
        src++;
       }

    return( sum );
}

float vRMS64( src, size )

   float   *src;
   NMR_INT size;
{
    float   sum;
    NMR_INT i;

    if (!src || size < 1) return( (float)0.0 );

    sum = 0.0;

    for( i = 0; i < size; i++ )
       {
        sum += (*src)*(*src);
        src++;
       }

    if (size) sum = sqrt( (double)(sum/size) );

    return( sum );
}

float vVariance64( src, size )

   float   *src;
   NMR_INT size;
{
    float   *rPtr, avg, var, d;
    NMR_INT i;

    if (!src || size < 2) return( (float)0.0 );

    avg = 0.0;
    var = 0.0;

    rPtr = src;

    for( i = 0; i < size; i++ ) avg += *rPtr++;

    avg  /= (float)size;
    rPtr  = src;

    for( i = 0; i < size; i++ )
       {
        d = *rPtr++ - avg; 
        var += d*d;
       }

    var /= (float)(size - 1);
 
    return( var );
}

float vStdDev64( src, size )

   float   *src;
   NMR_INT size;
{  
    float dev;

    if (!src || size < 2) return( (float)0.0 );

    dev = vVariance64( src, size );
    dev = sqrt( (double)dev );

    return( dev );
}

float vMedian64( src, size )

   float   *src;
   NMR_INT size;
{
    float *work, val;

    if (!src || size < 1) return( (float)0.0 );

    val = 0.0;

    return( val );
}

float vNoiseEst64( src, size )

   float   *src;
   NMR_INT size;
{
    float *work, val;

    if (!src || size < 1) return( (float)0.0 );

    val = 0.0;

    return( val );
}

float vNorm64( vec, length )

   float   *vec;
   NMR_INT length;
{
    float val;

    if (!vec || length < 1) return( (float)0.0 );

    val = 0.0;

    while( length-- )
       {
        val += (*vec)*(*vec);
        vec++;
       }

    val = sqrt( (double)val );

    return( val );
}

int applyNorm64( vec, ix1, ixn, kx1, kxn, normFunc, normMode, cC )

   NMR_INT ix1, ixn, kx1, kxn;
   int     normFunc, normMode;
   float   *vec, cC;
{
    NMR_INT kSize, iSize;
    float  val;

    if (!vec || ix1 > ixn || kx1 > kxn) return( 1 );

    val = 0.0;

    /* (void) fprintf( stderr,
                       "Applying Norm Func %d I: %d %d K: %d %d\n",
                       normFunc, (int)ix1, (int)ixn, (int)kx1, (int)kxn ); */

    iSize = 1 + ixn - ix1;
    kSize = 1 + kxn - kx1;

    switch( normFunc )
       {
        case NORM_FUNC_SUM:
           val = vSum64( vec + kx1, (NMR_INT)(1 + kxn - kx1) );
           break;
        case NORM_FUNC_AVG:
           val  = vSum64( vec + kx1, kSize )/(float)kSize;
           break;
        case NORM_FUNC_VARIANCE:
           val = vVariance64( vec + kx1, kSize );
           break;
        case NORM_FUNC_STDDEV:
           val = vStdDev64( vec + kx1, kSize );
           break;
        case NORM_FUNC_RMS:
           val = vRMS64( vec + kx1, kSize );
           break;
        case NORM_FUNC_MIN:
           val = vMin64( vec + kx1, kSize );
           break;
        case NORM_FUNC_MAX:
           val = vMax64( vec + kx1, kSize );
           break;
        case NORM_FUNC_MINABS:
           val = vMinAbs64( vec + kx1, kSize );
           break;
        case NORM_FUNC_MAXABS:
           val = vMaxAbs64( vec + kx1, kSize );
           break;
       case NORM_FUNC_AMIN:
           val = vAMin64( vec + kx1, kSize );
           break;
        case NORM_FUNC_AMAX:
           val = vAMax64( vec + kx1, kSize );
        case NORM_FUNC_RANGE:
           val = vMax64( vec + kx1, kSize ) - vMin64( vec + kx1, kSize );
           break;
        case NORM_FUNC_MEDIAN:
           val = vMedian64( vec + kx1, kSize );
           break;
        case NORM_FUNC_NOISE:
           val = vNoiseEst64( vec + kx1, kSize );
           break;
        case NORM_FUNC_NONE:
        default:
           return( 0 );
           break;
       };

    val *= cC;

    /* (void) fprintf( stderr, "Applying Norm Mode %d Val: %f\n", normMode, val ); */

    switch( normMode )
       {
        case NORM_MODE_ADD:
           (void) vsAdd64( vec + ix1, iSize, val );
           break;
        case NORM_MODE_SUBTRACT:
           (void) vsAdd64( vec + ix1, iSize, -val );
           break;
        case NORM_MODE_MULT:
           (void) vsMult64( vec + ix1, iSize, val );
           break;
        case NORM_MODE_DIVIDE:
           val = val == 0.0 ? 0.0 : 1.0/val;
           (void) vsMult64( vec + ix1, iSize, val );
           break;
        case NORM_MODE_REPLACE:
           (void) vsSet64( vec + ix1, iSize, val );
           break;
        default:
           break;
       };
 
    return( 0 );
}

/***/
/* getVec: get a dispersed vector from src.
/***/

int getVec64( src, dest, eCount, eJump )

   float    *src, *dest;
   NMR_INT  eCount, eJump;
{
    if (!src || !dest || eCount < 1) return( 0 );

    while( eCount-- )
       {
        *dest++  = *src;
        src     += eJump;
       }

    return( 0 );
}

/***/
/* putVec: put a vector back in dispersed form.
/***/

int putVec64( dest, src, eCount, eJump )

   float   *dest, *src;
   NMR_INT eCount, eJump;
{
    if (!src || !dest || eCount < 1) return( 0 );

    while( eCount-- )
       {
        *dest  = *src++;
        dest  += eJump;
       }

    return( 0 );
}

/***/
/* vvGetScale: find s to minimize a - s*b.
/***/

float vvGetScale64( a, b, n )

   float   *a, *b;
   NMR_INT n;
{
    float   top, bot, s;
    NMR_INT i;

    if (!a || !b || n < 1) return( (float)0.0 );

    top = 0.0;
    bot = 0.0;

    for( i = 0; i < n; i++ )
       {
        top += a[i]*b[i];
        bot += b[i]*b[i];
       }

    if (bot == 0.0) return( (float)0.0 );

    s = top/bot;

    return( s );
}

float vvScale64( a, b, n )

   float   *a, *b;
   NMR_INT n;
{
    float s;

    if (!a || !b || n < 1) return( (float)0.0 );

    s = vvGetScale64( a, b, n );

    (void) vsMult64( b, n, s );

    return( s );
}

float vvDot64( vec1, vec2, length )

   float   *vec1, *vec2;
   NMR_INT length;
{
    float val;

    if (!vec1 || !vec2 || length < 1) return( (float)0.0 );

    val = 0.0;

    while( length-- )
       {
        val += (*vec1++)*(*vec2++);
       }

    return( val );
}

/***/
/* vvGetScaleRI: find s to minimize a - s*b, complex version; returns phase.
/***/

float vvGetScaleRI64( a, b, n, sR, sI )

   float   *a, *b, *sR, *sI;
   NMR_INT n;
{
    float   *aR, *aI, *bR, *bI, topR, topI, bot, abR, bbR, abI, bbI, s;
    NMR_INT i;

    if (!a || !b || n < 1) return( (float)0.0 );

    abR = 0.0; abI = 0.0;
    bbR = 0.0; bbI = 0.0;

    aR = a; aI = a + n;
    bR = b; bI = b + n;

    *sR = 0.0;
    *sI = 0.0;

    s = 0.0;

    for( i = 0; i < n; i++ )
       {
        abR += (*aR)*(*bR) - (*aI)*(*bI);
        abI += (*aI)*(*bR) + (*aR)*(*bI);

        bbR += (*bR)*(*bR) - (*bI)*(*bI);
        bbI += 2.0*(*bR)*(*bI);

        aR++; aI++;
        bR++; bI++;
       }

    bot = bbR*bbR + bbI*bbI;

    if (bot == 0.0) return( s );

    *sR = (abR*bbR + abI*bbI)/bot;
    *sI = (abI*bbR - abR*bbI)/bot;

    s   = 180.0*atan2( -(double)(*sI), (double)(*sR) )/PI;

    return( s );
}

float vvRMS64( a, b, n )

   float   *a, *b;
   NMR_INT n;
{
    float   d, s;
    NMR_INT i;

    s = 0.0;

    if (!a || !b || n < 1) return( s );

    for( i = 0; i < n; i++ )
       {
        d  = a[i] - b[i];
        s += d*d;
       }

    s = sqrt( (double)(s/n) );

    return( s );
}

float vvPearsonR64( xList, yList, n )

   float   *xList, *yList;
   NMR_INT n;
{
    float   x, y, xd, yd, xMax, yMax, s, top, bot, sumX, sumX2, sumY, sumY2, sumXY, r;
    NMR_INT i;

    r = 0.0;

    if (!xList || !yList || n < 1) return( r );

    sumX  = 0.0;
    sumX2 = 0.0;
    sumY  = 0.0;
    sumY2 = 0.0;
    sumXY = 0.0;

    xMax = vMaxAbs64( xList, n );
    yMax = vMaxAbs64( yList, n );

    s = xMax > yMax ? xMax : yMax;
    s = s == 0.0 ? 1.0 : 1.0/s;
 
    for( i = 0; i < n; i++ )
       {
        x = s*(*xList++); y = s*(*yList++);

        sumX  += x;   sumY  += y;
        sumX2 += x*x; sumY2 += y*y;
        sumXY += x*y;
       }


    xd = sumX2 - sumX*sumX/(double)n;
    yd = sumY2 - sumY*sumY/(double)n;

    top = sumXY - sumX*sumY/(double)n;
    bot = sqrt( fabs(xd*yd) );

    if (bot == 0.0) return( r );

    r = top/bot;

    return( r );
}

int coAdd64( vec, length, aList, aCount )

   float   *vec, *aList;
   NMR_INT length, aCount;
{
    float   *aPtr, *vPtr, sum;
    NMR_INT i, j;

    if (!vec || length < 1 || aCount < 1) return( 0 );

    vPtr = vec;

    for( i = 0; i < length; i += aCount )
       {
        aPtr  = aList;
        sum   = 0.0;

        for( j = 0; j < aCount; j++ ) sum += (*vPtr++)*(*aPtr++);

        *vec++ = sum;
       }

    return( 0 );
}

int cSet64( cvec, count, cval )

   char    *cvec, cval;
   NMR_INT count;
{
    if (!cvec || count < 1) return( 0 );

    while( count-- ) *cvec++ = cval;

    return( 0 );
}

int rSet64( vec, length, val )

   float   *vec, val;
   NMR_INT length;
{
    if (!vec || length < 1) return( 0 );

    while( length-- ) *vec++ = val;

    return( 0 );
}

int rMove64( src, dest, length )

   float   *src, *dest;
   NMR_INT length;
{
    if (!src || !dest || length < 1) return( 0 );

    while( length-- ) *dest++ = *src++;

    return( 0 );
}

int vByteSwap64( vec, length )

   float   *vec;
   NMR_INT length;
{
    union fsw { float f; char s[4]; } in, out;

    if (!vec || length < 1) return( 0 );

    while( length-- )
       {
        in.f = *vec;

        out.s[0] = in.s[3];
        out.s[1] = in.s[2];
        out.s[2] = in.s[1];
        out.s[3] = in.s[0];

        *vec++ = out.f;
       }

    return( 0 );
}

int iByteSwap64( vec, length )

   int     *vec;
   NMR_INT length;
{
    union b2i { int i; char s[4]; } in, out;

    if (!vec || length < 1) return( 0 );

    while( length-- )
       {
        in.i = *vec;

        out.s[0] = in.s[3];
        out.s[1] = in.s[2];
        out.s[2] = in.s[1];
        out.s[3] = in.s[0];

        *vec++ = out.i;
       }

    return( 0 );
}

int vInterpol64( dest, src, inSize, outSize )

   float   *src, *dest;
   NMR_INT inSize, outSize;
{
    NMR_INT ix, xLoc0, xLoc1;
    float   dx, x0;

    if (!src || !dest || inSize < 1 || outSize < 1) return( 0 );

    for( ix = 0; ix < outSize; ix++ )
       {
        x0 = outSize > 1 ? (float)ix*(inSize - 1)/(outSize - 1) : 0;

        if (x0 < 0.0) x0 = 0.0;
        if (x0 > inSize - 1) x0 = inSize - 1;

        xLoc0 = (NMR_INT)x0;

        if (xLoc0 < 0) xLoc0 = 0;
        if (xLoc0 > inSize - 1) xLoc0 = inSize - 1;

        dx = x0 - xLoc0;

        if (dx < 0.0) dx = 0.0;
        if (dx > 1.0) dx = 1.0;

        xLoc1 = xLoc0 < inSize - 1 ? xLoc0 + 1 : xLoc0;

        dest[ix] = (1.0 - dx)*src[xLoc0] + dx*src[xLoc1];
       }

    return( 0 );
}

int bilinInterp2D64( inMatrix,  xSize1, ySize1, outMatrix, xSize2, ySize2 )

   float    *inMatrix, *outMatrix; 
   NMR_INT  xSize1, xSize2, ySize1, ySize2;
{
   NMR_INT loc00, loc10, loc01, loc11;

   float    *dest, x0, y0, dx, dy;
   NMR_INT  xLoc, yLoc, ix, iy;

   if (!inMatrix || !outMatrix || xSize1 < 1 || ySize1 < 1 || xSize2 < 1 || ySize2 < 1) return( 0 );

   dest = outMatrix;
 
   for( iy = 0; iy < ySize2; iy++ )
      {
       y0   = ySize2 > 1 ? (float)iy*(ySize1 - 1)/(ySize2 - 1) : 0;

       if (y0 < 0.0) y0 = 0.0;
       if (y0 > ySize1 - 1) y0 = ySize1 - 1;

       yLoc = (int)y0;

       if (yLoc < 0) yLoc = 0;
       if (yLoc > ySize1 - 1) yLoc = ySize1 - 1;

       dy = y0 - yLoc;

       if (dy < 0.0) dy = 0.0;
       if (dy > 1.0) dy = 1.0;

       for( ix = 0; ix < xSize2; ix++ )
          {
           x0   = xSize2 > 1 ? (float)ix*(xSize1 - 1)/(xSize2 - 1) : 0;

           if (x0 < 0.0) x0 = 0.0;
           if (x0 > xSize1 - 1) x0 = xSize1 - 1;

           xLoc = (int)x0;

           if (xLoc < 0) xLoc = 0;
           if (xLoc > xSize1 - 1) xLoc = xSize1 - 1;

           dx = x0 - xLoc;

           if (dx < 0.0) dx = 0.0;
           if (dx > 1.0) dx = 1.0;

           loc00 = (NMR_INT)xLoc + xSize1*yLoc;
           loc10 = xLoc < xSize1 - 1 ? loc00 + 1      : loc00;
           loc01 = yLoc < ySize1 - 1 ? loc00 + xSize1 : loc00;
           loc11 = xLoc < xSize1 - 1 ? loc01 + 1      : loc01;

           *dest++ = inMatrix[loc00]*(1.0 - dx)*(1.0 - dy)
                     + inMatrix[loc10]*dx*(1.0 - dy)
                       + inMatrix[loc01]*(1.0 - dx)*dy 
                         + inMatrix[loc11]*dx*dy;
          }
      }

   return( 0 );
}

int vAbsCor64( dest, src, bg, c, thresh, length )

   float   *dest, *src, *bg, c, thresh;
   NMR_INT length;
{
    float   f, d;
    NMR_INT i;

    if (!src || !dest || length < 1) return( 0 );

    for( i = 0; i < length; i++ )
       {
        f = bg[i]*c;

        if (bg[i] < thresh && f < src[i])
           {
            dest[i] = 0.0;
           }
        else
           {
            d = f == 0.0 ? 1.0 : src[i]/f;
            d = d <  0.0 ? -d  : d;
            dest[i] = d == 0.0 ? 0.0 : -log10( (double)d );
           }
       }

    return( 0 );
}

int phaseRI64( rdata, idata, ix1, ix3, length, p0, p1 )

   float     *rdata, *idata, p0, p1;
   NMR_INT   ix1, ix3, length;
{
    float   phi, tR, tI, tCOS, tSIN;
    NMR_INT i;

    if (!rdata || !idata || length < 1) return( 0 );
    if (p0 == 0.0 && p1 == 0.0) return( 0 );

    ix1--;
    ix3--;

    rdata += ix1;
    idata += ix1;

    p0 = 2.0*PI*p0/360.0;
    p1 = 2.0*PI*p1/360.0;

    for( i = ix1; i < ix3; i++ )
       {
        phi   = p0 + p1*i/length;

        tCOS  = cos( (double) phi );
        tSIN  = sin( (double) phi );

        tR    = *rdata;
        tI    = *idata;

        *rdata++ = tCOS*tR - tSIN*tI;
        *idata++ = tCOS*tI + tSIN*tR;
       }

    return( 0 );
}

/***/
/* phaseRIP: phase a separated real/imaginary data according to a phase vector.
/***/

int phaseRIP64( rdata, idata, pdata, length )

   float     *rdata, *idata, *pdata;
   NMR_INT   length;
{
    float phi, tR, tI, tCOS, tSIN;
    int i;

    if (!rdata || !idata || !pdata || length < 1) return( 0 );

    for( i = 0; i < length; i++ )
       {
        tCOS  = cos( (double)pdata[i] );
        tSIN  = sin( (double)pdata[i] );

        tR    = *rdata;
        tI    = *idata;

        *rdata++ = tCOS*tR - tSIN*tI;
        *idata++ = tCOS*tI + tSIN*tR;
       }

    return( 0 );
}

/***/
/* N-point smooth (moving average).
/***/

int nSmooth64( src, dest, length, window )

   float   *src, *dest;
   NMR_INT length, window;
{
    float   *rPtrI, *rPtrJ, sum;
    NMR_INT i, j, count;

    if (!src || !dest || length < 1 || window < 1) return( 0 );

/***/
/* Filter head of data:
/***/

    rPtrI = src;

    for( i = 0; i < window; i++ )
       {
        rPtrJ = rPtrI;
        sum   = 0.0;

        for( j = 0; j < i + window; j++ ) sum += *rPtrJ++;

        *dest++ = sum/(float)(1 + i + window);
       }

/***/
/* Filter interior of data:
/***/

    rPtrI = src;
    count = 1 + 2*window;

    for( i = window; i < length - window; i++ )
       {
        rPtrJ = rPtrI++;
        sum   = 0.0;

        for( j = 0; j < count; j++ ) sum += *rPtrJ++;

        *dest++ = sum/(float)count;
       }

/***/
/* Filter tail of data:
/***/

    for( i = length - window; i < length; i++ )
       {
        rPtrI = src + i - window;
        rPtrJ = rPtrI;
        sum   = 0.0;

        for( j = i - window; j < length; j++ ) sum += *rPtrJ++;

        *dest++ = sum/(float)(window + length - i);
       }

    return( 0 );
}

/***/
/* N-Point centering (subtraction of local average) of data in vec. 64-bit.
/***/

int nCenter64( src, dest, length, window )

   float   *src, *dest;
   NMR_INT length, window;
{
    float   *rPtrI, *rPtrJ,*rPtrK,  sum;
    NMR_INT i, j, count;

    if (!src || !dest || length < 1) return( 0 );

    if (window < 1)
       {
        for( i = 0; i < length; i++ ) *dest++ = *src++;

        return( 0 );
       }

/***/
/* Center head of data:
/***/

    rPtrI = src;
    rPtrK = src;

    for( i = 0; i < window; i++ )
       {
        rPtrJ = rPtrI;
        sum   = 0.0;

        for( j = 0; j < i + window; j++ ) sum += *rPtrJ++;

        *dest++ = *rPtrK++ - sum/(float)(1 + i + window);
       }

/***/
/* Center interior of data:
/***/

    rPtrI = src;
    count = 1 + 2*window;

    for( i = window; i < length - window; i++ )
       {
        rPtrJ = rPtrI++;
        sum   = 0.0;

        for( j = 0; j < count; j++ ) sum += *rPtrJ++;

        *dest++ = *rPtrK ++ - sum/(float)count;
       }

/***/
/* Center tail of data:
/***/

    for( i = length - window; i < length; i++ )
       {
        rPtrI = src + i - window;
        rPtrJ = rPtrI;
        sum   = 0.0;

        for( j = i - window; j < length; j++ ) sum += *rPtrJ++;

        *dest++ = *rPtrK++ - sum/(float)(window + length - i);
       }

    return( 0 );
}

int iFindLoc( i, iList, iCount )

   int *iList, i, iCount;
{
    int j;

    if (!iList || iCount < 1) return( -1 );

    for( j = 0; j < iCount; j++ ) if (iList[j] == i) return( j );

    return( -1 );
}

/***/
/* Eliminate zeros from v1, and return count of non-zero elements.
/***/

NMR_INT vecDeZero64( v1, n, tol )

   float   *v1, tol;
   NMR_INT n;
{
    float   *src, *dest;
    NMR_INT i, k;

    src  = v1;
    dest = v1;

    k = 0;

    for( i = 0; i < n; i++ )
       {
        if (*dest < -tol || *dest > tol)
           {
            *src++ = *dest;
            k++;
           }

        dest++;
       }

    return( k );
}

/***/
/* Bottom.
/***/

