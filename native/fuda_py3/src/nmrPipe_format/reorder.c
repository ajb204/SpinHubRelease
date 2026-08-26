/******************************************************************************/
/*                                                                            */
/*                   ---- NIH NMR Software System ----                        */
/*                        Copyright 1992 and 1993                             */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/

/***/
/* reorder: perform FT reorder and related procedures on single vectors.
/***/

#include <stdio.h>
#include <math.h>

#include "prec.h"
#include "vutil.h"
#include "shufdata.h"
#include "reorder.h"

#ifdef S_SPLINT_S
#include <imsl.h>
#else
typedef struct {
    float       re;
    float       im;
} f_complex;
#endif

#define ICLIP( M, LO, HI )  if (M < LO) M = LO; if (M > HI) M = HI

/***/
/* reorder: exchange left and right halves of vec1.
/***/

int reorder( vec1, length, invFlag )

    float *vec1;
    int   length, invFlag;
{
    float *vec2, rtemp, zero;
    int   i;

    if (length < 1 || !vec1) return( 0 );

/***/
/* Even Point Count:
/*  Left and right sides are of equal size.
/*
/* Odd Point Count:
/*  Forward mode: left side is 1 point larger than right  ([0 -1 -2] [2 1]).
/*  Inverse mode: left side is 1 point smaller than right ([2 1] [0 -1 -2]).
/***/

    if (length % 2)
       {
        if (invFlag)
           {
            i     = length/2;
            vec2  = vec1 + length - 1;
            vec1 += i - 1; 
            zero  = *vec2;

            while( i-- )
               {
                rtemp   = *vec1;
                *vec1-- = *(vec2 - 1);
                *vec2-- = rtemp;
               }

            *vec2 = zero;
           }
        else
           {
            i    = length/2;
            zero = *vec1;
            vec2 = vec1 + i + 1;

            while( i-- )
               {
                rtemp   = *(vec1 + 1);
                *vec1++ = *vec2;
                *vec2++ = rtemp;
               }

            *vec1 = zero;
           }
       }
    else
       { 
        length /= 2;
        vec2    = vec1 + length;

        while( length-- )
           {
            rtemp   = *vec1;
            *vec1++ = *vec2;
            *vec2++ = rtemp;
           }
       }

    return( 0 );
}

/***/
/* c_reorder: exchange left and right halves of complex vector vec1.
/***/

int c_reorder( vec1, length, invFlag )

    f_complex *vec1;
    int   length, invFlag;
{
    f_complex *vec2, rtemp, zero;
    int       i;

    if (length < 1 || !vec1) return( 0 );

/***/
/* Even Point Count:
/*  Left and right sides are of equal size.
/*
/* Odd Point Count:
/*  Forward mode: left side is 1 point larger than right  ([0 -1 -2] [2 1]).
/*  Inverse mode: left side is 1 point smaller than right ([2 1] [0 -1 -2]).
/***/

    if (length % 2)
       {
        if (invFlag)
           {
            i     = length/2;
            vec2  = vec1 + length - 1;
            vec1 += i - 1; 
            zero  = *vec2;

            while( i-- )
               {
                rtemp   = *vec1;
                *vec1-- = *(vec2 - 1);
                *vec2-- = rtemp;
               }

            *vec2 = zero;
           }
        else
           {
            i    = length/2;
            zero = *vec1;
            vec2 = vec1 + i + 1;

            while( i-- )
               {
                rtemp   = *(vec1 + 1);
                *vec1++ = *vec2;
                *vec2++ = rtemp;
               }

            *vec1 = zero;
           }
       }
    else
       { 
        length /= 2;
        vec2    = vec1 + length;

        while( length-- )
           {
            rtemp   = *vec1;
            *vec1++ = *vec2;
            *vec2++ = rtemp;
           }
       }

    return( 0 );
}

/***/
/* Reverse left half, then right half, of vec.
/***/

int reorder2( vec, length )

    float *vec;
    int   length;
{
    if (length < 1 || !vec) return( 0 );

    (void) vRev( vec, (length + 1)/2 );
  
    vec    += (length + 1)/2;
    length -= (length + 1)/2;

    (void) vRev( vec, length );

    return( 0 );
}

/***/
/* Left shift data in vec by n points, padded by zeroes.
/***/

int lShift( vec, length, n )

   float *vec;
   int   length, n;
{
    float *src, *dest;
    int   length2;

    if (!vec || !n || length < 1) return( 0 );

    if (n < 0)  return( rShift( vec, length, -n ));

    if (n > length) n = length;

    dest    = vec;
    src     = vec + n;
    length2 = length - n;
   
    while( length2-- ) *dest++ = *src++;

    dest = vec + length - n;

    while( n-- ) *dest++ = 0.0;

    return( 0 ); 
}

/***/
/* Circular Left shift data in vec by n points, inverting wrapped data.
/***/

int clShift( vec, length, n, work, invFlag )

   float *vec, *work;
   int   length, invFlag, n;
{
    float *src, *dest;
    int   inv1Flag, inv2Flag, length2;

    if (!vec || !work || !n || length < 1) return( 0 );

    if (n < 0) return( crShift( vec, length, -n, work, invFlag ));

    inv1Flag = 0;
    inv2Flag = invFlag;

    while( n > length )
       {
        n -= length;

        if (invFlag)
           {
            inv1Flag = !inv1Flag;
            inv2Flag = !inv2Flag;
           }
       }

    dest    = work;
    src     = vec + n;
    length2 = length - n;

    if (inv1Flag)
       while( length2-- ) *dest++ = -*src++;
    else
       while( length2-- ) *dest++ = *src++;

    dest = work + length - n;
    src  = vec;

    if (inv2Flag)
       while( n-- ) *dest++ = -*src++;
    else
       while( n-- ) *dest++ = *src++;

    dest = vec;
    src  = work;

    while( length-- ) *dest++ = *src++; 

    return( 0 );
}

/***/
/* Right shift data in vec by n points, padded by zeroes.
/***/

int rShift( vec, length, n )

   float *vec;
   int   length, n;
{
    float *src, *dest; 
    int   length2;

    if (!vec || !n || length < 1) return( 0 );

    if (n < 0)  return( lShift( vec, length, -n ));

    if (n > length) n = length;

    dest    = vec + length - 1;
    src     = vec + length - n - 1;
    length2 = length - n;

    while( length2-- ) *dest-- = *src--; 

    dest = vec;
    
    while( n-- ) *dest++ = 0.0;

    return( 0 );
}

/***/
/* Circular Right shift data in vec by n points.
/***/

int crShift( vec, length, n, work, invFlag )

   float *vec, *work;
   int   length, invFlag, n;
{
    float *src, *dest;
    int   inv1Flag, inv2Flag, length2;

    if (!vec || !work || !n || length < 1) return( 0 );

    if (n < 0) return( clShift( vec, length, -n, work, invFlag ));

    inv1Flag = 0;
    inv2Flag = invFlag;

    while( n > length ) 
       {
        n -= length;

        if (invFlag)
           {
            inv1Flag = !inv1Flag;
            inv2Flag = !inv2Flag;
           }
       }

    dest    = work + length - 1;
    src     = vec + length - n - 1;
    length2 = length - n;

    if (inv1Flag)
       while( length2-- ) *dest-- = -*src--;
    else
       while( length2-- ) *dest-- = *src--;

    dest = work;
    src  = vec + length - n;  

    if (inv2Flag)
       while( n-- ) *dest++ = -*src++;
    else
       while( n-- ) *dest++ = *src++;

    dest = vec;
    src  = work;

    while( length-- ) *dest++ = *src++;

    return( 0 );
}

/***/
/* zfTail: zero fill at end of data.
/***/

int zfTail( vec, length, pad )

   float *vec;
   int   length, pad;
{
    if (!vec || pad < 1 || length < 1) return( 0 );

    vec += length;

    while( pad-- ) *vec++ = 0.0;

    return( 0 );
}

/***/
/* zfRep: repeat data instead of zero fill at end of data.
/***/

int zfRep( vec, length, pad )

   float *vec;
   int   length, pad;
{
    float *src, *dest;
    int   i, j;

    if (!vec || pad < 1 || length < 1) return( 0 );

    dest = vec + length;
    src  = vec;

    j = 0;

    for( i = 0; i < pad; i++ )
       {
        *dest++ = *src++;

        j++;

        if (j == length)
           {
            src = vec;
            j   = 0;
           }
       }

    return( 0 );
}

/***/
/* zfMid: zero fill in middle of data.
/***/

int zfMid( vec, length, pad )

   float *vec;
   int   length, pad;
{
    float *endPtr, *midPtr;
    int i, lSize, rSize;

    if (!vec || pad < 0 || length < 0) return( 0 );

    endPtr = vec + length + pad - 1;
    midPtr = vec + length - 1;

    lSize  = length/2;
    rSize  = length - lSize;
 
    for( i = 0; i < rSize; i++ ) *endPtr-- = *midPtr--;
   
    midPtr = vec + lSize;

    for( i = 0; i < pad; i++ ) *midPtr++ = 0.0;
 
    return( 0 );
}

/***/
/* zfInter: inter-point zero-fill.
/***/

int zfInter( vec, length, zfCount, zfSize )

   float *vec;
   int   length, zfCount, zfSize;
{
    float *src, *dest;
    int   i, j, n;

    if (!vec || length < 1 || zfSize < 1) return( 0 );

    if (zfCount < 1) 
       {
        (void) zfTail( vec, length, zfSize - length );
        return( 0 );
       }

    n    = length*(zfCount + 1);
    src  = vec + length - 1;
    dest = vec + n - zfCount  - 1;

    for( i = 0; i < length; i++ )
       {
        *dest = *src--;
        dest -= zfCount + 1;
       }

    dest = vec + 1;

    for( i = 0; i < length; i++ )
       {
        for( j = 0; j < zfCount; j++ ) *dest++ = 0.0;
        dest++;
       }

    dest = vec + n;

    for( i = n; i < zfSize; i++ ) *dest++ = 0.0;

    return( 0 );
}

/***/
/* Reverse left half, then right half, of complex vector vec.
/***/

int c_reorder2( vec, length )

    f_complex *vec;
    int       length;
{
    if (!vec || length < 1) return( 0 );

    (void) c_reverse( vec, (length + 1)/2 );
  
    vec    += (length + 1)/2;
    length -= (length + 1)/2;

    (void) c_reverse( vec, length );

    return( 0 );
}

/***/
/* Reverse complex vector vec.
/***/

int c_reverse( vec, length )

   f_complex *vec;
   int       length;
{
    f_complex *rPtr, rtemp;

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

/***/
/* Replace real part of vec in origin-1 range ix1,ix3 with modulus.
/***/

int c_modulus( vec, ix1, ix3 )

   f_complex *vec;
   int       ix1, ix3;
{
    int i;

    ix1--;
    ix3--;

    vec += ix1;

    for( i = ix1; i < ix3; i++ )
       {
        vec->re = sqrt( (double) vec->re*vec->re + vec->im*vec->im );
        vec++;
       }

    return( 0 );
}

/***/
/* Phase correct complex vec in origin-1 range ix1,ix3.
/***/

int c_phase( vec, ix1, ix3, length, p0, p1 )

   int       ix1, ix3, length;
   float     p0, p1;
   f_complex *vec;
{
    float phi, tR, tI, tCOS, tSIN;
    int i;

    ix1--;
    ix3--;

    vec += ix1;

    p0 = 2.0*PI*p0/360.0;
    p1 = 2.0*PI*p1/360.0;

    for( i = ix1; i < ix3; i++ )
       {
        phi   = p0 + p1*i/length;

        tCOS  = cos( (double) phi );
        tSIN  = sin( (double) phi );

        tR    = vec->re;
        tI    = vec->im;

        vec->re = tCOS*tR - tSIN*tI;
        vec->im = tCOS*tI + tSIN*tR;

        vec++;
       }

    return( 0 );
}

/***/
/* Scale complex vector vec by inverse length (used for inverse FT).
/***/

int c_scale( vec, length )

   f_complex *vec;
   int       length;
{
    float c;

    if (!vec || length < 1) return( 0 );

    c = 1.0/(float)length;

    while( length-- )
       {
        vec->re *= c;
        vec->im *= c;
        vec++;
       }

    return( 0 );
}

/***/
/* Scale real vector vec by inverse length (used for inverse FT).
/***/

int r_scale( vec, length )

   float *vec;
   int   length;
{
    float c;

    if (!vec || length < 1) return( 0 );

    c = 1.0/(float)length;

    while( length-- )
       {
        *vec++ *= c;
       }

    return( 0 );
}

/***/
/* Append Mirror Image; assumes vec can hold 2*length elements.
/***/

int mirror( vec, length, leftFlag, invLFlag, invRFlag, htFlag )

   int   length, leftFlag, invLFlag, invRFlag, htFlag;
   float *vec;
{
    float *srcPtr, *destPtr;
    int   i;

/***/
/* Special case: mirror image as for Hilbert Transform;
/*   Original data in center.
/***/

    if (!vec || length < 1) return( 0 );

    if (htFlag)
       {
        if (!invLFlag && !invRFlag)
           (void) htMirror2( vec, length );
        else
           (void) htMirror( vec, length );

        return( 0 );
       }
    
    if (leftFlag)
       {
        srcPtr  = vec + length - 1;
        destPtr = vec + length;

        for( i = 0; i < length; i++ ) *destPtr++ = *srcPtr--;
       }
    else
       {
        srcPtr  = vec;
        destPtr = vec + length;

        for( i = 0; i < length; i++ ) *destPtr++ = *srcPtr++;

        srcPtr  = vec + length;
        destPtr = vec + length - 1;

        for( i = 0; i < length; i++ ) *destPtr-- = *srcPtr++;
       }

    if (invLFlag)
       {
        srcPtr = vec;

        for( i = 0; i < length; i++ )
           {
            *srcPtr = -*srcPtr;
            srcPtr++;
           }
       }

    if (invRFlag)
       {
        srcPtr = vec + length;

        for( i = 0; i < length; i++ )
           {
            *srcPtr = -*srcPtr;
            srcPtr++;
           }
       }

    return( 0 );
}

/***/
/* lpMirror: Mirror image trick for Linear Prediction of zero-phase data.
/*           Assumes vec can hold 2*length-1 elements.
/***/

int lpMirror( vec, length, imagFlag )

   float *vec;
   int   length, imagFlag;
{
    float *srcPtr, *destPtr;
    int   i;

/***/
/* Data start as 1...N (N Points)
/* Data return as N...2 1 2...N (2N - 1 Points)
/* Imaginary reflection is negated.
/***/

    if (!vec || length < 1) return( 0 );

    srcPtr  = vec + length - 1;
    destPtr = vec + 2*length - 2;

    for( i = 0; i < length; i++ ) *destPtr-- = *srcPtr--;

    srcPtr  = vec + 2*length - 2;
    destPtr = vec;
    length--;

    if (imagFlag)
       {
        for( i = 0; i < length; i++ ) *destPtr++ = -(*srcPtr--);
       }
    else
       {
        for( i = 0; i < length; i++ ) *destPtr++ = *srcPtr--;
       }

    return( 0 );
}

/***/
/* htMirror: Mirror Image trick for Hilbert Transform;
/***/

int htMirror( vec, length )

   float *vec;
   int   length;
{
    float *srcPtr, *destPtr;
    int   i;

/***/
/* Input:                       (Left Half)(Right Half) 
/* Output: (-Reverse Right Half)(Left Half)(Right Half)(-Reverse Left Half)
/***/

    if (!vec || length < 1) return( 0 );

    destPtr = vec + 3*length/2 - 1;
    srcPtr  = vec + length - 1;

    for( i = 0; i < length; i++ ) *destPtr-- = *srcPtr--;

    destPtr = vec;
    srcPtr  = vec + length;

    for( i = 0; i < length/2; i++ ) *destPtr++ = -(*srcPtr++);

    destPtr = vec + 3*length/2;
    srcPtr  = vec + length/2;

    for( i = 0; i < length/2; i++ ) *destPtr++ = -(*srcPtr++);

    return( 0 );
}

int htMirror2( vec, length )

   float *vec;
   int   length;
{
    float *srcPtr, *destPtr;
    int   i;

/***/
/* Input:                      (Left Half)(Right Half)
/* Output: (Reverse Right Half)(Left Half)(Right Half)(Reverse Left Half)
/***/

    if (!vec || length < 1) return( 0 );

    destPtr = vec + 3*length/2 - 1;
    srcPtr  = vec + length - 1;

    for( i = 0; i < length; i++ ) *destPtr-- = *srcPtr--;

    destPtr = vec;
    srcPtr  = vec + length;

    for( i = 0; i < length/2; i++ ) *destPtr++ = *srcPtr++;

    destPtr = vec + 3*length/2;
    srcPtr  = vec + length/2;

    for( i = 0; i < length/2; i++ ) *destPtr++ = *srcPtr++;

    return( 0 );
}

/***/
/* htMirrorC: Mirror Image trick for Hilbert Transform: complex input;
/***/

int htMirrorC( vec, length )

   f_complex *vec;
   int       length;
{
    f_complex *srcPtr, *destPtr;
    int       i;

/***/
/* Input:                       (Left Half)(Right Half)
/* Output: (-Reverse Right Half)(Left Half)(Right Half)(-Reverse Left Half)
/***/

    if (!vec || length < 1) return( 0 );

    destPtr = vec + 3*length/2 - 1;
    srcPtr  = vec + length - 1;

    for( i = 0; i < length; i++ ) *destPtr-- = *srcPtr--;

    destPtr = vec;
    srcPtr  = vec + length;

    for( i = 0; i < length/2; i++ )
       {
        destPtr->re = -srcPtr->re;
        destPtr->im = 0.0;
        destPtr++;
        srcPtr++;
       }

    destPtr = vec + 3*length/2;
    srcPtr  = vec + length/2;

    for( i = 0; i < length/2; i++ )
       {
        destPtr->re = -srcPtr->re;
        destPtr->im = 0.0;
        destPtr++;
        srcPtr++;
       }

    return( 0 );
}

/***/
/* ri2c: interleave elements from separate real and imag vectors.
/***/

int ri2c( rdata, idata, cdata, length )

   float *rdata, *idata, *cdata;
   int   length;
{

/***/
/* Start from end, to allow in-place use as ri2c( rdata, idata, rdata, length )
/***/

    if (!rdata || !idata || !cdata || length < 1) return( 0 );

    cdata += 2*length - 1;
    rdata += length   - 1;
    idata += length   - 1;

    while( length-- )
       {
        *cdata-- = *idata--;
        *cdata-- = *rdata--;
       }

    return( 0 );
}

/***/
/* c2ri: separate interleaved elements into real and imag vectors.
/***/

int c2ri( cdata, rdata, idata, length )

   float *rdata, *idata, *cdata;
   int   length;
{

    if (!rdata || !idata || !cdata || length < 1) return( 0 );

    while( length-- )
       {
        *rdata++ = *cdata++;
        *idata++ = *cdata++;
       }

    return( 0 );
}

/***/
/* vR2I: convert elements of vec from real to integer.
/***/

int vR2I( fVec, length )

   float *fVec;
   int   length;
{
    union r2i { int i; float f; char s[4]; } *vec;
    int i;

    if (!vec || length < 1) return( 0 );

    vec = (union r2i *)fVec;

    while( length-- )
       {
        i      = (int)vec->f;
        vec->i = i;
        vec++;
       }

    return( 0 );
}

/***/
/* vI2R: convert elements of vec from integer to real.
/***/

int vI2R( fVec, length )

   float *fVec;
   int   length;
{
    union r2i { int i; float f; char s[4]; } *vec;
    float f;

    if (!fVec || length < 1) return( 0 );

    vec = (union r2i *)fVec;

    while( length-- )
       {
        f      = (float)vec->i;
        vec->f = f;
        vec++;
       }

    return( 0 );
}

/***/
/* vAlt2: sign-alternate interleaved complex data in vec (+ + - -).
/***/
 
int vAlt2( vec, length )
 
   float *vec;
   int   length;
{
    if (!vec || length < 1) return( 0 );

    length /= 4;
    vec    += 2;
 
    while( length-- )
       {
        *vec = -(*vec);
        vec++;

        *vec = -(*vec);
        vec += 3;
       }
 
    return( 0 );
}

/***/
/* nextPower2: returns power of 2 nearest to N.
/***/

int nextPower2( n )

   int n;
{
    int m;

    m = 1;

    while( n > m ) m *= 2;

    return( m );
}

/***/
/* roundReg: attempts to round region coords to achieve a final size;
/*           If n > 0: so that size is multiple of n.
/*           If n < 0: so that size is power of n.
/***/

int roundReg( ix1, ix3, xSize, n )

   int *ix1, *ix3, xSize, n;
{
    int extra, inSize, outSize, kx1, kx3, adjust1, adjust2;

/***/
/* Skip trivial rounding.
/* Find the desired final size after rounding (avoid empty loop for lint).
/***/

    if (n > -2 && n < 2) return( 0 );

    inSize = 1 + *ix3 - *ix1;
    extra  = 0;

    if (n > 1)
       {
        for( outSize = 0; outSize < inSize; outSize +=  n ) extra++;
       }
    else if (n < -1)
       {
        for( outSize = 1; outSize < inSize; outSize *= -n ) extra++;
       }
    else
       {
        outSize = inSize;
       }

/***/
/* Find left and right adjustments which center new region.
/* Find second-order adjustments which place region in range 1 to xSize.
/* Clip region bounds.
/***/

    extra   = outSize - inSize;
    adjust1 = extra/2;
    adjust2 = extra - adjust1;

    kx1 = *ix1 - adjust1;
    kx3 = *ix3 + adjust2;

    adjust1 = kx1 < 1     ? 1 - kx1     : 0;
    adjust2 = kx3 > xSize ? kx3 - xSize : 0;

    *ix1 = kx1 - adjust2;
    *ix3 = kx3 + adjust1;

    ICLIP( *ix1, 1, xSize );
    ICLIP( *ix3, 1, xSize );

    return( 0 );
}

/***/
/* hadamard: sequentially ordered Hadamard transform, power of 2 size only.
/***/

int hadamard( rdata, work, size, invFlag )

   float *rdata, *work;
   int   size, invFlag;
{
    float *buffer, *d, d1, d2, scale;
    int   i, j, k, jump1, jump2;

    if (!rdata || size < 2) return( 0 );
    if (size != nextPower2( size )) return( 1 );

    jump1 = 2;
    jump2 = 1;
    scale = 1.0/(float)size;

    while( jump1 <= size ) 
       {
        for( i = 0; i < size; i++ ) work[i] = rdata[i];

        for( i = 0; i < size; i += jump1 ) 
           {
            k = i;
            d = rdata + i;

            for( j = 0; j < jump1; j += 2) 
               {
                d1 = work[k];
                d2 = work[k+jump2];

                if (j & 0x2) 
                   {
                    d[j]   = d1 - d2;
                    d[j+1] = d1 + d2;
                   }
                else 
                   {
                    d[j]   = d1 + d2;
                    d[j+1] = d1 - d2;
                   }

                k++;
               }
           }

        jump2   = jump1;
        jump1 <<= 1;
       }

    if (invFlag) for( i = 0; i < size; i++ ) rdata[i] *= scale;

    return( 0 );
}

/***/
/* Left shift data in vec by n points, padded by zeroes. 64-bit.
/***/

int lShift64( vec, length, n )

   float   *vec;
   NMR_INT length, n;
{
    float   *src, *dest;
    NMR_INT length2;

    if (!vec || !n || length < 1) return( 0 );

    if (n < 0)  return( rShift64( vec, length, -n ));

    if (n > length) n = length;

    dest    = vec;
    src     = vec + n;
    length2 = length - n;
   
    while( length2-- ) *dest++ = *src++;

    dest = vec + length - n;

    while( n-- ) *dest++ = 0.0;

    return( 0 ); 
}

/***/
/* Circular Left shift data in vec by n points, inverting wrapped data. 64-bit.
/***/

int clShift64( vec, length, n, work, invFlag )

   float   *vec, *work;
   NMR_INT length, n;
   int     invFlag;
{
    float   *src, *dest;
    int     inv1Flag, inv2Flag;
    NMR_INT length2;

    if (!vec || !work || !n || length < 1) return( 0 );

    if (n < 0) return( crShift64( vec, length, -n, work, invFlag ));

    inv1Flag = 0;
    inv2Flag = invFlag;

    while( n > length )
       {
        n -= length;

        if (invFlag)
           {
            inv1Flag = !inv1Flag;
            inv2Flag = !inv2Flag;
           }
       }

    dest    = work;
    src     = vec + n;
    length2 = length - n;

    if (inv1Flag)
       while( length2-- ) *dest++ = -*src++;
    else
       while( length2-- ) *dest++ = *src++;

    dest = work + length - n;
    src  = vec;

    if (inv2Flag)
       while( n-- ) *dest++ = -*src++;
    else
       while( n-- ) *dest++ = *src++;

    dest = vec;
    src  = work;

    while( length-- ) *dest++ = *src++; 

    return( 0 );
}

/***/
/* Right shift data in vec by n points, padded by zeroes. 64-bit.
/***/

int rShift64( vec, length, n )

   float   *vec;
   NMR_INT length, n;
{
    float   *src, *dest; 
    NMR_INT length2;

    if (!vec || !n || length < 1) return( 0 );

    if (n < 0)  return( lShift64( vec, length, -n ));

    if (n > length) n = length;

    dest    = vec + length - 1;
    src     = vec + length - n - 1;
    length2 = length - n;

    while( length2-- ) *dest-- = *src--; 

    dest = vec;
    
    while( n-- ) *dest++ = 0.0;

    return( 0 );
}

/***/
/* Circular Right shift data in vec by n points. 64-bit.
/***/

int crShift64( vec, length, n, work, invFlag )

   float   *vec, *work;
   NMR_INT length, n;
   int     invFlag;
{
    float   *src, *dest;
    int     inv1Flag, inv2Flag;
    NMR_INT length2;

    if (!vec || !work || !n || length < 1) return( 0 );

    if (n < 0) return( clShift64( vec, length, -n, work, invFlag ));

    inv1Flag = 0;
    inv2Flag = invFlag;

    while( n > length ) 
       {
        n -= length;

        if (invFlag)
           {
            inv1Flag = !inv1Flag;
            inv2Flag = !inv2Flag;
           }
       }

    dest    = work + length - 1;
    src     = vec + length - n - 1;
    length2 = length - n;

    if (inv1Flag)
       while( length2-- ) *dest-- = -*src--;
    else
       while( length2-- ) *dest-- = *src--;

    dest = work;
    src  = vec + length - n;  

    if (inv2Flag)
       while( n-- ) *dest++ = -*src++;
    else
       while( n-- ) *dest++ = *src++;

    dest = vec;
    src  = work;

    while( length-- ) *dest++ = *src++;

    return( 0 );
}

/***/
/* nextPower264: returns power of 2 nearest to N. 64-bit.
/***/

NMR_INT nextPower264( n )

   NMR_INT n;
{
    NMR_INT m;

    m = 1;

    while( n > m ) m *= 2;

    return( m );
}

/***/
/* hadamard: sequentially ordered Hadamard transform, power of 2 size only. 64-bit.
/***/

int hadamard64( rdata, work, size, invFlag )

   float   *rdata, *work;
   NMR_INT size;
   int     invFlag;
{
    float   *buffer, *d, d1, d2, scale;
    NMR_INT i, j, k, jump1, jump2;

    if (!rdata || size < 2) return( 0 );
    if (size != nextPower264( size )) return( 1 );

    jump1 = 2;
    jump2 = 1;
    scale = 1.0/(float)size;

    while( jump1 <= size ) 
       {
        for( i = 0; i < size; i++ ) work[i] = rdata[i];

        for( i = 0; i < size; i += jump1 ) 
           {
            k = i;
            d = rdata + i;

            for( j = 0; j < jump1; j += 2) 
               {
                d1 = work[k];
                d2 = work[k+jump2];

                if (j & 0x2) 
                   {
                    d[j]   = d1 - d2;
                    d[j+1] = d1 + d2;
                   }
                else 
                   {
                    d[j]   = d1 + d2;
                    d[j+1] = d1 - d2;
                   }

                k++;
               }
           }

        jump2   = jump1;
        jump1 <<= 1;
       }

    if (invFlag) for( i = 0; i < size; i++ ) rdata[i] *= scale;

    return( 0 );
}

/***/
/* Bottom.
/***/
