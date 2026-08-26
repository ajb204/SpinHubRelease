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
/* getStat: some useful simple statistics.
/***/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "prec.h"
#include "getstat.h"
#include "vutil.h"

#define FABS(X) ((X) > 0.0 ? (X) : (-(X)))

static int   qTestSpecVal();
static float getDevNum();

int getNorm64( vec, length, norm )

   float   *vec, *norm;
   NMR_INT length;
{
    *norm = 0.0;

    if (!vec || length < 1) return( 1 );

    while( length-- )
       {
        *norm += (*vec)*(*vec);
        vec++;
       }

   return( 0 );
}

int getSumSq64( vec, length, dev )

   float   *vec, *dev;
   NMR_INT length;
{
    NMR_INT i;

    *dev = 0.0;

    if (!vec || length < 1) return( 1 );

    for( i = 0; i < length; i++ )
       {
        *dev += (*vec)*(*vec);
        vec++;
       }

    return( 0 );
}

int getSumAbs64( vec, length, dev )

   float   *vec, *dev;
   NMR_INT length;
{
    NMR_INT i;

    *dev = 0.0;

    if (!vec || length < 1) return( 1 );

    for( i = 0; i < length; i++ )
       {
        *dev += FABS(*vec);
        vec++;
       }

    return( 0 );
}

int getSum64( vec, length, sum )

   float   *vec, *sum;
   NMR_INT length;
{
    *sum = 0.0;

    if (!vec || length < 1) return( 1 );

    while( length-- )
       {
        *sum += *vec++;
       }

    return( 0 );
}

int getAvg64( vec, length, avg )

   float   *vec, *avg;
   NMR_INT length;
{
    NMR_INT i;
    float   r;

    *avg = 0.0;

    if (!vec || length < 1) return( 1 );

    r = 1.0/(float)length;

    for( i = 0; i < length; i++ )
       {
        *avg += (*vec++)*r;
       }

    return( 0 );       
}

int getStdDev64( vec, length, dev )

   float   *vec, *dev;
   NMR_INT length;
{
    float    r, avg;
    NMR_INT  i;

    *dev = 0.0;

    if (!vec || length < 1) return( 1 );

    if (length == 1)
       {
        *dev = *vec;
       }
    else
       {
        (void) getAvg64( vec, length, &avg );

        for( i = 0; i < length; i++ )
           {
            r     = avg - *vec++; 
            *dev += r*r;
           }

        *dev = sqrt( (double) *dev/(float)(length - 1) );
       }

    return( 0 );
} 

int getVariance64( vec, length, dev )

   float   *vec, *dev;
   NMR_INT length;
{
    float    r, avg;
    NMR_INT  i;

    *dev = 0.0;

    if (!vec || length < 1) return( 1 );

    if (length == 1)
       {
        *dev = *vec;
       }
    else
       {
        (void) getAvg64( vec, length, &avg );

        for( i = 0; i < length; i++ )
           {
            r     = avg - *vec++;
            *dev += r*r;
           }

        *dev /= (float)length;
       }

    return( 0 );
}

int getRMS64( vec, length, dev )

   float   *vec, *dev;
   NMR_INT length;
{
    NMR_INT i;

    *dev = 0.0;

    if (!vec || length < 1) return( 1 );

    if (length == 1)
       {
        *dev = *vec;
       }
    else
       {
        for( i = 0; i < length; i++ )
           {
            *dev += (*vec)*(*vec);
            vec++;
           }

        *dev = sqrt( (double)(*dev)/(double)length );
       }

    return( 0 );
}

int getMin64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    *val = 0.0; 

    if (!vec || length < 1) return( 1 );

    *val = *vec;

    while( length-- )
       {
        if (*vec < *val) *val = *vec;
        vec++;
       }

    return( 0 );
}

int getMax64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    *val = 0.0; 

    if (!vec || length < 1) return( 1 );

    *val = *vec;

    while( length-- )
       {
        if (*vec > *val) *val = *vec;
        vec++;
       }

    return( 0 );
}

int getAMin64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    *val = 0.0; 

    if (!vec || length < 1) return( 1 );

    *val = *vec;

    while( length-- )
       {
        if (FABS(*vec) < FABS(*val)) *val = *vec;
        vec++;
       }

    return( 0 );
}

int getAMax64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    *val = 0.0; 

    if (!vec || length < 1) return( 1 );

    *val = *vec;

    while( length-- )
       {
        if (FABS(*vec) > FABS(*val)) *val = *vec;
        vec++;
       }

    return( 0 );
}

int getMinAbs64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{  
    *val = 0.0;
    
    if (!vec || length < 1) return( 1 );
    
    *val = FABS(*vec);
    
    while( length-- )
       {
        if (FABS(*vec) < *val) *val = FABS(*vec);
        vec++;
       }

    return( 0 );
}

int getMaxAbs64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    *val = 0.0;

    if (!vec || length < 1) return( 1 );

    *val = FABS(*vec);

    while( length-- )
       {
        if (FABS(*vec) > *val) *val = FABS(*vec);
        vec++;
       }

    return( 0 );
}

int getRange64( vec, length, val )

   float   *vec, *val;
   NMR_INT length;
{
    float minVal, maxVal;

    *val = 0.0; 

    if (!vec || length < 1) return( 1 );

    minVal = *vec;
    maxVal = *vec;

    while( length-- )
       {
        if (*vec < minVal) minVal = *vec;
        if (*vec > maxVal) maxVal = *vec;
        vec++;
       }

    *val = maxVal - minVal;

    return( 0 );
}

/***/
/* medianW: replace data with local median.
/***/

int medianW64( vec, length, nw, work )

   float   *vec, *work;
   NMR_INT length, nw;
{
    NMR_INT i, j, k, first, last;
    float   *w1, *w2;

    if (nw < 1) return( 0 );
    if (nw > length) return( 0 );

    w1 = work;
    w2 = w1 + length;

    for( i = 0; i < length; i++ )
       {
        first = i - nw;
        last  = i + nw;
        k     = 0;

        if (first < 0) first = 0;
        if (last  > length - 1) last = length - 1;

        for( j = first; j <= last; j++ ) w1[k++] = vec[j];

        (void) getMedian( w1, 1 + last - first, &w2[i] ); 
       }

    for( i = 0; i < length; i++ ) vec[i] = w2[i];

   return( 0 );
}

/***/
/* getMedian: estimate of the median; input data is sorted on return.
/***/

int getMedian64( ra, n, rmed )

   float   *ra, *rmed;
   NMR_INT n;
{
    NMR_INT l, j, ir, i;
    float   rra;

    *rmed = 0.0;

    if (!ra || n < 1) return( 1 );

    ra--;

    if (n < 2)
       {
        *rmed = ra[(n+1)/2];
        return( 0 );
       }
 
    l  = (n >> 1)+1;
    ir = n;

    for (;;)
       {
        if (l > 1)
            rra=ra[--l];
        else 
           {
            rra=ra[ir];
            ra[ir]=ra[1];

            if (--ir == 1)
               {
                ra[1] = rra;
                *rmed = ra[(n+1)/2];
                return( 0 );
               }
        }

        i=l;
        j=l << 1;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j+1]) ++j;
            if (rra < ra[j]) {
                ra[i]=ra[j];
                j += (i=j);
            }
            else j=ir+1;
        }
        ra[i]=rra;
    }

   return( 0 );
}

/***/
/* getCenterOfMass: return coords of center of mass.
/***/

int getCenterOfMass64( rdata, sizeList, centerList, dimCount )

   float *rdata, *centerList;
   int   *sizeList, dimCount;
{
    float   *rPtr, sumX, sumY, sumZ, sumA, sumI;
    NMR_INT ix, iy, iz, ia;

    (void) vsSet( centerList, dimCount, (float)0.0 );

    if (dimCount < 1 || dimCount > 4) return( 1 );

    sumX = 0.0;
    sumY = 0.0;
    sumZ = 0.0;
    sumA = 0.0;
    sumI = 0.0;

    rPtr = rdata;

    if (dimCount == 1)
       {
        if (sizeList[0] < 1) return( 1 );

        for( ix = 0; ix < sizeList[0]; ix++ )
           {
            sumX += (*rPtr)*ix;
            sumI += *rPtr++;
           }

        centerList[0] = 1.0 + (sumI == 0.0 ? 0.0 : sumX/sumI); 
       }
    else if (dimCount == 2)
       {
        if (sizeList[0] < 1) return( 1 );
        if (sizeList[1] < 1) return( 1 );

        for( iy = 0; iy < sizeList[1]; iy++ )
           {
            for( ix = 0; ix < sizeList[0]; ix++ )
               {
                sumX += (*rPtr)*ix;
                sumY += (*rPtr)*iy;
                sumI += *rPtr++;
               }
           }

        centerList[0] = 1.0 + (sumI == 0.0 ? 0.0 : sumX/sumI); 
        centerList[1] = 1.0 + (sumI == 0.0 ? 0.0 : sumY/sumI); 
       }
    else if (dimCount == 3)
       {
        if (sizeList[0] < 1) return( 1 );
        if (sizeList[1] < 1) return( 1 );
        if (sizeList[2] < 1) return( 1 );

        for( iz = 0; iz < sizeList[2]; iz++ )
           {
            for( iy = 0; iy < sizeList[1]; iy++ )
               {
                for( ix = 0; ix < sizeList[0]; ix++ )
                   {
                    sumX += (*rPtr)*ix;
                    sumY += (*rPtr)*iy;
                    sumZ += (*rPtr)*iz;
                    sumI += *rPtr++;
                   }
               }
           }

        centerList[0] = 1.0 + (sumI == 0.0 ? 0.0 : sumX/sumI);
        centerList[1] = 1.0 + (sumI == 0.0 ? 0.0 : sumY/sumI);
        centerList[2] = 1.0 + (sumI == 0.0 ? 0.0 : sumZ/sumI);
       }
    else
       {
        if (sizeList[0] < 1) return( 1 );
        if (sizeList[1] < 1) return( 1 );
        if (sizeList[2] < 1) return( 1 );
        if (sizeList[3] < 1) return( 1 );

        for( ia = 0; ia < sizeList[3]; iz++ )
           {
            for( iz = 0; iz < sizeList[2]; iz++ )
               {
                for( iy = 0; iy < sizeList[1]; iy++ )
                   {
                    for( ix = 0; ix < sizeList[0]; ix++ )
                       {
                        sumX += (*rPtr)*ix;
                        sumY += (*rPtr)*iy;
                        sumZ += (*rPtr)*iz;
                        sumA += (*rPtr)*ia;
                        sumI += *rPtr++;
                       }
                   }
               }
           }

        centerList[0] = 1.0 + (sumI == 0.0 ? 0.0 : sumX/sumI);
        centerList[1] = 1.0 + (sumI == 0.0 ? 0.0 : sumY/sumI);
        centerList[2] = 1.0 + (sumI == 0.0 ? 0.0 : sumZ/sumI);
        centerList[3] = 1.0 + (sumI == 0.0 ? 0.0 : sumA/sumI);
       }


    return( 0 );
}

/***/
/* getHist: return a histogram of yArray in hist; returns max histogram value.
/***/

int getHist64( yArray, yLength, xArray, hist, histLength )

   float   *yArray, *xArray, *hist;
   NMR_INT yLength, histLength;
{
    float   maxVal, minVal, range;
    NMR_INT i, histLoc, histMax;

    maxVal = *yArray;
    minVal = *yArray;

    for( i = 0; i < yLength; i++ )
       {
        if (maxVal < *(yArray + i)) maxVal = *(yArray + i);
        if (minVal > *(yArray + i)) minVal = *(yArray + i);
       }

    if (maxVal != minVal)
       range = maxVal - minVal;
    else
       range = 1.0;

    for( i = 0; i < histLength; i++ )
       {
        *(xArray + i) = minVal + range*(float)i/histLength;
       }

    for( i = 0; i < histLength; i++ )
       {
        *(hist + i) = 0.0;
       }

    histMax = 0;

    for( i = 0; i < yLength; i++ )
       {
        histLoc = histLength*(*(yArray++) - minVal)/range;

        if (histLoc < 0) histLoc = 0;
        if (histLoc > histLength - 1) histLoc = histLength - 1;

        *(hist + histLoc) += 1.0;

        if (*(hist + histLoc) > histMax) histMax = *(hist + histLoc);
       }

    return( histMax );
}

int getIMin3( i, j, k )

   int i, j, k;
{
    int minVal;

    minVal = i;

    if (j < minVal) minVal = j;
    if (k < minVal) minVal = k;

    return( minVal );
}

int getIMax3( i, j, k )

   int i, j, k;
{
    int maxVal;

    maxVal = i;

    if (j > maxVal) maxVal = j;
    if (k > maxVal) maxVal = k;

    return( maxVal );
}

int getIMin4( i, j, k, l )

   int i, j, k, l;
{
    int minVal;

    minVal = i;

    if (j < minVal) minVal = j;
    if (k < minVal) minVal = k;
    if (l < minVal) minVal = l;

    return( minVal );
}

int getIMax4( i, j, k, l )

   int i, j, k, l;
{
    int maxVal;

    maxVal = i;

    if (j > maxVal) maxVal = j;
    if (k > maxVal) maxVal = k;
    if (l > maxVal) maxVal = l;

    return( maxVal );
}

/***/
/* minMax: return highest and lowest values in vec.
/***/

int minMax64( vec, length, minVal, maxVal )

   float   *vec, *minVal, *maxVal;
   NMR_INT length;

{
    *minVal = *vec;
    *maxVal = *vec;

    while( length-- )
       {
        *minVal = (*vec < *minVal) ? *vec : *minVal;
        *maxVal = (*vec > *maxVal) ? *vec : *maxVal;
        vec++;
       }

    return( 0 );
}

/***/
/* minMax2: return highest and lowest values in vec without first initializing maxVal and minVal.
/***/

int minMax264( vec, length, minVal, maxVal )

   float   *vec, *minVal, *maxVal;
   NMR_INT length;

{
    while( length-- )
       {
        *minVal = (*vec < *minVal) ? *vec : *minVal;
        *maxVal = (*vec > *maxVal) ? *vec : *maxVal;
        vec++;
       }

    return( 0 );
}

/***/
/* minMaxJ: return highest and lowest values in dispersed vectors of matrix.
/***/

int minMaxJ64( matrix, size, count, jump, minVal, maxVal )

   float   *matrix, *minVal, *maxVal;
   NMR_INT size, count, jump; 
{
    *minVal = *matrix;
    *maxVal = *matrix;

    while( count-- )
       {
        (void) minMax264( matrix, size, minVal, maxVal );
        matrix += jump;
       }

    return( 0 );
}

/***/
/* getNoise: estimate for standard deviatation of 1D spectral noise.
/***/

int getNoise64( rdata,       /* Spectral data to analyze.                     */
                work,        /* Work array for analysis.                      */
                length,      /* Number of spectral data points.               */
                window,      /* Length of averaging window.                   */
                fraction,    /* Fraction of points to consider as noise.      */
                iterCount,   /* Use additional iterations to refine estimate. */
                noise )      /* On return, the estimated noise value.         */

   float   *rdata, *work, *noise, fraction;
   NMR_INT length, window;
   int     iterCount;
{
   NMR_INT loc, newLength;

   static float lastFrac = 0.0, devNum = 0.0;

/***/
/* Find the number of standard deviations associated with the given fraction.
/***/

    *noise = 0.0;

    if (!rdata || !work || length < 1) return( 0 );

    if (fraction != lastFrac)
       {
        devNum   = getDevNum( fraction );
        lastFrac = fraction;
       }

/***/
/* Center the data (subtract neighborhood average).
/* Sort the adjusted data by absolute value.
/* Estimate the std deviation of the noise from the last intensity.
/***/

    (void) nCenter64( rdata, work, length, window );

    newLength = vecDeZero64( work, length, (float)1.0e-16 );
    if (!newLength) return( 0 );

    loc = fraction*newLength;
       
    (void) qsort( (void *)work, (size_t)newLength, sizeof(float), qTestSpecVal );

    *noise = devNum <= 0.0 ? 0.0 : FABS( work[loc]/devNum );

    return( 0 );
}

/***/
/* getBase: estimate of zero-order baseline for 1D spectral noise.
/***/

int getBase64( rdata,       /* Spectral data to analyze; destroyed.          */
               length,      /* Number of spectral data points.               */
               fraction,    /* Fraction of points to consider as noise.      */
               base )       /* On return, the estimated baseline value.      */

   float   *rdata, *base, fraction;
   NMR_INT length;
{
   NMR_INT count;

/***/
/* Sort the adjusted data by absolute value.
/* Estimate the baseline as the average.
/***/

    count = fraction*length;

    (void) qsort( (void *)rdata, (size_t)length, sizeof(float), qTestSpecVal );
    (void) getAvg64( rdata, count, base );

    return( 0 );
}

/***/
/* getDevNum: estimate the number of standard deviations which encompasses
/*            the given fraction of data; this is a crude calculation of
/*            integrated error function.
/***/

static float getDevNum( fraction )

   float fraction;
{
    float x, dx, s, r;
    int   i;

    if (fraction < 0.0) fraction = 0.0;
    if (fraction > 1.0) fraction = 1.0;

    s  = 0.00;
    x  = 0.00;
    dx = 0.01;

    r = 2.0*dx/sqrt( 2.0*3.1415 );

    for( i = 0; i < 500; i++ )
       {
        s += r*exp( -0.5*x*x );
        if (s >= fraction) return( x );
        x += dx;
       }

    return( x ); 
}

/***/
/* qTestSpecVal: function used by sorting module.
/***/

static int qTestSpecVal( val1, val2 )

   float *val1, *val2;
{
    float r1, r2;

    r1 = FABS( *val1 );
    r2 = FABS( *val2 );

    if ( r1 > r2 )
       return( 1 );
    else if ( r1 < r2 )
       return( -1 );
    else
       return( 0 );
}

/***/
/* Bottom.
/***/
