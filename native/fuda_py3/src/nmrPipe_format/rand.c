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
/* gaussrand: returns a normally distributed random number with
/*            zero mean and unit variance. 
/*            Taken from numerical recipes.
/***/

static int    currentISeed = 0;   /* Current random number seed.              */
static int    gValid       = 0;   /* Whether gset is valid from last call.    */
static int    sValid       = 0;   /* Whether random number has been seeded.   */
static double gsetD        = 0.0; /* Prev float rand val inside unit circle.  */
static float  gsetR        = 0.0; /* Prev double rand val inside unit circle. */

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/time.h>

#include "rand.h"
#include "prec.h"

float gaussRand( )
{
    float  fac, r, v1, v2;

/***/
/* Find two (0,1) random numbers inside the unit circle.
/* These define a random angle, which is transformed to two normal deviates.
/* One is saved for use in the next call.
/***/

    if (gValid) 
       {
        gValid = 0;
        return( gsetR );
       }
    else
       {
        do 
           {
            v1 = 2.0*RANDOM1 - 1.0;
            v2 = 2.0*RANDOM1 - 1.0;
            r  = v1*v1 + v2*v2;
           } 
        while( r >= 1.0 );

        fac    = sqrt( (double) -2.0*log( (double) r)/r);
        gsetR  = v1*fac;
        gValid = 1;

        return( v2*fac );
       } 
}

/***/
/* gaussRandD: returns a normally distributed random number with
/*             zero mean and unit variance. 
/*             Taken from numerical recipes.
/***/

#include <math.h>
#include <stdlib.h>

#include "rand.h"

double gaussRandD( )
{
    double  fac, r, v1, v2;

/***/
/* Find two (0,1) random numbers inside the unit circle.
/* These define a random angle, which is transformed to two normal deviates.
/* One is saved for use in the next call.
/***/

    if (gValid) 
       {
        gValid = 0;
        return( gsetD );
       }
    else
       {
        do 
           {
            v1 = 2.0*DRANDOM1 - 1.0;
            v2 = 2.0*DRANDOM1 - 1.0;
            r  = v1*v1 + v2*v2;
           } 
        while( r >= 1.0 );

        fac    = sqrt( (double) -2.0*log( (double) r)/r);
        gsetD  = v1*fac;
        gValid = 1;

        return( v2*fac );
       } 
}

int vAddGNoise64( v, length, sigma )

   float   *v, sigma;
   NMR_INT length; 
{
    if (!v || length < 1 || sigma == 0.0) return( 0 );

    while( length-- ) *v++ += sigma*gaussRand();

    return( 0 );
}

#ifdef USE_OLD_RAND
int sRandInit( iseed )

   int iseed;
{
   int prevISeed;

   if (sValid) return( 0 );

   prevISeed = sRand( iseed );

   return( prevISeed );
}

int sRand( iseed )

   int iseed;
{
    int prevISeed;

    prevISeed = currentISeed;

    if (iseed == 0) iseed = (int) time( (time_t *)NULL );

    (void) srand( (unsigned)(iseed) );
    (void) srand48( (long)(iseed) );

    currentISeed = iseed;
    gValid       = 0;
    sValid       = 1;

    return( prevISeed );
}

#include "drand48.c"

#endif
