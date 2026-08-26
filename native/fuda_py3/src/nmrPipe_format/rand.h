/***/
/* Random number generation.
/***/

#include <stdlib.h>

#if defined (WINNT) || defined(WINXP) || defined(SGI)
#define USE_OLD_RAND
#endif

#ifdef USE_OLD_RAND
#ifndef RAND_MAX 
#define RAND_MAX  2147483648.0 
#endif
#endif

#ifdef  IBM_FROM_HELL
#define RAND_MAX 32768.0
#endif

#ifdef CUBE
#undef  RAND_MAX
#define RAND_MAX 32768.0
#endif

#ifdef USE_OLD_RAND
#define SRAND(ISEED)     (void)sRand(ISEED)
#define RANDOM1          (((float)rand())/(float)(RAND_MAX+1.0))
#define DRANDOM1         drand48()
#else
#define SRAND(ISEED)     (void)srandom((unsigned int)(ISEED))
#define RANDOM1          (((float)random())/(float)(RAND_MAX+1.0))
#define DRANDOM1         (((double)random())/(double)(RAND_MAX+1.0))
#endif

/***/
/* Generates normally distributed random numbers with unit variance:
/***/

double gaussRandD();
float  gaussRand();
int    vAddGNoise64();

#define vAddGNoise( V, N, SIGMA ) vAddGNoise64( V, (NMR_INT)((NMR_INT)N), SIGMA )

int drand48Dummy();

#ifdef USE_OLD_RAND
int    sRand(), sRandInit();
double drand48();

double drand48();
long   irand48();
long   krand48();
void   lcong48();
long   lrand48();
long   mrand48();
void   srand48();

unsigned short *seed48();
#endif

