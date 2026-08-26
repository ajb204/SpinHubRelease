#include "prec.h"

int getMin64();
int getMax64();
int getAMin64();
int getAMax64();
int getMinAbs64();
int getMaxAbs64();
int getAvg64();
int getSum64();
int getSumSq64();
int getSumAbs64();
int getRMS64();
int getNorm64();
int getVariance64();
int getStdDev64();
int getRange64();
int getMedian64();
int getCenterOfMass64();

#define getMin( VEC, N, VPTR )      getMin64(    VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getMax( VEC, N, VPTR )      getMax64(    VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getAMin( VEC, N, VPTR )     getAMin64(   VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getAMax( VEC, N, VPTR )     getAMax64(   VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getMinAbs( VEC, N, VPTR )   getMinAbs64( VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getMaxAbs( VEC, N, VPTR )   getMaxAbs64( VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getAvg( VEC, N, VPTR )      getAvg64(    VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getSum( VEC, N, VPTR )      getSum64(    VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getSumSq( VEC, N, VPTR )    getSumSq64(  VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getSumAbs( VEC, N, VPTR )   getSumAbs64( VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getRMS( VEC, N, VPTR )      getRMS64(    VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getNorm( VEC, N, VPTR )     getNorm64(   VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getStdDev( VEC, N, VPTR )   getStdDev64( VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getRange( VEC, N, VPTR )    getRange64(  VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getMedian( VEC, N, VPTR )   getMedian64( VEC, (NMR_INT)((NMR_INT)N), VPTR )
#define getVariance( VEC, N, VPTR ) getVariance64( VEC, (NMR_INT)((NMR_INT)N), VPTR )

int medianW64();
int getHist64();
int getNoise64();
int getBase64();

#define medianW( V, N, NW, W )     medianW64( V, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)NW), W )
#define getHist( Y, NY, X, H, NH ) getHist64( Y, (NMR_INT)((NMR_INT)NY), X, H, (NMR_INT)((NMR_INT)NH) )

#define getNoise( RD, WORK, N, NW, FRAC, NI, SPTR ) getNoise64( RD, WORK, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)NW), (float)(FRAC), NI, SPTR )

#define getBase( RD, N, FRAC, BPTR ) getBase64( RD, (NMR_INT)((NMR_INT)N), (float)(FRAC), BPTR )

int minMax64();
int minMax264();
int minMaxJ64();

#define minMax(  VEC, N, MINPTR, MAXPTR ) minMax64(  VEC, (NMR_INT)((NMR_INT)N), MINPTR, MAXPTR )
#define minMax2( VEC, N, MINPTR, MAXPTR ) minMax264( VEC, (NMR_INT)((NMR_INT)N), MINPTR, MAXPTR )

#define minMaxJ( VEC, N, COUNT, JUMP, MINPTR, MAXPTR ) minMaxJ64( VEC, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)COUNT), (NMR_INT)((NMR_INT)JUMP), MINPTR, MAXPTR )

int getIMax3();
int getIMax4();
int getIMin3();
int getIMin4();

