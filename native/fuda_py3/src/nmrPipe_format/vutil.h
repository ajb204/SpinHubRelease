/***/
/* vutil.h: vector utilities.
/***/

#include "prec.h"

int iFindLoc();

int cSet64();
int iSet64();
int rSet64();

#define cSet( VEC, N, VAL ) cSet64( VEC, (NMR_INT)((NMR_INT)N), (char)(VAL) )
#define iSet( VEC, N, VAL ) iSet64( VEC, (NMR_INT)((NMR_INT)N), (int)(VAL) )

int vByteSwap64();
int iByteSwap64();

#define vByteSwap( VEC, N ) vByteSwap64( VEC, (NMR_INT)((NMR_INT)N) )
#define iByteSwap( VEC, N ) iByteSwap64( VEC, (NMR_INT)((NMR_INT)N) )

int getVec64();
int putVec64();

#define getVec( src, dest, eCount, eJump ) getVec64( src, dest, (NMR_INT)((NMR_INT)eCount), (NMR_INT)((NMR_INT)eJump) )
#define putVec( dest, src, eCount, eJump ) putVec64( dest, src, (NMR_INT)((NMR_INT)eCount), (NMR_INT)((NMR_INT)eJump) )

float vMax64();
float vMin64();
float vAMin64();
float vAMax64();
float vMinAbs64();
float vMaxAbs64();
float vRMS64();
float vVariance64();
float vStdDev64();
float vMedian64();
float vNoiseEst64();
float vRMS64();
float vSum64();
float vSumPos64();
float vSumNeg64();
float vNorm64();
float vSumSq64();

#define vMax( VEC, N )      vMax64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vMin( VEC, N )      vMin64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vAMin( VEC, N )     vAMin64(   VEC, (NMR_INT)((NMR_INT)N) )
#define vAMax( VEC, N )     vAMax64(   VEC, (NMR_INT)((NMR_INT)N) )
#define vMinAbs( VEC, N )   vMinAbs64( VEC, (NMR_INT)((NMR_INT)N) )
#define vMaxAbs( VEC, N )   vMaxAbs64( VEC, (NMR_INT)((NMR_INT)N) )
#define vRMS( VEC, N )      vRMS64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vSum( VEC, N )      vSum64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vSumPos( VEC, N )   vSumPos64( VEC, (NMR_INT)((NMR_INT)N) )
#define vSumNeg( VEC, N )   vSumNeg64( VEC, (NMR_INT)((NMR_INT)N) )
#define vNorm( VEC, N )     vNorm64(   VEC, (NMR_INT)((NMR_INT)N) )
#define vSumSq( VEC, N )    vSumSq64(  VEC, (NMR_INT)((NMR_INT)N) )

#define vAvg( VEC, N )      ((N) < 1 ? 0.0 : vSum64(    VEC, (NMR_INT)((NMR_INT)N) )/((float)(N)))

#define vVariance( VEC, N ) vVariance64(  VEC, (NMR_INT)((NMR_INT)N) )

float vvGetScale64();
float vvRMS64();
float vvScale64();
float vvDot64();
float vvPearsonR64();

#define vvGetScale( A, B, N ) vvGetScale64(  A, B, (NMR_INT)((NMR_INT)N) )
#define vvRMS( A, B, N )      vvRMS64(       A, B, (NMR_INT)((NMR_INT)N) )
#define vvScale( A, B, N )    vvScale64(     A, B, (NMR_INT)((NMR_INT)N) )
#define vvDot( A, B, N )      vvDot64(       A, B, (NMR_INT)((NMR_INT)N) )
#define vvPearsonR( A, B, N ) vvPearsonR64(  A, B, (NMR_INT)((NMR_INT)N) )

float vvGetScaleRI64();

#define vvGetScaleRI( A, B, N, SR, SI )  vvGetScaleRI64( A, B, (NMR_INT)((NMR_INT)N), SR, SI )

int vAlt64();
int vNeg64();
int vAbs64();
int vSign();
int vRev64();
int vSort64();
int vASort64();
int vFixNaN64();

#define vAlt( VEC, N )    vAlt64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vNeg( VEC, N )    vNeg64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vAbs( VEC, N )    vAbs64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vSign( VEC, N )   vSign64(   VEC, (NMR_INT)((NMR_INT)N) )
#define vRev( VEC, N )    vRev64(    VEC, (NMR_INT)((NMR_INT)N) )
#define vFixNaN( VEC, N ) vFixNaN64( VEC, (NMR_INT)((NMR_INT)N) )

#define vSort( VEC, N )   vSort64(  VEC, (NMR_INT)((NMR_INT)N) )
#define vASort( VEC, N )  vASort64( VEC, (NMR_INT)((NMR_INT)N) )

int vsAdd64();
int vsAddGRand64();
int vsMult64();
int vsSet64();
int vsClipAbove64();
int vsClipBelow64();

#define vsSet( VEC, N, VAL )  vsSet64(  VEC, (NMR_INT)((NMR_INT)N), (float)(VAL) )
#define vsAdd( VEC, N, VAL )  vsAdd64(  VEC, (NMR_INT)((NMR_INT)N), (float)(VAL) )
#define vsSub( VEC, N, VAL )  vsAdd64(  VEC, (NMR_INT)((NMR_INT)N), (float)(-1.0*(VAL)) )
#define vsMult( VEC, N, VAL ) vsMult64( VEC, (NMR_INT)((NMR_INT)N), (float)(VAL) )

#define vsAddGRand( VEC, N, VAL )  vsAddGRand64(  VEC, (NMR_INT)((NMR_INT)N), (float)(VAL) )

int vsMultThresh64();

#define vsMultThresh( VEC, N, VAL, T ) vsMultThresh64( VEC, (NMR_INT)((NMR_INT)N), (float)(VAL), (float)(T) )

int vqMult64();

#define vqMult( VEC, N, VALR, VALI ) vqMult64( VEC, (NMR_INT)((NMR_INT)N), (float)(VALR), (float)(VALI) )

int vvCopy64();
int vvAdd64();
int vvMult64();
int vvSub64();
int vvMask64();
int vvIMask64();
int vvSwap64();
int vvSetSign64();

#define vvCopy( DEST, SRC, N )     vvCopy64(     DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvAdd( DEST, SRC, N )      vvAdd64(      DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvMult( DEST, SRC, N )     vvMult64(     DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvDiv( DEST, SRC, N )      vvDiv64(      DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvSub( DEST, SRC, N )      vvSub64(      DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvMask( DEST, SRC, N )     vvMask64(     DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvIMask( DEST, SRC, N )    vvIMask64(    DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvSwap( DEST, SRC, N )     vvSwap64(     DEST, SRC, (NMR_INT)((NMR_INT)N) )
#define vvSetSign( DEST, SRC, N )  vvSetSign64(  DEST, SRC, (NMR_INT)((NMR_INT)N) )

int vvCopyOff64();

#define vvCopyOff( D, S, N, DOFF, SOFF ) vvCopyOff( D, S, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)DOFF), (NMR_INT)((NMR_INT)SOFF) )

int vsvAdd64();
int vsvMult64();
int vvsMult64();

#define vsvAdd(  DEST, SRC, N, VAL ) vsvAdd64(  DEST, SRC, (NMR_INT)((NMR_INT)N), (float)(VAL) )
#define vsvMult( DEST, SRC, N, VAL ) vsvMult64( DEST, SRC, (NMR_INT)((NMR_INT)N), (float)(VAL) )
#define vvsMult( DEST, SRC, N, VAL ) vvsMult64( DEST, SRC, (NMR_INT)((NMR_INT)N), (float)(VAL) )

#define rMove( SRC, DEST, N ) rMove64( SRC, DEST, (NMR_INT)((NMR_INT)N) )
#define rSet( VEC, N, VAL )   rSet64( VEC, (NMR_INT)((NMR_INT)N), (float)(VAL) )

int vvProd64(), vvvProd64(), vvvvProd64();

int coAdd64();
int vAbsCor64();
int vInterpol64();
int bilinInterp2D64();

#define vvProd( DEST, SRCX, SRCY, NX, NY ) vvProd64( DEST, SRCX, SRCY, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )

#define coAdd( VEC, N, A, NA ) coAdd64( VEC, (NMR_INT)((NMR_INT)N), A, (NMR_INT)((NMR_INT)NA) )

#define vAbsCor( DEST, SRC, BG, C, T, N ) vAbsCor64( DEST, SRC, BG, (float)(C), (float)(T), (NMR_INT)((NMR_INT)N) )

#define bilinInterp2D( SRC, NX_SRC, NY_SRC, DEST, NX_DEST, NY_DEST ) bilinInterp2D64( SRC, (NMR_INT)((NMR_INT)NX_SRC), (NMR_INT)((NMR_INT)NY_SRC), DEST, (NMR_INT)((NMR_INT)NX_DEST), (NMR_INT)((NMR_INT)NY_DEST) )

int nCenter64();
int nSmooth64();

#define nCenter( SRC, DEST, N, W ) nCenter64( SRC, DEST, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)W) )
#define nSmooth( SRC, DEST, N, W ) nSmooth64( SRC, DEST, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)W) )

int phaseRI64();
int phaseRIP64();

#define phaseRI( RD, ID, IX1, IXN, N, P0, P1 ) phaseRI64(  RD, ID, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)IXN), (NMR_INT)((NMR_INT)N), (float)(P0), (float)(P1) )
#define phaseRIP( RD, ID, PD, N )              phaseRIP64( RD, ID, PD, (NMR_INT)((NMR_INT)N) )

int dx();
int integ();
int dxV();
int dx2D();

int dx64();
int integ64();

int matInv64();
int svdLSMat();

int applyNorm64();

#define applyNorm( RD, IX1, IXN, KX1, KXN, NORMFN, NORMMODE, C ) applyNorm64( RD, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)IXN), (NMR_INT)((NMR_INT)KX1), (NMR_INT)((NMR_INT)KXN), NORMFN, NORMMODE, C )

NMR_INT vecDeZero64();

#define NORM_FUNC_NONE     0
#define NORM_FUNC_SUM      1
#define NORM_FUNC_AVG      2
#define NORM_FUNC_VARIANCE 3
#define NORM_FUNC_STDDEV   4
#define NORM_FUNC_RMS      5
#define NORM_FUNC_MIN      6
#define NORM_FUNC_MAX      7
#define NORM_FUNC_MINABS   8
#define NORM_FUNC_MAXABS   9
#define NORM_FUNC_AMIN    10
#define NORM_FUNC_AMAX    11
#define NORM_FUNC_RANGE   12
#define NORM_FUNC_MEDIAN  13
#define NORM_FUNC_NOISE   14 

#define NORM_MODE_NONE     0
#define NORM_MODE_ADD      1
#define NORM_MODE_SUBTRACT 2
#define NORM_MODE_MULT     3
#define NORM_MODE_DIVIDE   4
#define NORM_MODE_REPLACE  5

