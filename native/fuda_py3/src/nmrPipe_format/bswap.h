
#include "prec.h"

int vDbl2Flt_64(); 
int byteSwap2_64();
int byteSwap3_64();
int byteSwap4_64();
int byteSwap8_64();

int byteSwapN_64();

#define byteSwap2( BUFF, N ) byteSwap2_64( BUFF, (NMR_INT)((NMR_INT)N) )
#define byteSwap3( BUFF, N ) byteSwap3_64( BUFF, (NMR_INT)((NMR_INT)N) )
#define byteSwap4( BUFF, N ) byteSwap4_64( BUFF, (NMR_INT)((NMR_INT)N) )
#define byteSwap8( BUFF, N ) byteSwap8_64( BUFF, (NMR_INT)((NMR_INT)N) )
#define vDbl2Flt( BUFF, N )  vDbl2Flt_64( BUFF, (NMR_INT)((NMR_INT)N) )

#define byteSwapN( BUFF, N, M ) byteSwapN_64( BUFF, (NMR_INT)((NMR_INT)N), (NMR_INT)((NMR_INT)M) )


