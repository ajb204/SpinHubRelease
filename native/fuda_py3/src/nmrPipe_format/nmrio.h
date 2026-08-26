/***/
/* Some definitions for I/O using NMRPipe-format files.
/***/

#include "prec.h"

int rdNMRi();
int rdNMRi3D();
int rdNMRi4D();
int rdNMRiU();
int rdNMRiZ();
int rdNMRiZU();

int getNMRMinMax();
int scaleNMR();
int scaleNMR2();
int scaleNMR3D();

int wrNMRi();
int wrNMRi3D();
int wrNMRi4D();
int wrNMRiU();

int truncNMR();

int rdNMRPipeData64(), wrNMRPipeData64();

#define rdNMRPipeData( NAME, A, OFF, N ) rdNMRPipeData64( NAME, A, (NMR_INT)((NMR_INT)OFF), (NMR_INT)((NMR_INT)N) )
#define wrNMRPipeData( NAME, A, OFF, N ) wrNMRPipeData64( NAME, A, (NMR_INT)((NMR_INT)OFF), (NMR_INT)((NMR_INT)N) )

#define rdNMRi( VEC, NAME, IX1, NX, IY1, NY ) rdNMRi64( VEC, NAME, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)IY1), (NMR_INT)((NMR_INT)NY) )
#define wrNMRi( VEC, NAME, IX1, NX, IY1, NY ) wrNMRi64( VEC, NAME, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)IY1), (NMR_INT)((NMR_INT)NY) )

#define rdNMRiU( VEC, FP, IX1, NX, IY1, NY, FD ) rdNMRiU64( VEC, FP, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)IY1), (NMR_INT)((NMR_INT)NY), FD )
#define wrNMRiU( VEC, FP, IX1, NX, IY1, NY, FD ) wrNMRiU64( VEC, FP, (NMR_INT)((NMR_INT)IX1), (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)IY1), (NMR_INT)((NMR_INT)NY), FD )
