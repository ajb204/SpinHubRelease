
int mv_read64();
int mv_write64();
int mv_move64();

#define mv_read(  FP,  EC, OFF, VC, VJ, VEC ) mv_read64(  FP,  (NMR_INT)((NMR_INT)EC), (NMR_INT)((NMR_INT)OFF), (NMR_INT)((NMR_INT)VC), (NMR_INT)((NMR_INT)VJ), VEC )
#define mv_write( FP,  EC, OFF, VC, VJ, VEC ) mv_write64( FP,  (NMR_INT)((NMR_INT)EC), (NMR_INT)((NMR_INT)OFF), (NMR_INT)((NMR_INT)VC), (NMR_INT)((NMR_INT)VJ), VEC )
#define mv_move(  SRC, EC, OFF, VC, VJ, VEC ) mv_move64(  SRC, (NMR_INT)((NMR_INT)EC), (NMR_INT)((NMR_INT)OFF), (NMR_INT)((NMR_INT)VC), (NMR_INT)((NMR_INT)VJ), VEC ) 

