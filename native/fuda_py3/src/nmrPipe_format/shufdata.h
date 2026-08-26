/***/
/* shufdata.h: definitions for some useful data shuffling.
/***/

#define SHUF_NULL       0
#define SHUF_RI2C       1
#define SHUF_C2RI       2
#define SHUF_EXLR       3
#define SHUF_ROLR       4
#define SHUF_SWAP       5
#define SHUF_NEG_ALL    6 
#define SHUF_NEG_REAL   7
#define SHUF_NEG_IMAG   8
#define SHUF_NEG_LEFT   9
#define SHUF_NEG_RIGHT 10
#define SHUF_NEG_ALT   11 
#define SHUF_ABS       12 
#define SHUF_BYTESWAP  13
#define SHUF_SIGNUM    14
#define SHUF_R2I       15
#define SHUF_I2R       16
#define SHUF_RR2RI     17
#define SHUF_RI2RR     18
#define SHUF_RCPI      19
#define SHUF_RIIR      20
#define SHUF_RI2RI     21
#define SHUF_SORT      22
#define SHUF_ASORT     23
#define SHUF_FIXNAN    24

#define MC_MAG    1   /* Functions return real-only.    */
#define MC_POW    2
#define MC_ABS    3
#define MC_ADD    4
#define MC_MULT   5
#define MC_DIV    6

#define MC_INV    7   /* Functions retain complex data. */
#define MC_LOG    8
#define MC_EXP    9
#define MC_SEXP   10
#define MC_SQRT   11 
#define MC_SQUARE 12 
#define MC_ATAN   13 
#define MC_UNDER  14 

#define BINNING_NONE     0
#define BINNING_SUM      1
#define BINNING_AVG      2
#define BINNING_RMS      3
#define BINNING_MIN      4
#define BINNING_MAX      5
#define BINNING_MINABS   6
#define BINNING_MAXABS   7
#define BINNING_AMIN     8
#define BINNING_AMAX     9
#define BINNING_MEDIAN   10
#define BINNING_STDDEV   11
#define BINNING_VARIANCE 12
#define BINNING_RANGE    13

int shufData();
