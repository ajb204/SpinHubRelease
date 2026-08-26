/***/
/* varian.h: some definitions for conversion of Varian data.
/***/

static int VHDR_BEGIN = 8;
static int VHDR_MID   = 7;

#ifdef CRAY
#define VARDATA(ISIGN) \
  (i2rFlag ? (signFlag?ISIGN:1)*getVInt( varMat ) : \
             (signFlag?ISIGN:1)*getVFloat( varMat ))

typedef char VMATTYPE;
#define VARJUMP(J)     varMat += VINTSIZE*J
#define VARINCR        varMat += VINTSIZE;
#define VARSWAP4       SWAP4( varMat )
#else
#define VARDATA(ISIGN) \
  (i2rFlag ? (signFlag?ISIGN:1)*varMat->i : (signFlag?ISIGN:1)*varMat->r)

typedef union b2f { int i; float r; char c[sizeof(int)]; } VMATTYPE;
#define VARJUMP(J)     varMat += J
#define VARINCR        varMat++
#define VARSWAP4       SWAP4( varMat->c )
#endif

#define ISCALE     4294967296
#define ISCALE64K  281474976710656
#define VHDRSIZE   32

#ifdef CRAY
#define VINTSIZE   4
#define VFLOATSIZE 2
#define VSHORTSIZE 2
#else
#define VINTSIZE   sizeof(int)
#define VFLOATSIZE sizeof(float)
#define VSHORTSIZE sizeof(short)
#endif

struct VHdrBuff
   {
    int     v1, v2, v3, v4;
    char    v5[10];
    short   v6;
    char    v7[4];
   };

struct VBlockHdrBuff 
   {
    short scale;
    short status;
    short index;
    short mode;
    int   ctcount;
    float lpval;
    float rpval;
    float lvl;
    float tlt;
   };
