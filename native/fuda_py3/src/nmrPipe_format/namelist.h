/***/
/* namelist: various lists mapping keywords to values.
/***/

#include "namelist2.h"
#include "dimloc.h"

#define KEYWORD_BAD -1  /* The keyword was not found in the given list.   */
#define KEYWORD_OK   0  /* The keyword was found successfuly.             */
#define KEYWORD_NUM  1  /* A number was specified, rather than a keyword. */
#define KEYWORD_NONE 2  /* No keyword was specified.                      */

#define WHEN_NEVER     0
#define WHEN_SOMETIMES 1
#define WHEN_ALWAYS    2

#ifdef WHEN_LIST
#ifdef S_SPLINT_S
   static struct NameVal whenList[1];
#else
   static struct NameVal whenList[] =
      {
       "Never",           (float)WHEN_NEVER,
       "Sometimes",       (float)WHEN_SOMETIMES,
       "Always",          (float)WHEN_ALWAYS,
       (char *)NULL,      (float)0.0
      };
#endif
#endif
 
#ifdef MODELIST
#ifdef S_SPLINT_S
   static struct NameVal modeList[1];
#else
   static struct NameVal modeList[] =
      {
       "States",          0.0,
       "States-TPPI",     0.0,
       "Complex",         0.0,
       "Sim",             0.0,
       "DQD",             0.0,
       "Simultaneous",    0.0,
       "Real",            1.0,
       "TPPI",            1.0,
       "Seq",             2.0,
       "Sequential",      2.0,
       "Bruk",            2.0,
       "Bruker",          2.0,
       "Rance",           3.0,
       "SE",              3.0,
       "Grad",            4.0,
       "Gradient",        4.0,
       "Rance-Kay",       4.0,
       "Echo-AntiEcho",   4.0,
       "States-N",        0.0,
       "States-TPPI-N",   0.0,
       "Complex-N",       0.0,
       "Sim-N",           0.0,
       "DQD-N",           0.0,
       "Simultaneous-N",  0.0,
       "Real-N",          1.0,
       "TPPI-N",          1.0,
       "Seq-N",           2.0,
       "Sequential-N",    2.0,
       "Bruk-N",          2.0,
       "Bruker-N",        2.0,
       "Rance",           3.0,
       "SE-N",            3.0,
       "Grad-N",          4.0,
       "Gradient-N",      4.0,
       "Rance-Kay-N",     4.0,
       "Echo-AntiEcho-N", 4.0,
       (char *)NULL,      0.0
      };
#endif
#endif

#ifdef AQLIST
#ifdef S_SPLINT_S
   static struct NameVal aq2DList[1];
#else
   static struct NameVal aq2DList[] =
      {
       "Real",        0.0,
       "Magnitude",   0.0,
       "Real",        1.0,
       "TPPI",        1.0,
       "Complex",     2.0,
       "States",      2.0,
       "States-TPPI", 2.0,
       "Image",       3.0,
       "Array",       4.0,
       (char *)NULL,  0.0
      };
#endif
#endif

#ifdef FTLIST
#ifdef S_SPLINT_S
   static struct NameVal ftList[1];
#else
   static struct NameVal ftList[] =
      {
       "Time",        0.0,
       "Array",       0.0,
       "Freq",        1.0,
       "Frequency",   1.0,
       (char *)NULL,  0.0
      };
#endif
#endif

#ifdef FTDOMLIST
#ifdef S_SPLINT_S
   static struct NameVal ftDomList[1];
#else
   static struct NameVal ftDomList[] =
      {
       "Freq",        0.0,
       "Frequency",   0.0,
       "Spectra",     0.0,
       "Spectral",    0.0,
       "Space",       1.0,
       "Spatial",     1.0,
       "Dist",        1.0,
       "Distance",    1.0,
       (char *)NULL,  0.0
      };
#endif
#endif

#ifdef ACQMETHLIST
#ifdef S_SPLINT_S
   static struct NameVal acqMethList[1];
#else
   static struct NameVal acqMethList[] =
      {
       "FT",          0.0,
       "Fourier",     0.0,
       "Dir",         1.0,
       "Direct",      1.0,
       (char *)NULL,  0.0
      };
#endif
#endif

#ifdef SASE_LIST
#ifdef S_SPLINT_S
   static struct NameVal saseList[1];
#else
   static struct NameVal saseList[] =
      {
       "Exp",         SASE_SIG_EXP,
       "Gauss",       SASE_SIG_GAUSS,
       "CT",          SASE_SIG_CT,
       (char *)NULL,  0.0
      };
#endif
#endif

#ifdef DSPLIST
#ifdef S_SPLINT_S
   static struct NameVal dspList[1];
#else
   static struct NameVal dspList[] =
      {
       "States",         AQ_MOD_QSIM,
       "States-TPPI",    AQ_MOD_QSIM,
       "Complex",        AQ_MOD_QSIM,
       "Sim",            AQ_MOD_QSIM,
       "DQD",            AQ_MOD_DQD,
       "Simultaneous",   AQ_MOD_QSIM,
       "Real",           AQ_MOD_QF,
       "TPPI",           AQ_MOD_QF,
       "Seq",            AQ_MOD_QSEQ,
       "Sequential",     AQ_MOD_QSEQ,
       "Bruk",           AQ_MOD_QSEQ,
       "Bruker",         AQ_MOD_QSEQ,
       "States-N",       AQ_MOD_QSIM,
       "States-TPPI-N",  AQ_MOD_QSIM,
       "Complex-N",      AQ_MOD_QSIM,
       "Sim-N",          AQ_MOD_QSIM,
       "DQD-N",          AQ_MOD_DQD,
       "Simultaneous-N", AQ_MOD_QSIM,
       "Real-N",         AQ_MOD_QF,
       "TPPI-N",         AQ_MOD_QF,
       "Seq-N",          AQ_MOD_QSEQ,
       "Sequential-N",   AQ_MOD_QSEQ,
       "Bruk-N",         AQ_MOD_QSEQ,
       "Bruker-N",       AQ_MOD_QSEQ,
       (char *)NULL,     0.0
      };
#endif
#endif

#ifdef ALTLIST 
#ifdef S_SPLINT_S
   static struct NameVal altList[1];
#else
   static struct NameVal altList[] =
      {
       "States",            0.0,
       "States-TPPI",       2.0,
       "Complex",           0.0,
       "Sim",               0.0,
       "DQD",               0.0,
       "Simultaneous",      0.0,
       "Real",              0.0,
       "TPPI",              0.0,
       "Seq",               1.0,
       "Sequential",        1.0,
       "Bruk",              1.0,
       "Bruker",            1.0,
       "Rance",             0.0,
       "SE",                0.0,
       "Grad",              0.0,
       "Gradient",          0.0,
       "Rance-Kay",         0.0,
       "Echo-AntiEcho",     0.0,
       "States-N",         16.0,
       "States-TPPI-N",    18.0,
       "Complex-N",        16.0,
       "Sim-N",            16.0,
       "DQD-N",            16.0,
       "Simultaneous-N",   16.0,
       "Real-N",           16.0,
       "TPPI-N",           16.0,
       "Seq-N",            17.0,
       "Sequential-N",     17.0,
       "Bruk-N",           17.0,
       "Bruker-N",         17.0,
       "Rance-N",          16.0,
       "SE-N",             16.0,
       "Grad-N",           16.0,
       "Gradient-N",       16.0,
       "Rance-Kay-N",      16.0,
       "Echo-AntiEcho-N",  16.0,
       (char *)NULL,        0.0
      };
#endif
#endif

#ifdef UNITS_LIST
#ifdef S_SPLINT_S
   static struct NameVal unitLabelList[1];
#else
   static struct NameVal unitLabelList[] =
      {
       LAB_PTS,    LAB_PTS_ID,
       LAB_HZ,     LAB_HZ_ID,
       LAB_PPM,    LAB_PPM_ID,
       LAB_PCT,    LAB_PCT_ID,
       LAB_WN,     LAB_WN_ID,
       LAB_CM,     LAB_CM_ID,
       LAB_MM,     LAB_MM_ID,
       LAB_MIC,    LAB_MIC_ID,
       LAB_NM,     LAB_NM_ID,
       LAB_WNM,    LAB_WNM_ID,
       LAB_INCH,   LAB_INCH_ID,
       LAB_C1,     LAB_C1_ID,
       LAB_C2,     LAB_C2_ID,
       LAB_PIX,    LAB_PIX_ID,
       LAB_HEIGHT, LAB_HEIGHT_ID,
       LAB_USER,   LAB_USER_ID,
       (char *)NULL,        0.0
      };
#endif
#endif

#ifdef APLIST
#ifdef S_SPLINT_S
   static struct NameVal apList[1];
#else
   static struct NameVal apList[] =
      {
       "None",        (float)APOD_NULL,
       "SP",          (float)APOD_SP,
       "EM",          (float)APOD_EM,
       "GM",          (float)APOD_GM,
       "TM",          (float)APOD_TM,
       "ZE",          (float)APOD_ZE,
       "TRI",         (float)APOD_TRI,
       "GMB",         (float)APOD_GMB,
       "JMOD",        (float)APOD_JMOD,
       "TMC",         (float)APOD_TMC,
       "FILE",        (float)APOD_FILE,
       (char *)NULL,  (float)0.0
      };
#endif
#endif

#ifdef SZ_LIST
#ifdef S_SPLINT_S
   static struct NameVal szList[1];
#else
   static struct NameVal szList[] =
      {
       "Small",        0.0,
       "Medium",       1.0,
       "Large",        2.0,
       "ExtraLarge",   3.0,
       (char *)NULL,   0.0
      };
#endif
#endif

#ifdef GXLIST
#ifdef S_SPLINT_S
   static struct NameVal gxList[1];
#else
   static struct NameVal gxList[] =
      {
       "clear",         (float)GXclear,
       "and",           (float)GXand,
       "andReverse",    (float)GXandReverse,
       "copy",          (float)GXcopy,
       "andInverted",   (float)GXandInverted,
       "noop",          (float)GXnoop,
       "xor",           (float)GXxor,
       "or",            (float)GXor,
       "nor",           (float)GXnor,
       "equiv",         (float)GXequiv,
       "invert",        (float)GXinvert,
       "orReverse",     (float)GXorReverse,
       "copyInverted",  (float)GXcopyInverted,
       "orInverted",    (float)GXorInverted,
       "nand",          (float)GXnand,
       "set",           (float)GXset,
       (char *)NULL,    0.0
      };
#endif
#endif

#ifdef VIS_LIST

#define DEFAULT_VISUAL_CLASS -666

#ifdef S_SPLINT_S
   static struct NameVal visList[1];
#else
   static struct NameVal visList[] =
      {
       "StaticGray",    (float)StaticGray,
       "GrayScale",     (float)GrayScale,
       "StaticColor",   (float)StaticColor,
       "PseudoColor",   (float)PseudoColor,
       "TrueColor",     (float)TrueColor,
       "DirectColor",   (float)DirectColor,
       "Default",       (float)DEFAULT_VISUAL_CLASS,
       (char *)NULL,    0.0
      };
#endif
#endif

#ifdef DIM_LIST
#ifdef S_SPLINT_S
   static struct NameVal dimNameList[1];
#else
   static struct NameVal dimNameList[] =
      {
       "BAD_DIM",       (float)BAD_DIM,
       "NULL_DIM",      (float)NULL_DIM,
       "NULL",          (float)NULL_DIM,
       "NONE",          (float)NULL_DIM,
       "CUR_XDIM",      (float)CUR_XDIM,
       "CUR_YDIM",      (float)CUR_YDIM,
       "CUR_ZDIM",      (float)CUR_ZDIM,
       "CUR_ADIM",      (float)CUR_ADIM,
       "X_AXIS",        (float)CUR_XDIM,
       "Y_AXIS",        (float)CUR_YDIM,
       "Z_AXIS",        (float)CUR_ZDIM,
       "A_AXIS",        (float)CUR_ADIM,
       "X",             (float)CUR_XDIM,
       "Y",             (float)CUR_YDIM,
       "Z",             (float)CUR_ZDIM,
       "A",             (float)CUR_ADIM,
       "ABS_XDIM",      (float)ABS_XDIM,
       "ABS_YDIM",      (float)ABS_YDIM,
       "ABS_ZDIM",      (float)ABS_ZDIM,
       "ABS_ADIM",      (float)ABS_ADIM,
       "CUR_HDIM",      (float)CUR_HDIM,
       "CUR_VDIM",      (float)CUR_VDIM,
       (char *)NULL,    0.0
      };
#endif
#endif

#ifdef PTR_CODE_LIST
#ifdef S_SPLINT_S
   static struct NameVal ptrPrefixList[1];
   static struct NameVal ptrNameList[1];
   static struct NameVal ptrSizeList[1];
#else
   static struct NameVal ptrPrefixList[] =
      {
       "ptr",          (float)PTR_PTR,
       "chr",          (float)PTR_CHAR,
       "byt",          (float)PTR_BYTE,
       "int",          (float)PTR_INT,
       "shr",          (float)PTR_SHORT,
       "lng",          (float)PTR_LONG,
       "flt",          (float)PTR_FLT,
       "dbl",          (float)PTR_DBL,
       "uch",          (float)PTR_UCHAR,
       "uin",          (float)PTR_UINT,
       "ush",          (float)PTR_USHORT,
       "uln",          (float)PTR_ULONG,
       (char *)NULL,   0.0
      };

   static struct NameVal ptrNameList[] =
      {
       "ptr",          (float)PTR_PTR,
       "char",         (float)PTR_CHAR,
       "byte",         (float)PTR_BYTE,
       "int",          (float)PTR_INT,
       "short",        (float)PTR_SHORT,
       "long",         (float)PTR_LONG,
       "float",        (float)PTR_FLT,
       "double",       (float)PTR_DBL,
       "uchar",        (float)PTR_UCHAR,
       "uint",         (float)PTR_UINT,
       "ushort",       (float)PTR_USHORT,
       "ulong",        (float)PTR_ULONG,
       (char *)NULL,   0.0
      };

   static struct NameVal ptrSizeList[] =
      {
       "ptr",          (float)sizeof(char *),
       "chr",          (float)sizeof(char),
       "byt",          (float)sizeof(unsigned char),
       "int",          (float)sizeof(int),
       "shr",          (float)sizeof(short),
       "lng",          (float)sizeof(long),
       "flt",          (float)sizeof(float),
       "dbl",          (float)sizeof(double),
       "uch",          (float)sizeof(unsigned char),
       "uin",          (float)sizeof(unsigned int),
       "ush",          (float)sizeof(unsigned short),
       "uln",          (float)sizeof(unsigned long),
       (char *)NULL,   0.0
      };
#endif
#endif

#ifdef GDB_LIST
#ifdef S_SPLINT_S
   static struct NameVal gdbObjTypeList[1];
   static struct NameVal gdbVarTypeList[1];
#else
   static struct NameVal gdbObjTypeList[] =
      {
       "null",          (float)GDB_NULL,
       "object",        (float)GDB_OBJECT,
       "data",          (float)GDB_DATA,
       "cell",          (float)GDB_CELL,
       "entry",         (float)GDB_ENTRY,
       "variable",      (float)GDB_VARIABLE,
       "text",          (float)GDB_TEXT,
       "table",         (float)GDB_TABLE,
       "dbase",         (float)GDB_DBASE,
       "ptr",           (float)GDB_PTR,
       "root",          (float)GDB_ROOT,
       (char *)NULL,    0.0
      };

   static struct NameVal gdbVarTypeList[] =
      {
       "null",           (float)GDB_ITEM_NULL,
       "text",           (float)GDB_ITEM_TXT,
       "int",            (float)GDB_ITEM_INT,
       "float",          (float)GDB_ITEM_FLT,
       "double",         (float)GDB_ITEM_DBL,
       "ptr",            (float)GDB_ITEM_PTR,
       (char *)NULL,     0.0
      };
#endif
#endif

char *getNameStrByVal();

int  strStrArgD();
int  fltStrArgD();
int  getNameByVal();
int  getValByName();
int  intStrArgD();


