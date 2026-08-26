/***/
/* Internal definitions for conversion of JEOL Delta Format NMR Data.
/* Based on their v 1.1 specifications.
/* Copyright 2006-2015 Frank Delaglio
/***/

#define JMAXDIM        8

#define JEOL_UNDEF     0
#define JEOL_TEXT      1
#define JEOL_1BIT      2
#define JEOL_2BIT      3
#define JEOL_4BIT      4
#define JEOL_6BIT      5
#define JEOL_UINT1     6
#define JEOL_UINT2     7
#define JEOL_UINT4     8
#define JEOL_INT2      9
#define JEOL_INT4     10
#define JEOL_UINT8    11
#define JEOL_FLOAT    12 
#define JEOL_DOUBLE   13 
#define JEOL_JTIME    14 
#define JEOL_JUNIT    15 
#define JEOL_EUNIT    16
#define JEOL_JVAL     17

#ifdef S_SPLINT_S
static struct NameVal jeolItemTypeList[1];
#else
static struct NameVal jeolItemTypeList[] =
   {
    "Undefined",               (float)JEOL_UNDEF,
    "Text",                    (float)JEOL_TEXT,
    "1-Bit Boolean",           (float)JEOL_1BIT,
    "2-Bit Unsigned Integer",  (float)JEOL_2BIT,
    "4-Bit Unsigned Integer",  (float)JEOL_4BIT,
    "6-Bit Unsigned Integer",  (float)JEOL_6BIT,
    "1-Byte Unsigned Integer", (float)JEOL_UINT1,
    "2-Byte Unsigned Integer", (float)JEOL_UINT2,
    "4-Byte Unsigned Integer", (float)JEOL_UINT4,
    "8-Byte Unsigned Integer", (float)JEOL_UINT8,
    "4-Byte Float",            (float)JEOL_FLOAT,
    "8-Byte Float",            (float)JEOL_DOUBLE,
    "JEOL Time Structure",     (float)JEOL_JTIME,
    "JEOL Unit Structure",     (float)JEOL_JUNIT,
    "JEOL Extended Units",     (float)JEOL_EUNIT,
    (char *)NULL,              (float)0.0
   };
#endif

#define JEOL_BIG_ENDIAN      0
#define JEOL_LITTLE_ENDIAN   1

#ifdef S_SPLINT_S
static struct NameVal jeolEndianModeList[1];
#else
static struct NameVal jeolEndianModeList[] = 
   {
    "Big-Endian",    (float)JEOL_BIG_ENDIAN,
    "Little-Endian", (float)JEOL_LITTLE_ENDIAN,
    (char *)NULL,    (float)0.0
   };
#endif

#define JEOL_FORMAT_1D         1
#define JEOL_FORMAT_2D         2
#define JEOL_FORMAT_3D         3
#define JEOL_FORMAT_4D         4
#define JEOL_FORMAT_5D         5
#define JEOL_FORMAT_6D         6
#define JEOL_FORMAT_7D         7
#define JEOL_FORMAT_8D         8
#define JEOL_FORMAT_9          9
#define JEOL_FORMAT_10        10 
#define JEOL_FORMAT_11        11 
#define JEOL_FORMAT_SMALL2D   12
#define JEOL_FORMAT_SMALL3D   13
#define JEOL_FORMAT_SMALL4D   14

#define JEOL_DATATYPE_DOUBLE   0
#define JEOL_DATATYPE_FLOAT    1

#ifdef S_SPLINT_S
static struct NameVal jeolDataTypeList[1];
static struct NameVal jeolDataFmtList[1];
#else
static struct NameVal jeolDataTypeList[] =
   {
    "64-bit Float",  (float)JEOL_DATATYPE_DOUBLE,
    "32-bit Float",  (float)JEOL_DATATYPE_FLOAT,
    (char *)NULL,    (float)0.0
   };

static struct NameVal jeolDataFmtList[] =
   {
    "1D",          (float)JEOL_FORMAT_1D,
    "2D",          (float)JEOL_FORMAT_2D,
    "3D",          (float)JEOL_FORMAT_3D,
    "4D",          (float)JEOL_FORMAT_4D,
    "5D",          (float)JEOL_FORMAT_5D,
    "6D",          (float)JEOL_FORMAT_6D,
    "7D",          (float)JEOL_FORMAT_7D,
    "8D",          (float)JEOL_FORMAT_8D,
    "Non-NMR 9",   (float)JEOL_FORMAT_8D,
    "Non-NMR 10",  (float)JEOL_FORMAT_8D,
    "Non-NMR 11",  (float)JEOL_FORMAT_8D,
    "Small 2D",    (float)JEOL_FORMAT_SMALL2D,
    "Small 3D",    (float)JEOL_FORMAT_SMALL3D,
    "Small 4D",    (float)JEOL_FORMAT_SMALL4D,
    (char *)NULL,  (float)0.0
   };
#endif

#define JEOL_SRC_NONE          0
#define JEOL_SRC_GSX           1
#define JEOL_SRC_ALPHA         2
#define JEOL_SRC_ECLIPSE       3
#define JEOL_SRC_MASS_SPEC     4
#define JEOL_SRC_COMPILER      5
#define JEOL_SRC_OTHER_NMR     6
#define JEOL_SRC_UNKNOWN       7
#define JEOL_SRC_GEMINI        8
#define JEOL_SRC_UNITY         9
#define JEOL_SRC_ASPECT       10
#define JEOL_SRC_UX           11
#define JEOL_SRC_FELIX        12
#define JEOL_SRC_LAMBDA       13
#define JEOL_SRC_GE_1280      14
#define JEOL_SRC_GE_OMEGA     15
#define JEOL_SRC_CHEMAGNETICS 16
#define JEOL_SRC_CDFF         17
#define JEOL_SRC_GALACTIC     18
#define JEOL_SRC_TRIAD        19
#define JEOL_SRC_GENERIC_NMR  20
#define JEOL_SRC_GAMMA        21
#define JEOL_SRC_JCAMP_DX     22
#define JEOL_SRC_AMX          23
#define JEOL_SRC_DMX          24
#define JEOL_SRC_ECA          25
#define JEOL_SRC_ALICE        26
#define JEOL_SRC_NMRPIPE      27
#define JEOL_SRC_SIMPSON      28

#ifdef S_SPLINT_S
static struct NameVal jeolDataSrcList[1];
#else
static struct NameVal jeolDataSrcList[] =
   {
    ".",             (float)JEOL_SRC_NONE,
    "GSX",           (float)JEOL_SRC_GSX,
    "Alpha",         (float)JEOL_SRC_ALPHA,
    "Eclipse",       (float)JEOL_SRC_ECLIPSE,
    "Mass Spec",     (float)JEOL_SRC_MASS_SPEC,
    "Compiler",      (float)JEOL_SRC_COMPILER,
    "Other NMR",     (float)JEOL_SRC_OTHER_NMR,
    "Unkown",        (float)JEOL_SRC_UNKNOWN,
    "Gemini",        (float)JEOL_SRC_GEMINI,
    "Unity",         (float)JEOL_SRC_UNITY,
    "Aspect",        (float)JEOL_SRC_ASPECT,
    "UX",            (float)JEOL_SRC_UX,
    "FELIX",         (float)JEOL_SRC_FELIX,
    "Lambda",        (float)JEOL_SRC_LAMBDA,
    "GE 1280",       (float)JEOL_SRC_GE_1280,
    "GE Omega",      (float)JEOL_SRC_GE_OMEGA,
    "Chemagnetics",  (float)JEOL_SRC_CHEMAGNETICS,
    "CDFF",          (float)JEOL_SRC_CDFF,
    "Galactic",      (float)JEOL_SRC_GALACTIC,
    "TRIAD",         (float)JEOL_SRC_TRIAD,
    "Generic NMR",   (float)JEOL_SRC_GENERIC_NMR,
    "Gamma",         (float)JEOL_SRC_GAMMA,
    "JCAMP-DX",      (float)JEOL_SRC_JCAMP_DX,
    "AMX",           (float)JEOL_SRC_AMX,
    "DMX",           (float)JEOL_SRC_DMX,
    "ECA",           (float)JEOL_SRC_ECA,
    "ALICE",         (float)JEOL_SRC_ALICE,
    "NMRPipe",       (float)JEOL_SRC_NMRPIPE,
    "SIMPSON",       (float)JEOL_SRC_SIMPSON,
    (char *)NULL,    (float)0.0
   };
#endif

#define JEOL_AXISTYPE_NONE         0
#define JEOL_AXISTYPE_REAL         1
#define JEOL_AXISTYPE_TPPI         2
#define JEOL_AXISTYPE_COMPLEX      3
#define JEOL_AXISTYPE_REAL_COMPLEX 4
#define JEOL_AXISTYPE_ENVELOPE     5

#ifdef S_SPLINT_S
static struct NameVal jeolAxisTypeList[1];
#else
static struct NameVal jeolAxisTypeList[] =
   {
    ".",           (float)JEOL_AXISTYPE_NONE,
    "Real",        (float)JEOL_AXISTYPE_REAL,
    "TPPI",        (float)JEOL_AXISTYPE_TPPI,
    "Complex",     (float)JEOL_AXISTYPE_COMPLEX,
    "Magnitude",   (float)JEOL_AXISTYPE_REAL_COMPLEX,
    "Envelope",    (float)JEOL_AXISTYPE_ENVELOPE,
    (char *)NULL,  (float)0.0
   };
#endif

#define JEOL_SCALE_YOTTA -8
#define JEOL_SCALE_ZETTA -7
#define JEOL_SCALE_EXA   -6
#define JEOL_SCALE_PETA  -5
#define JEOL_SCALE_TERA  -4
#define JEOL_SCALE_GIGA  -3
#define JEOL_SCALE_MEGA  -2
#define JEOL_SCALE_KILO  -1
#define JEOL_SCALE_NONE   0 
#define JEOL_SCALE_MILLI  1
#define JEOL_SCALE_MICRO  2
#define JEOL_SCALE_NANO   3
#define JEOL_SCALE_PICO   4
#define JEOL_SCALE_FEMTO  5
#define JEOL_SCALE_ATTO   6
#define JEOL_SCALE_ZEPTO  7

#ifdef S_SPLINT_S
static struct NameVal jeolUnitScaleList[1];
#else
static struct NameVal jeolUnitScaleList[] =
   {
    "Yotta", (float)JEOL_SCALE_YOTTA,
    "Zetta", (float)JEOL_SCALE_ZETTA,
    "Exa",   (float)JEOL_SCALE_EXA,
    "Peta",  (float)JEOL_SCALE_PETA,
    "Tera",  (float)JEOL_SCALE_TERA,
    "Giga",  (float)JEOL_SCALE_GIGA,
    "Mega",  (float)JEOL_SCALE_MEGA,
    "Kilo",  (float)JEOL_SCALE_KILO,
    "None",  (float)JEOL_SCALE_NONE,
    "Milli", (float)JEOL_SCALE_MILLI,
    "Micro", (float)JEOL_SCALE_MICRO,
    "Nano",  (float)JEOL_SCALE_NANO,
    "Pico",  (float)JEOL_SCALE_PICO,
    "Femto", (float)JEOL_SCALE_FEMTO,
    "Atto",  (float)JEOL_SCALE_ATTO,
    "Zepto", (float)JEOL_SCALE_ZEPTO,
    (char *)NULL,  (float)0.0
   };
#endif

#define JEOL_SIUNIT_NONE          0
#define JEOL_SIUNIT_ABUNDANCE     1
#define JEOL_SIUNIT_AMPERE        2
#define JEOL_SIUNIT_CANDELA       3
#define JEOL_SIUNIT_CELSIUS       4
#define JEOL_SIUNIT_COULOMB       5
#define JEOL_SIUNIT_DEGREE        6
#define JEOL_SIUNIT_ELECTRONVOLT  7
#define JEOL_SIUNIT_FARAD         8
#define JEOL_SIUNIT_SIEVERT       9
#define JEOL_SIUNIT_GRAM         10
#define JEOL_SIUNIT_GRAY         11
#define JEOL_SIUNIT_HENRY        12
#define JEOL_SIUNIT_HZ           13
#define JEOL_SIUNIT_KELVIN       14
#define JEOL_SIUNIT_JOULE        15
#define JEOL_SIUNIT_LITER        16
#define JEOL_SIUNIT_LUMEN        17
#define JEOL_SIUNIT_LUX          18
#define JEOL_SIUNIT_METER        19
#define JEOL_SIUNIT_MOLE         20
#define JEOL_SIUNIT_NEWTON       21
#define JEOL_SIUNIT_OHM          22
#define JEOL_SIUNIT_PASCAL       23
#define JEOL_SIUNIT_PERCENT      24
#define JEOL_SIUNIT_POINT        25
#define JEOL_SIUNIT_PPM          26
#define JEOL_SIUNIT_RADIAN       27
#define JEOL_SIUNIT_SECONDS      28 
#define JEOL_SIUNIT_SIEMENS      29 
#define JEOL_SIUNIT_STERADIAN    30 
#define JEOL_SIUNIT_TESLA        31
#define JEOL_SIUNIT_VOLT         32
#define JEOL_SIUNIT_WATT         33
#define JEOL_SIUNIT_WEBER        34
#define JEOL_SIUNIT_DECIBEL      35
#define JEOL_SIUNIT_DALTON       36
#define JEOL_SIUNIT_THOMPSON     37
#define JEOL_SIUNIT_GENERIC      38

#ifdef S_SPLINT_S
static struct NameVal jeolUnitTypeList[1];
#else
static struct NameVal jeolUnitTypeList[] =
   {
    "None",           (float)JEOL_SIUNIT_NONE,
    "abundance",      (float)JEOL_SIUNIT_ABUNDANCE,
    "ampere",         (float)JEOL_SIUNIT_AMPERE,
    "candela",        (float)JEOL_SIUNIT_CANDELA,
    "celsius",        (float)JEOL_SIUNIT_CELSIUS,
    "coulomb",        (float)JEOL_SIUNIT_COULOMB,
    "degree",         (float)JEOL_SIUNIT_DEGREE,
    "electronvolt",   (float)JEOL_SIUNIT_ELECTRONVOLT,
    "farad",          (float)JEOL_SIUNIT_FARAD,
    "sievert",        (float)JEOL_SIUNIT_SIEVERT,
    "gram",           (float)JEOL_SIUNIT_GRAM,
    "gray",           (float)JEOL_SIUNIT_GRAY,
    "henry",          (float)JEOL_SIUNIT_HENRY,
    "Hz",             (float)JEOL_SIUNIT_HZ,
    "kelvin",         (float)JEOL_SIUNIT_KELVIN,
    "joule",          (float)JEOL_SIUNIT_JOULE,
    "liter",          (float)JEOL_SIUNIT_LITER,
    "lumen",          (float)JEOL_SIUNIT_LUMEN,
    "lux",            (float)JEOL_SIUNIT_LUX,
    "meter",          (float)JEOL_SIUNIT_METER,
    "mole",           (float)JEOL_SIUNIT_MOLE,
    "Newton",         (float)JEOL_SIUNIT_NEWTON,
    "ohm",            (float)JEOL_SIUNIT_OHM,
    "pascal",         (float)JEOL_SIUNIT_PASCAL,
    "percent",        (float)JEOL_SIUNIT_PERCENT,
    "point",          (float)JEOL_SIUNIT_POINT,
    "PPM",            (float)JEOL_SIUNIT_PPM,
    "radian",         (float)JEOL_SIUNIT_RADIAN,
    "sec",            (float)JEOL_SIUNIT_SECONDS,
    "siemens",        (float)JEOL_SIUNIT_SIEMENS,
    "steradian",      (float)JEOL_SIUNIT_STERADIAN,
    "tesla",          (float)JEOL_SIUNIT_TESLA,
    "volt",           (float)JEOL_SIUNIT_VOLT,
    "watt",           (float)JEOL_SIUNIT_WATT,
    "weber",          (float)JEOL_SIUNIT_WEBER,
    "decibel",        (float)JEOL_SIUNIT_DECIBEL,
    "dalton",         (float)JEOL_SIUNIT_DALTON,
    "thompson",       (float)JEOL_SIUNIT_THOMPSON,
    "generic",        (float)JEOL_SIUNIT_GENERIC,
    (char *)NULL,  (float)0.0
   };
#endif

#define JEOL_RULER_RANGED     0
#define JEOL_RULER_LISTED_OLD 1
#define JEOL_RULER_SPARSE     2
#define JEOL_RULER_LISTED     3

#ifdef S_SPLINT_S
static struct NameVal jeolRulerTypeList[1];
#else
static struct NameVal jeolRulerTypeList[] =
   {
    "Ranged",             (float)JEOL_RULER_RANGED,
    "Listed (Old Type)",  (float)JEOL_RULER_LISTED_OLD,
    "Sparse",             (float)JEOL_RULER_SPARSE,
    "Listed",             (float)JEOL_RULER_LISTED,
    (char *)NULL,         (float)0.0
   };
#endif

#define JEOL_DAYFRAC (86400.0/65535.0)

#ifdef S_SPLINT_S
static struct NameVal monthList[1];
#else
static struct NameVal monthList[] = 
   {
    "Jan", (float)1,
    "Feb", (float)2,
    "Mar", (float)3,
    "Apr", (float)4,
    "May", (float)5,
    "Jun", (float)6,
    "Jul", (float)7,
    "Aug", (float)8,
    "Sep", (float)9,
    "Oct", (float)10,
    "Nov", (float)11,
    "Dec", (float)12,
    (char *)NULL, (float)0.0
   };
#endif

struct jTime
   { 
    int year;
    int month;
    int day;
    int dayFrac;
    int hour;
    int min;
    int sec;
   };

#define JEOL_NUM_EXTUNITS 5

struct jUnit
   {
    int unitType;
    int unitExp;
    int scaleType;
   };

struct jExtendedUnit
   {
    struct jUnit jUnitInfo[JEOL_NUM_EXTUNITS];
    int unitExp;
   };

#define JEOL_PARMVAL_NONE -1 /* Not JEOL's definition. */
#define JEOL_PARMVAL_STR   0
#define JEOL_PARMVAL_INT   1
#define JEOL_PARMVAL_FLT   2
#define JEOL_PARMVAL_Z     3
#define JEOL_PARMVAL_INF   4

#define JEOL_INF_NEG       1
#define JEOL_INF_MINUS1    2
#define JEOL_INF_ZERO      3
#define JEOL_INF_PLUS1     4
#define JEOL_INF_POS       5

#define JEOL_JVAL_STRLEN   16

#ifdef S_SPLINT_S
struct NameVal jeolParmValList[1];
struct NameVal jeolInfList[1];
#else
struct NameVal jeolParmValList[] =
   {
    "None",       (float)JEOL_PARMVAL_NONE,
    "String",     (float)JEOL_PARMVAL_STR,
    "Integer",    (float)JEOL_PARMVAL_INT,
    "Float",      (float)JEOL_PARMVAL_FLT,
    "Complex",    (float)JEOL_PARMVAL_Z,
    "Infinity",   (float)JEOL_PARMVAL_INF,
    (char *)NULL, (float)0.0
   };

struct NameVal jeolInfList[] = 
   {
    "Negative Infinity",  (float)JEOL_INF_NEG,
    "Negative 1",         (float)JEOL_INF_MINUS1,
    "Zero",               (float)JEOL_INF_ZERO,
    "Positive 1",         (float)JEOL_INF_PLUS1,
    "Positive Infinity",  (float)JEOL_INF_POS,
    (char *)NULL,         (float)0.0
   };
#endif

struct jComplex
   {
    double re;
    double im;
   };

struct jVal
   {
    char   s[NAMELEN+1];
    double d;
    int    i;
    int    inf;
    struct jComplex z;
   };

struct DeltaHdrInfo
   {
    char   fileID[NAMELEN+1];
    int    endianMode;
    int    majorVersion;
    int    minorVersion;
    int    dimCount;
    int    dimExists[JMAXDIM];
    int    dataType;
    int    dataFormat;
    int    instrument;
    int    translate[JMAXDIM];
    int    axisType[JMAXDIM];
    struct jUnit unitList[JMAXDIM];
    char   title[NAMELEN+1];
    int    axisRanged[JMAXDIM];
    int    sizeList[JMAXDIM];
    int    offsetStart[JMAXDIM];
    int    offsetStop[JMAXDIM];
    double axisStart[JMAXDIM];
    double axisStop[JMAXDIM];
    struct jTime creationTime;
    struct jTime revisionTime;
    char   nodeName[NAMELEN+1];
    char   site[NAMELEN+1];
    char   author[NAMELEN+1];
    char   comment[NAMELEN+1];
    char   axisTitles[JMAXDIM*(NAMELEN+1)];
    double baseFreq[JMAXDIM];
    double zeroPoint[JMAXDIM];
    int    reversed[JMAXDIM];
    char   reserved1[NAMELEN+1];
    int    annotationFlag;
    int    reserved2[7];
    int    historyUsed;
    int    historyLength;
    int    paramStart;
    int    paramLength;
    int    listStart[JMAXDIM];
    int    listLength[JMAXDIM];
    int    dataStart;
    int    dataLength;
    int    contextStart;
    int    contextLength;
    int    annoteStart;
    int    annoteLength;
    int    totalSize;
    int    unitLoc[JMAXDIM];
    struct jExtendedUnit eUnits[2];
   };

struct DeltaParmHdr
   {
    int parmSize;
    int loID;
    int hiID;
    int totalSize;
   };

struct DeltaParmVal
   {
    char   class[NAMELEN+1];
    int    unitScale;
    struct jUnit units[2];
    struct jVal  val;
    int    valType;
    char   name[NAMELEN+1];
   };

struct DeltaHdrItem
  {
   char   *dest;
   char   *itemName;
   int    byteOffset;
   int    bitOffset;
   int    byteCount;
   int    itemCount; 
   int    itemType;
   struct NameVal *nameList;
  };

static struct DeltaHdrInfo  deltaHdr;
static struct DeltaParmHdr  deltaParmHdr;
static struct DeltaParmVal  deltaParmVal;
static struct jVal          deltaValVal;

#ifdef S_SPLINT_S
static struct DeltaHdrItem deltaHdrList1_1[1];
static struct DeltaHdrItem deltaParmHdrList1_1[1];
static struct DeltaHdrItem deltaParmValList1_1[1];
#else
static struct DeltaHdrItem deltaHdrList1_1[] = 
   {
    (char *)deltaHdr.fileID,          "FileIdentifier",         0, 0,    8, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)&deltaHdr.endianMode,     "Endian",                 8, 0,    1, 1,       JEOL_UINT1,  jeolEndianModeList,
    (char *)&deltaHdr.majorVersion,   "MajorVersion",           9, 0,    1, 1,       JEOL_UINT1,  (struct NameVal *)NULL,
    (char *)&deltaHdr.minorVersion,   "MinorVersion",          10, 0,    2, 1,       JEOL_UINT2,  (struct NameVal *)NULL,
    (char *)&deltaHdr.dimCount,       "DataDimensionNumber",   12, 0,    1, 1,       JEOL_UINT1,  (struct NameVal *)NULL,
    (char *)deltaHdr.dimExists,       "DataDimensionExist",    13, 0,    1, JMAXDIM, JEOL_1BIT,   (struct NameVal *)NULL,
    (char *)&deltaHdr.dataType,       "DataType",              14, 0,    1, 1,       JEOL_2BIT,   jeolDataTypeList,
    (char *)&deltaHdr.dataFormat,     "DataFormat",            14, 2,    0, 1,       JEOL_6BIT,   jeolDataFmtList,
    (char *)&deltaHdr.instrument,     "Instrument",            15, 0,    1, 1,       JEOL_UINT1,  jeolDataSrcList,
    (char *)deltaHdr.translate,       "Translate",             16, 0,    8, JMAXDIM, JEOL_UINT1,  (struct NameVal *)NULL,
    (char *)deltaHdr.axisType,        "DataAxisType",          24, 0,    8, JMAXDIM, JEOL_UINT1,  jeolAxisTypeList, 
    (char *)deltaHdr.unitList,        "DataUnits",             32, 0,   16, JMAXDIM, JEOL_JUNIT,  (struct NameVal *)NULL,
    (char *)deltaHdr.title,           "Title",                 48, 0,  124, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.axisRanged,      "DataAxisRanged",       172, 0,    4, JMAXDIM, JEOL_4BIT,   jeolRulerTypeList,
    (char *)deltaHdr.sizeList,        "DataPoints",           176, 0,   32, JMAXDIM, JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)deltaHdr.offsetStart,     "DataOffsetStart",      208, 0,   32, JMAXDIM, JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)deltaHdr.offsetStop,      "DataOffsetStop",       240, 0,   32, JMAXDIM, JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)deltaHdr.axisStart,       "DataAxisStart",        272, 0,   64, JMAXDIM, JEOL_DOUBLE, (struct NameVal *)NULL,
    (char *)deltaHdr.axisStop,        "DataAxisStop",         336, 0,   64, JMAXDIM, JEOL_DOUBLE, (struct NameVal *)NULL,
    (char *)&deltaHdr.creationTime,   "CreationTime",         400, 0,    4, 1,       JEOL_JTIME,  (struct NameVal *)NULL,
    (char *)&deltaHdr.revisionTime,   "RevisionTime",         404, 0,    4, 1,       JEOL_JTIME,  (struct NameVal *)NULL,
    (char *)deltaHdr.nodeName,        "NodeName",             408, 0,   16, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.site,            "Site",                 424, 0,  128, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.author,          "Author",               552, 0,  128, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.comment,         "Comment",              680, 0,  128, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.axisTitles,      "DataAxisTitles",       808, 0,  256, JMAXDIM, JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)deltaHdr.baseFreq,        "BaseFreq",            1064, 0,   64, JMAXDIM, JEOL_DOUBLE, (struct NameVal *)NULL,
    (char *)deltaHdr.zeroPoint,       "ZeroPoint",           1128, 0,   64, JMAXDIM, JEOL_DOUBLE, (struct NameVal *)NULL,
    (char *)deltaHdr.reversed,        "Reversed",            1192, 0,    8, JMAXDIM, JEOL_1BIT,   (struct NameVal *)NULL,
    (char *)deltaHdr.reserved1,       "Reserved1",           1200, 0,    3, 1,       JEOL_UNDEF,  (struct NameVal *)NULL,
    (char *)&deltaHdr.annotationFlag, "AnnotationValid",     1203, 0,    1, 1,       JEOL_1BIT,   (struct NameVal *)NULL,
    (char *)deltaHdr.reserved2,       "Reserved2",           1203, 1,    0, 7,       JEOL_1BIT,   (struct NameVal *)NULL,
    (char *)&deltaHdr.historyUsed,    "HistoryUsed",         1204, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.historyLength,  "HistoryLength",       1208, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.paramStart,     "ParamStart",          1212, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.paramLength,    "ParamLength",         1216, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)deltaHdr.listStart,      "ListStart",            1220, 0,   32, JMAXDIM, JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)deltaHdr.listLength,     "ListLength",           1252, 0,   32, JMAXDIM, JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.dataStart,      "DataStart",           1284, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.dataLength,     "DataLength",          1288, 0,    8, 1,       JEOL_UINT8,  (struct NameVal *)NULL,
    (char *)&deltaHdr.contextStart,   "ContextStart",        1296, 0,    8, 1,       JEOL_UINT8,  (struct NameVal *)NULL,
    (char *)&deltaHdr.contextLength,  "ContextLength",       1304, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.annoteStart,    "AnnoteStart",         1308, 0,    8, 1,       JEOL_UINT8,  (struct NameVal *)NULL,
    (char *)&deltaHdr.annoteLength,   "AnnoteLength",        1316, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaHdr.totalSize,      "TotalSize",           1320, 0,    8, 1,       JEOL_UINT8,  (struct NameVal *)NULL,
    (char *)deltaHdr.unitLoc,         "UnitLocation",        1328, 0,    8, JMAXDIM, JEOL_UINT1,  (struct NameVal *)NULL,
    (char *)deltaHdr.eUnits,          "ExtendedUnits",       1336, 0,   24, 2,       JEOL_EUNIT,  (struct NameVal *)NULL,
    (char *)NULL,                     (char *)NULL,            -1, 0,    0, 0,       0,           (struct NameVal *)NULL,
   };

static struct DeltaHdrItem deltaParmHdrList1_1[] = 
   {
    (char *)&deltaParmHdr.parmSize,   "ParamSize",              0, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaParmHdr.loID,       "LowIndex",               4, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL, 
    (char *)&deltaParmHdr.hiID,       "HighIndex",              8, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)&deltaParmHdr.totalSize,  "ParamTotalSize",        12, 0,    4, 1,       JEOL_UINT4,  (struct NameVal *)NULL,
    (char *)NULL,                     (char *)NULL,            -1, 0,    0, 0,       0,           (struct NameVal *)NULL,
   };

/* Parameter value is interpreted last, so its type will be known: */

static struct DeltaHdrItem deltaParmValList1_1[] =
   {
    (char *)deltaParmVal.name,        "Name",             36, 0,   28, 1,       JEOL_TEXT,   (struct NameVal *)NULL,
    (char *)&deltaParmVal.valType,    "ValueType",        32, 0,    4, 1,       JEOL_UINT4,  jeolParmValList,
    (char *)&deltaParmVal.val,        "Value",            16, 0,   16, 1,       JEOL_JVAL,   (struct NameVal *)NULL,
    (char *)deltaParmVal.class,       "Class",             0, 0,    4, 1,       JEOL_UNDEF,  (struct NameVal *)NULL,
    (char *)deltaParmVal.units,       "Units",             6, 0,   10, 5,       JEOL_JUNIT,  (struct NameVal *)NULL,
    (char *)&deltaParmVal.unitScale,  "UnitScale",         4, 0,    2, 1,       JEOL_INT2,   (struct NameVal *)NULL,
    (char *)NULL,                     (char *)NULL,       -1, 0,    0, 0,       0,           (struct NameVal *)NULL,
   };
#endif

static int    uc2int2();
static int    uc2int4();
static int    uc2uint2();
static int    uc2uint4();
static int    uc2uint8();
static float  uc2flt4();
static double uc2flt8();

static int    trimString();
static int    getSMXSize();
static int    openNMR();
static int    closeNMR();

static int    workAlloc();
static int    workFree();

static int    initDeltaVars();
static int    getDeltaArgs();
static int    getDeltaWordSize();
static int    getDeltaAQ2D();

static int    isDeltaQuad();
static int    isDeltaTimeDomain();
static int    isDeltaPPM();
static int    isDeltaHz();

static float  getDeltaSW();
static float  getDeltaObs();
static float  getDeltaCar();
static float  getDeltaOrig();
static float  getDeltaTemperature();

static int    getDeltaDFSize();
static float  getDeltaDF();

static double getDeltaUnitVal();

static char   *space2us();

