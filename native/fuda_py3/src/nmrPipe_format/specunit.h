/***/
/* specunit.h: definitions for spectral unit conversions.
/***/

/***/
/* Nota bene: the following label defines must stay in numerical order:
/***/

#define C_CM_SEC 2.99792e+10
#define C_NM_SEC 2.99792e+17

#define LAB_PTS     "Pts"
#define LAB_HZ      "Hz"
#define LAB_PPM     "ppm"
#define LAB_PCT     "%"
#define LAB_WN      "WN"
#define LAB_CM      "cm"
#define LAB_MM      "mm"
#define LAB_MIC     "mic"
#define LAB_NM      "nm"
#define LAB_WNM     "WL-nm"
#define LAB_INCH    "in"
#define LAB_C1      "c1"
#define LAB_C2      "c2"
#define LAB_PIX     "pix"
#define LAB_HEIGHT  "Hi"
#define LAB_SEC     "sec"
#define LAB_USER    "user"

#define LAB_NONE_ID   -1 

#define LAB_PTS_ID     1
#define LAB_HZ_ID      2
#define LAB_PPM_ID     3
#define LAB_PCT_ID     4
#define LAB_WN_ID      5
#define LAB_CM_ID      6
#define LAB_MM_ID      7
#define LAB_MIC_ID     8
#define LAB_NM_ID      9
#define LAB_WNM_ID    10 
#define LAB_INCH_ID   11
#define LAB_C1_ID     12
#define LAB_C2_ID     13
#define LAB_PIX_ID    14
#define LAB_HEIGHT_ID 15
#define LAB_SEC_ID    16
#define LAB_USER_ID   17

#define LAB_ID_COUNT  17

#define CLIP_NONE     -1
#define CLIP_FOLD      0
#define CLIP_TRUNC     1

#define ALL_LAB_LIST \
 LAB_PTS,  LAB_HZ, LAB_PPM, LAB_PCT, \
 LAB_WN,   LAB_CM, LAB_MM,  LAB_MIC, LAB_NM,     LAB_WNM, \
 LAB_INCH, LAB_C1, LAB_C2,  LAB_PIX, LAB_HEIGHT, LAB_SEC, \
 LAB_USER

#define SPEC_LAB_LIST LAB_PTS, LAB_HZ,  LAB_PPM, LAB_PCT, \
                      LAB_WN,  LAB_CM,  LAB_MM,  LAB_MIC, \
                      LAB_NM,  LAB_WNM, LAB_SEC, LAB_USER

#define GRAPH_LAB_LIST \
 LAB_INCH, LAB_PCT, LAB_CM, LAB_MM, LAB_C1, LAB_C2, LAB_PIX

static char *validSpecUnits[]  = {SPEC_LAB_LIST,  (char *)NULL};
static char *validGraphUnits[] = {GRAPH_LAB_LIST, (char *)NULL};
static char *allValidUnits[]   = {ALL_LAB_LIST,   (char *)NULL};

float spec2rPnt(), spec2rPnt2();
float iPnt2spec(), rPnt2spec();
float specWidth2rPnt(), rPnt2specWidth();

int   getSpecUnits(), getSpecLabel();
int   spec2iPnt(), updateOrigin(), hasSpecLabel();
int   *specList2Pnt();

int   specArgD(), wideArgD(), specFltArgD(), wideFltArgD();

int   badSpecUnits;
