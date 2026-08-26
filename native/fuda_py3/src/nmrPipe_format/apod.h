/***/
/* Definitions for apodization functions:
/***/

int apodInit(), apodInit2(), apodFree(), apodFn();
int apodEM(), apodGM(), apodGMB(), apodTM(), apodSB(), apodTRI();
int apodSP(), apodJMOD(), apodTMC(), apodTRI(), apodZE();

/***/
/* Add to the list, but do not change existing values:
/***/

#define APOD_NULL  0
#define APOD_SP    1
#define APOD_EM    2
#define APOD_GM    3
#define APOD_TM    4
#define APOD_ZE    5
#define APOD_TRI   6
#define APOD_GMB   7
#define APOD_JMOD  8
#define APOD_TMC   9
#define APOD_FILE  101

#define APOD_COUNT 11 

