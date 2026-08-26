/***/
/* version.h: some convenient definitions for version display.
/***/

#include "prec.h"

#define NMRPIPE_REV_TXT "2020.219.15.07"
#define NMRPIPE_VER_TXT "Version 10.9"

#define FPR_VERSION (void)fprintf( stderr, " ** NMRPipe System %s Rev %s %d-bit **\n", NMRPIPE_VER_TXT, NMRPIPE_REV_TXT, (int)(8*sizeof(void *)) )
