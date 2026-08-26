/******************************************************************************/
/*                                                                            */
/*                   ---- NIH NMR Software System ----                        */
/*                         Copyright 1992 - 1995                              */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "prec.h"

#define FMT "%m %e %Y %k %M %S"

/***/
/* getTime: returns date and time; NULL pointers are not changed.
/***/

int getTime( month, day, year, hour, min, sec )

   int *month, *day, *year, *hour, *min, *sec;
{
    time_t iclock;
    struct tm *currentTime;
    char   ctemp[NAMELEN+1];
    int    imonth, iday, iyear, ihour, imin, isec;

#ifdef Y2KTEST
    (void) strcpy( ctemp, "11 17 1998 10 00 00" );
#else
    (void) time( &iclock );

    currentTime = localtime( &iclock );

    (void) strftime( ctemp, NAMELEN, FMT, currentTime );
#endif

    (void) sscanf( ctemp,
                   "%d %d %d %d %d %d",
                   &imonth, &iday, &iyear, &ihour, &imin, &isec );

    if (month) *month = imonth;
    if (day)   *day   = iday;
    if (year)  *year  = iyear;
    if (hour)  *hour  = ihour;
    if (min)   *min   = imin;
    if (sec)   *sec   = isec;

    return( 0 );
}
