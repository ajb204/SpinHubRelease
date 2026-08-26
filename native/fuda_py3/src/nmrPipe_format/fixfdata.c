/******************************************************************************/
/*                                                                            */
/*                   ---- NIH NMR Software System ----                        */
/*                        Copyright 1992 and 1993                             */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/


/***/
/* fixfdata: rough adjustment of older file headers for forward-compatiblity.
/***/

#include <stdio.h>
#include "fdatap.h"

int fixfdata( fdata )

   float fdata[FDATASIZE];
{
    int veryOldFlag, i, count;

/***/
/* Abort on null header:
/***/

    veryOldFlag = 0;
    count       = 0;

    for( i = 0; i < FDATASIZE; i++ )
       {
        if (fdata[i] == 0.0) count++;
       }

    if (count == FDATASIZE) return( 1 );

/***/
/* Fix null format and byte order.
/* Fix null SIZES and their QUADFLAGS.
/* Fix missing DIMCOUNT.
/* Fix missing DIMORDER.
/* Fix missing FILECOUNT.
/*
/* For very old data, fix all QUADFLAGS and LABELS
/***/

    fdata[FD2DVIRGIN] = 1.0;

    if ((int)fdata[FDFLTFORMAT] == 0) fdata[FDFLTFORMAT] = FDFMTCONS;
    if ((int)fdata[FDFLTORDER]  == 0) fdata[FDFLTORDER]  = FDORDERCONS;

    if ((int)fdata[FDSPECNUM] < 1)
       {
        fdata[FDSPECNUM]    = 1;
        fdata[FDF1QUADFLAG] = 1;
       }

    if ((int)fdata[FDF3SIZE] < 1)
       {
        fdata[FDF3SIZE]     = 1;
        fdata[FDF3QUADFLAG] = 1;
       }

    if ((int)fdata[FDF4SIZE] < 1)
       {
        fdata[FDF4SIZE]     = 1;
        fdata[FDF4QUADFLAG] = 1;
       }

    if ((int)fdata[FDREALSIZE] < 1 && (int)fdata[FDF2APOD] > 0)
       {
        fdata[FDREALSIZE] = fdata[FDF2APOD];
       }

    if ((int)fdata[FDF2APOD] < 1 && (int)fdata[FDREALSIZE] > 0)
       {
        fdata[FDF2APOD] = fdata[FDREALSIZE];
       }

    if ((int)fdata[FDF2APOD] < 1)
       {
        fdata[FDREALSIZE] = fdata[FDSIZE];
        fdata[FDF2APOD]   = fdata[FDSIZE];
       }

    if ((int)fdata[FDDIMCOUNT] == 0)
       {
        veryOldFlag = 1;

        if ((int)fdata[FDSPECNUM] > 1)
           fdata[FDDIMCOUNT] = 2;
        else
           fdata[FDDIMCOUNT] = 1;
       }

    if ((int)fdata[FDDIMORDER1] == 0 || (int)fdata[FDDIMORDER2] == 0)
       {
        if ((int)fdata[FDTRANSPOSED] == 0)
           {
            fdata[FDDIMORDER1] = 2;
            fdata[FDDIMORDER2] = 1;
            fdata[FDDIMORDER3] = 3;
            fdata[FDDIMORDER4] = 4;
           }
        else
           {
            fdata[FDDIMORDER1] = 1;
            fdata[FDDIMORDER2] = 2;
            fdata[FDDIMORDER3] = 3;
            fdata[FDDIMORDER4] = 4;
           }   
       }

    if ((int)fdata[FDDIMORDER3] == 0) fdata[FDDIMORDER3] = 3;
    if ((int)fdata[FDDIMORDER4] == 0) fdata[FDDIMORDER4] = 4;

    if ((int)fdata[FDFILECOUNT] == 0)
       {
        fdata[FDFILECOUNT] = fdata[FDF3SIZE]*fdata[FDF4SIZE];
       }

    if (veryOldFlag)
       {
        if ((int)fdata[FDQUADFLAG] == 1)
           {
            fdata[FDF2QUADFLAG] = 1;
            fdata[FDF1QUADFLAG] = 1;
           }

        if (fdata[FDF2OBS] == 0.0) fdata[FDF2OBS] = 1.0;
        if (fdata[FDF1OBS] == 0.0) fdata[FDF1OBS] = fdata[FDF2OBS];

        if (fdata[FDF2SW] == 0.0) fdata[FDF2SW] = 1.0;
        if (fdata[FDF1SW] == 0.0) fdata[FDF1SW] = fdata[FDF2SW];

        if ((int)fdata[FDTRANSPOSED] == 0)
           {
            (void) txt2Flt( "X", fdata + FDF2LABEL, SIZE_F2LABEL );
            (void) txt2Flt( "Y", fdata + FDF1LABEL, SIZE_F1LABEL );
            (void) txt2Flt( "Z", fdata + FDF3LABEL, SIZE_F3LABEL );
            (void) txt2Flt( "A", fdata + FDF4LABEL, SIZE_F4LABEL );
           }
        else
           {
            (void) txt2Flt( "Y", fdata + FDF2LABEL, SIZE_F2LABEL );
            (void) txt2Flt( "X", fdata + FDF1LABEL, SIZE_F1LABEL );
            (void) txt2Flt( "Z", fdata + FDF3LABEL, SIZE_F3LABEL );
            (void) txt2Flt( "A", fdata + FDF4LABEL, SIZE_F4LABEL );
           }
       }

    return( 0 );
}
