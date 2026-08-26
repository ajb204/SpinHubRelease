/******************************************************************************/
/*                   ---- NIH NMR Software System ----                        */
/*                        Copyright 1992 to 2015                              */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/******************************************************************************/

/***/
/* fdataio: procedures for I/O of NMRPipe format headers.
/***/

#include <stdio.h>

#include "memory.h"
#include "dataio.h"
#include "nmrio.h"
#include "fdatap.h"
#include "prec.h"

#define FPR (void)fprintf

/***/
/* rdFDATA: read a file header, given a file name.
/***/

int rdFDATA( inName, fdata )

   char  *inName;
   float fdata[FDATASIZE];
{
   FILE_UNIT( inUnit );

    if (dataOpen( inName, &inUnit, FB_READ )) return( 1 );

    if (rdFDATAU( inUnit, fdata ))
       {
        (void) dataClose( inUnit );
        return( 1 );       
       }

    (void) dataClose( inUnit );

    return( 0 );
}

/***/
/* rdFDATAS: read a file header, given a file name; auto-swap version.
/***/

int rdFDATAS( inName, fdata, swapDone )

   char  *inName;
   float fdata[FDATASIZE];
   int   *swapDone;
{
   FILE_UNIT( inUnit );

   int error;

    if (dataOpen( inName, &inUnit, FB_READ )) return( 1 );

    error = readHdr( inUnit, fdata, swapDone );

    (void) dataClose( inUnit );

    return( error );
}

/***/
/* rdFDATAU: read file header, given a file unit.
/***/

int rdFDATAU( inUnit, fdata )

   FILE_UNIT( inUnit ); 

   float fdata[FDATASIZE];
{
    int swapDone, error;

    error    = 0;
    swapDone = 0;

    if ((error = dataPos( inUnit, 0 ))) return( error );

    if ((error = readHdr( inUnit, fdata, &swapDone )))
       {
        return( error );
       }
 
    return( 0 ); 
}

/***/
/* wrFDATA: write a file header, given a file name.
/***/

int wrFDATA( outName, fdata )

   char  *outName;
   float fdata[FDATASIZE];
{
    FILE_UNIT( outUnit );

    if (dataOpen( outName, &outUnit, FB_WRITE )) return( 1 );

    if (wrFDATAU( outUnit, fdata ))
       {
        (void) dataClose( outUnit );
        return( 1 );
       }

    (void) dataClose( outUnit );

    return( 0 );
}

/***/
/* wrFDATAU: write a file header, given a file unit.
/***/

int wrFDATAU( outUnit, fdata )

   FILE_UNIT( outUnit ); 

   float fdata[FDATASIZE];
{
    if (dataPos( outUnit, 0 ))
       {
        return( 1 );
       }

    if (dataWrite( outUnit, fdata, sizeof(float)*FDATASIZE ))
       {
        return( 1 );
       }
 
    return( 0 ); 
}

/***/
/* getNMRsize: interpret file sizes given a file name.
/***/

int getNMRsize( inName, size, specnum, quadState )

   int  *size, *specnum, *quadState;
   char *inName;
{
    float fdata[FDATASIZE];
    int   dimCount;

    if (rdFDATA( inName, fdata )) return( 1 );
    
    *size    = fdata[FDSIZE];
    *specnum = fdata[FDSPECNUM];
    dimCount = fdata[FDDIMCOUNT];

    if (dimCount == 1) *specnum = 1; 
    if (*specnum < 1)  *specnum = 1;
 
    if (fdata[FDQUADFLAG] == 1.0)
       *quadState = 1;
    else
       *quadState = 2;
 
    return( 0 );
}

/***/
/* getNMRsizeND: interpret file size of an ND file list.
/***/

int getNMRsizeND( inList, inCount, sizeList, dimCount, sDimCount, quadState )

   int  *sizeList, *sDimCount, *quadState, inCount, dimCount;
   char **inList;
{
    float  fdata[FDATASIZE];
    char   *cPtr;
    int    i;

    cPtr = *inList;

    if (rdFDATA( cPtr, fdata )) return( 1 );

    *sDimCount  = fdata[FDDIMCOUNT];
    sizeList[0] = fdata[FDSIZE];
    sizeList[1] = fdata[FDSPECNUM];

    switch( dimCount )
       {
        case 1:
           break;
        case 2:
           if (sizeList[1] <= 1 && inCount > 1) sizeList[1] = inCount;
           break; 
        case 3:
           sizeList[2] = inCount;
           break;
        case 4:
           sizeList[2] = fdata[FDF3SIZE];
           sizeList[3] = fdata[FDF4SIZE];
           break;
        default:
           (void) printf( "NMR I/O Error: Max 4D exceeded.\n" );
           return( 1 );
       }

    if (fdata[FDQUADFLAG] == 1.0)
       *quadState = 1;
    else
       *quadState = 2;

    for( i = 0; i < dimCount; i++ )
       {
        if (sizeList[i] < 1) sizeList[i] = 1;
       }

    return( 0 );
}

/***/
/* Bottom.
/***/
