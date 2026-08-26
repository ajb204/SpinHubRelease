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
/* vectorio: vector-oriented I/O procedures.
/***/

#include <stdio.h>

#include "dataio.h"
#include "vectorio.h"

/***/
/* mv_read: read contiguous vectors in non-contiguous matrix into a contiguous workspace.
/***/

#define FPR (void)fprintf

int mv_read64( dataUnit, eCount1, offset, vCount1, vJump1, vector1 )

   FILE_UNIT( dataUnit );

   int   offset, eCount1, vCount1, vJump1;
   float *vector1;
{
    NMR_INT i;
    int     error;

/***/
/* Force large contiguous block if possible:
/***/

    if (vJump1 == eCount1)
       {
        eCount1 *= vCount1;
        vCount1  = 1;
       }

    for( i = 0; i < vCount1; i++ )
       {
        error = dataPos( dataUnit, sizeof(float)*offset );

        if (error)
           {
            FPR( stderr, "Error positioning file, mv_read\n" );
            goto shutdown;
           }

        error = dataRead( dataUnit, vector1, sizeof(float)*eCount1 );

        if (error)
           {
            FPR( stderr, "Error reading file, mv_read\n" );
            goto shutdown;
           }

        vector1 += eCount1;
        offset  += vJump1;
       }

/***/
/* Exit point:
/***/

shutdown:

    return( error );
}

/***/
/* Write contiguous vectors in contiguous matrix as contiguous vectors in non-contiguous matrix.
/***/

int mv_write64( dataUnit, eCount1, offset, vCount1, vJump1, vector1 )

   FILE_UNIT( dataUnit );

   NMR_INT offset, eCount1, vCount1, vJump1;
   float *vector1;
{
    NMR_INT i;
    int     error;

/***/
/* Force large contiguous block if possible:
/***/

    if (vJump1 == eCount1)
       {
        eCount1 *= vCount1;
        vCount1  = 1;
       }

    for( i = 0; i < vCount1; i++ )
       {
        error = dataPos( dataUnit, sizeof(float)*offset );

        if (error)
           {
            FPR( stderr, "Error positioning file, mv_write\n" );
            goto shutdown;
           }

        error = dataWrite( dataUnit, vector1, sizeof(float)*eCount1 );

        if (error)
           {
            FPR( stderr, "Error writing file, mv_write\n" );
            goto shutdown;
           }

        vector1 += eCount1;
        offset  += vJump1;
       }

/***/
/* Exit point:
/***/

shutdown:

    return( error );
}

/***/
/* Move contiguous vectors from non-contiguous matrix into a contiguous workspace.
/***/

int mv_move64( dataPtr, eCount, offset, vCount, vJump, dest )

   NMR_INT offset, eCount, vCount, vJump;
   float   *dataPtr, *dest;
{
   NMR_INT i, j;
   float   *src;

   dataPtr += offset;

   for( i = 0; i < vCount; i++ )
      {
       src = dataPtr;

       for( j = 0; j < eCount; j++ ) *dest++ = *src++;

       dataPtr += vJump;
      }

   return( 0 );
}

/***/
/* Bottom.
/***/
