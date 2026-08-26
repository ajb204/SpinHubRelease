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
/* tp2DMMM: in-memory 2D transpose of 2D NMR data matrix.
/*        (Miserable Memory-hungry Method).
/*
/* matrix: matrix to transpose; on return, will contain
/*         transposed result.
/*
/* xSize: length of an all-real row in matrix.
/*
/* ySize: length of an all-real column in matrix.
/*
/* zSize: number of 2D planes in matrix; 
/*
/*        1 = all-real.
/*        2 = complex.
/*        4 = hyper-complex.
/***/

#include <stdlib.h>
#include <stddef.h>

#include "memory.h"
#include "prec.h"
#include "vutil.h"

int tp2DMMM( matrix, xSize, ySize, zSize )

   float *matrix;
   int    xSize, ySize, zSize;

{
    float   *xPtr, *yPtr, *zPtr, *destPtr, *work;
    int     ix, iy, iz, error;
    NMR_INT total, xInc, yInc;

/***/
/* Initialize:
/***/

    error   = 0;
    work    = (float *)NULL;
    destPtr = (float *)NULL;
    xPtr    = (float *)NULL;
    yPtr    = (float *)NULL;
    zPtr    = (float *)NULL;

    total = (NMR_INT)xSize*ySize*zSize;

    if (!(work = fltAlloc( "tp2DMMM", total )))
       {
        error = 1;
        goto shutdown;
       }

/***/
/* Move original data into work array.
/* Adjust for hypercomplex format: { V1 V2 V3 V4 } -> { V1 V3 V2 V4 }
/***/

    xPtr    = matrix;
    destPtr = work;
    xInc    = (NMR_INT)2*xSize;
    yInc    = (NMR_INT)4*xSize;

    if (zSize == 4)
       {
        for( iy = 0; iy < ySize; iy++ )
           {
            (void) vvCopy( destPtr, xPtr, xSize );
            (void) vvCopy( destPtr + xInc, xPtr + xSize, xSize ); 
            (void) vvCopy( destPtr + xSize, xPtr + xInc, xSize );
            (void) vvCopy( destPtr + xSize + xInc, xPtr + xSize + xInc, xSize );

            destPtr += yInc;
            xPtr    += yInc;
           }
       }
    else
       {
        (void) vvCopy( destPtr, xPtr, total );
       }

/***/
/* Move data in work matrix back to original matrix, in transposed order:
/***/

    yInc = (NMR_INT)xSize*zSize;

    destPtr = matrix;
    xPtr    = work;

    for( ix = 0; ix < xSize; ix++ )
       {
        zPtr = xPtr;

        for( iz = 0; iz < zSize; iz++ )
           {
            yPtr = zPtr;

            for( iy = 0; iy < ySize; iy++ )
               {
                *destPtr++ = *yPtr;
                yPtr += yInc;
               }

            zPtr += xSize;
           }

        xPtr++;
       }

/***/
/* Exit point:
/***/

shutdown:

#ifdef S_SPLINT_S
    (void) free( (void *)work );
#else
    (void) deAlloc( "tp2DMMM", work, total );
#endif

    return( error );
}

/***/
/* Bottom.
/***/
