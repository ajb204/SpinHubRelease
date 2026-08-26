/***/
/* Read single-file NMRPipe-format data. 
/***/

#include <stdio.h>
#include <string.h>

#include "fdatap.h"
#include "dataio.h"
#include "inquire.h"
#include "memory.h"
#include "prec.h"
#include "readnmr.h"

/***/
/* Allocate and read entire matrix from single-file data:
/*   Space is allocated and data returned in matPtr.
/*   Effective dimension count is returned in dimCountPtr.
/*   Sizes of each dimension are returned in sizeList.
/*   Total number of points is returned in totalPts.
/*   Quad size is returned in qSizePtr (1 Real, 2 Complex). 
/*   If the data is complex, sizeList[XLOC] will be doubled.
/***/

int readNMR( char *inName, float fdata[FDATASIZE], float **matPtr, int *sizeList, NMR_INT *totalPts, int *qSizePtr, int *dimCountPtr )
{
    int     i, dimCount, pipeFlag, cubeFlag, quadSize, xQuadSize, yQuadSize, error;

    float   *rPtr;
    UNIT    inUnit;
    NMR_INT n;

    inUnit       = UNIT_NULL;
    error        = 0;

    *matPtr      = (float *)NULL;
    *totalPts    = 0;
    *qSizePtr    = 0;
    *dimCountPtr = 0;

    for( i = 0; i < MAXDIM; i++ ) sizeList[i] = 0;

    if (!fileExists( inName )) return( 1 );

    if ((error = dataOpen( inName, &inUnit, FB_READ ))) return( 2 );

    if ((error = rdFDATAU( inUnit, fdata ))) return( 3 );

    for( i = 0; i < MAXDIM; i++ ) sizeList[i] = 1;

    dimCount = getParm( fdata, FDDIMCOUNT, NULL_DIM );
    pipeFlag = getParm( fdata, FDPIPEFLAG, NULL_DIM );
    cubeFlag = getParm( fdata, FDCUBEFLAG, NULL_DIM );

    quadSize  = getQuad( fdata, FDQUADFLAG, NULL_DIM );
    xQuadSize = getQuad( fdata, NDQUADFLAG, CUR_XDIM );
    yQuadSize = getQuad( fdata, NDQUADFLAG, CUR_YDIM );

    if (dimCount == 1)
       {
        sizeList[XLOC] = getParm( fdata, NDSIZE, CUR_XDIM );
       }
    else if (dimCount == 2)
       {
        sizeList[XLOC] = getParm( fdata, NDSIZE, CUR_XDIM );
        sizeList[YLOC] = getParm( fdata, NDSIZE, CUR_YDIM );
       }
    else if (dimCount == 3)
       {
        sizeList[XLOC] = getParm( fdata, NDSIZE, CUR_XDIM );
        sizeList[YLOC] = getParm( fdata, NDSIZE, CUR_YDIM );

        if (pipeFlag)
           sizeList[ZLOC] = getParm( fdata, NDSIZE, CUR_ZDIM );
        else
           dimCount = 2;
       } 
    else if (dimCount == 4)
       {
        sizeList[XLOC] = getParm( fdata, NDSIZE, CUR_XDIM );
        sizeList[YLOC] = getParm( fdata, NDSIZE, CUR_YDIM );

        if (pipeFlag)
           {
            sizeList[ZLOC] = getParm( fdata, NDSIZE, CUR_ZDIM );
            sizeList[ALOC] = getParm( fdata, NDSIZE, CUR_ADIM );
           }
        else if (cubeFlag)
           {
            sizeList[ZLOC] = getParm( fdata, NDSIZE, CUR_ZDIM );
            dimCount = 3;
           }
        else
           {
            dimCount = 2;
           }
       }

    if (quadSize == 2) sizeList[XLOC] *= 2;

    if (dimCount == 2 && sizeList[YLOC] == 1) dimCount = 1;

    n = 1;

    for( i = 0; i < dimCount; i++ ) n *= sizeList[i];

    if (!(rPtr = fltAlloc( "nmr", n ))) return( 4 );

    *matPtr      = rPtr;
    *dimCountPtr = dimCount;
    *totalPts    = n;
    *qSizePtr    = quadSize;

    if ((error = dataRead( inUnit, rPtr, sizeof(float)*n ))) return( 5 );

    return( 0 );
}

