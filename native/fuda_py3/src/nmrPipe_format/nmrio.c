/******************************************************************************/
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
/* nmrio: procedures for I/O with NMR format files.
/***/

#include <stdio.h>
#include <unistd.h>

#include "getstat.h"
#include "memory.h"
#include "dataio.h"
#include "fdatap.h"
#include "nmrio.h"
#include "vectorio.h"
#include "inquire.h"
#include "testsize.h"
#include "prec.h"

#define FPR (void)fprintf

/***/
/* rdNMRPipeData: Read spectral data from an NMRPipe-format file.
/*                Accounts for byte-swap. 
/***/

int  rdNMRPipeData64( inName,    /* Name of NMRPipe format file to read.     */
                      array,     /* On return, will contain data read.       */
                      offset,    /* Additional offset past header in points. */
                      length )   /* Number of points to read.                */

   NMR_INT  offset, length;
   char    *inName;
   float   *array;
{
    float  fdata[FDATASIZE]; 
    int    error;

    FILE_UNIT( inUnit );

/***/
/* Open the file.
/* Read the header (sets byteswap state).
/* Read the requested data.
/* Close the file.
/***/

    if ((error = dataOpen( inName, &inUnit, FB_READ )))
       {
        FPR( stderr, "NMR I/O Error Opening File %s\n", inName );
        return( error );
       }

    if ((error = rdFDATAU( inUnit, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Header in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = dataPos( inUnit, sizeof(float)*(FDATASIZE + offset) )))
       {
        FPR( stderr, "NMR I/O Error Positioning in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = dataRead( inUnit, array, sizeof(float)*length )))
       {
        FPR( stderr, "NMR I/O Error Reading Data in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    (void) dataClose( inUnit );

    return( 0 );
}

/***/
/* wrNMRPipeData: Read spectral data from an NMRPipe-format file.
/*                Accounts for byte-swap. 
/***/

int  wrNMRPipeData64( inName,    /* Name of NMRPipe format file to to write. */
                      array,     /* Data to write.                           */
                      offset,    /* Additional offset past header in points. */
                      length )   /* Number of points to write.               */

   NMR_INT offset, length;
   char    *inName;
   float   *array;
{
    float  fdata[FDATASIZE]; 
    int    error;

    FILE_UNIT( inUnit );

/***/
/* Open the file.
/* Read the header (sets byteswap state).
/* Write the requested data.
/* Close the file.
/***/

    if ((error = dataOpen( inName, &inUnit, FB_READWRITE )))
       {
        FPR( stderr, "NMR I/O Error Opening File %s\n", inName );
        return( error );
       }

    if ((error = rdFDATAU( inUnit, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Header in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = dataPos( inUnit, sizeof(float)*(FDATASIZE + offset) )))
       {
        FPR( stderr, "NMR I/O Error Positioning in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = dataWrite( inUnit, array, sizeof(float)*length )))
       {
        FPR( stderr, "NMR I/O Error Reading Data in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    (void) dataClose( inUnit );

    return( 0 );
}

/***/
/* truncNMR: reduce file size to match prediction from header.
/***/

int truncNMR( inName, vFlag )

   char *inName;
   int  vFlag;
{
    int     zSize, aSize, iz, ia, dimCount, pipeFlag, eCount, error;
    char    thisName[NAMELEN+1];
    float   fdata[FDATASIZE];

    struct  FileSize predSizeInfo, trueSizeInfo;
    NMR_INT predBytes, origBytes;

    eCount    = 0;
    predBytes = 0;
    origBytes = 0;

    (void) sprintf( thisName, inName, 1, 1 );

    if ((error = rdFDATA( thisName, fdata )))
       {
        FPR( stderr, "NMR I/O Truncate Error Reading Header in File %s\n", thisName );
        return( 1 );
       }

    dimCount = getParmI( fdata, FDDIMCOUNT, NULL_DIM );
    pipeFlag = getParmI( fdata, FDPIPEFLAG, NULL_DIM );

    zSize = dimCount < 3 ? 1 : getParmI( fdata, NDSIZE, CUR_ZDIM );
    aSize = dimCount < 4 ? 1 : getParmI( fdata, NDSIZE, CUR_ADIM );

    if (zSize < 1) zSize = 1;
    if (aSize < 1) aSize = 1;

    (void) getNMRBytes( fdata, NMR_SINGLE_FILE, &predSizeInfo );

    predBytes = predSizeInfo.iTotalBytes;

    if (isAZFmt( inName ))
       {
        for( ia = 1; ia <= aSize; ia++ )
           {
            for( iz = 1; iz <= zSize; iz++ )
               {
                (void) sprintf( thisName, inName, ia, iz );
                if (dataTrunc( thisName, predBytes, &origBytes, vFlag )) eCount++;
               }
           }
       }
    else if (isZFmt( inName ))
       {
        if (dimCount == 4) zSize = aSize;

        for( iz = 1; iz <= zSize; iz++ )
           {
            (void) sprintf( thisName, inName, iz );
            if (dataTrunc( thisName, predBytes, &origBytes, vFlag )) eCount++;
           }
       }
    else
       {
        if (dataTrunc( inName, predBytes, &origBytes, vFlag )) eCount++;
       }

    if (eCount)
       {
        FPR( stderr, "NMR I/O Truncate Error: %d Errors while truncating %s\n", eCount, inName );
        return( 1 );
       }

    return( 0 );
}

/***/
/* Read 2D file region, given a file name.
/***/

int  rdNMRi64( array,     /* On return, will contain data read.              */
               inName,    /* Name of NMR format file to read.                */
               xStart,    /* Origin 1 start of region to read, X-Axis.       */
               xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to read, Y-Axis.       */
               yLength )  /* Length of region to read, Y-Axis, 1->SPECNUM.   */

   NMR_INT xStart, xLength, yStart, yLength;
   char    *inName;
   float   *array;
{
    float  fdata[FDATASIZE]; 
    int    error;

    FILE_UNIT( inUnit );

/***/
/* Open the file.
/* Read the header.
/* Read the requested region.
/* Close the file.
/***/

    if ((error = dataOpen( inName, &inUnit, FB_READ )))
       {
        FPR( stderr, "NMR I/O Error Opening File %s\n", inName );
        return( error );
       }

    if ((error = rdFDATAU( inUnit, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Header in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = rdNMRiU64( array, inUnit, xStart, xLength, yStart, yLength, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Data from %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    (void) dataClose( inUnit );

    return( 0 );
}

/***/
/* rdNMRiZ: read 2D file region, given a file name; non-folding version.
/***/

int  rdNMRiZ( array,     /* On return, will contain data read.              */
              inName,    /* Name of NMR format file to read.                */
              xStart,    /* Origin 1 start of region to read, X-Axis.       */
              xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE. */
              yStart,    /* Origin 1 start of region to read, Y-Axis.       */
              yLength )  /* Length of region to read, Y-Axis, 1->SPECNUM.   */

   int   xStart, xLength, yStart, yLength;
   char  *inName;
   float *array;
{
    float  fdata[FDATASIZE]; 
    int    error;

    FILE_UNIT( inUnit );

/***/
/* Open the file.
/* Read the header.
/* Read the requested region.
/* Close the file.
/***/

    if ((error = dataOpen( inName, &inUnit, FB_READ )))
       {
        FPR( stderr, "NMR I/O Error Opening File %s\n", inName );
        return( error );
       }

    if ((error = rdFDATAU( inUnit, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Header in File %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    if ((error = rdNMRiZU( array, inUnit, xStart, xLength, yStart, yLength, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Data from %s\n", inName );
        (void) dataClose( inUnit );
        return( error );
       }

    (void) dataClose( inUnit );

    return( 0 );
}

/***/
/* rdNMRi3D: read 3D file region, given a file name list.
/*           Folds along Z Axis if Z Region exceeds list bounds. 
/***/

int  rdNMRi3D( array,     /* On return, will contain data read.              */
               inList,    /* List of NMR format files to read.               */
               fileCount, /* Number of files listed in inList.               */
               xStart,    /* Origin 1 start of region to read, X-Axis.       */
               xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to read, Y-Axis.       */
               yLength,   /* Length of region to read, Y-Axis, 1->SPECNUM.   */
               zStart,    /* Origin 1 start of region to read, Z-Axis.       */
               zLength )  /* Length of region to read, Z-Axis.               */

   int   fileCount, xStart, xLength, yStart, yLength, zStart, zLength;
   char  **inList;
   float *array;
{
    int  iz, error;
    char *inName;

/***/
/* For each file in the specifed Z-Axis limits:
/*   Read the requested region.
/***/
    
    for( iz = zStart - 1; iz < (zStart + zLength - 1); iz++ )
       {
        if (iz < 0)
           inName = *(inList + fileCount + (iz % fileCount));
        else
           inName = *(inList + (iz % fileCount));
         
        if ((error = rdNMRi( array, inName, xStart, xLength, yStart, yLength )))
           {
            FPR( stderr, "NMR I/O Error Reading Data from %s\n", inName );
            return( error );
           }

        array += xLength*yLength;
       }

    return( 0 );
}

/***/
/* rdNMRi4D: read 4D file region, given a file name list.
/***/

int  rdNMRi4D( array,     /* On return, will contain data read.              */
               inList,    /* List of NMR format files to read.               */
               zSize,     /* Length of Z-Axis dimension of spectrum.         */
               aSize,     /* Length of A-Axis dimension of spectrum.         */
               xStart,    /* Origin 1 start of region to read, X-Axis.       */
               xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to read, Y-Axis.       */
               yLength,   /* Length of region to read, Y-Axis, 1->SPECNUM.   */
               zStart,    /* Origin 1 start of region to read, Z-Axis.       */
               zLength,   /* Length of region to read, Z-Axis.               */
               aStart,    /* Origin 1 start of region to read, A-Axis.       */
               aLength )  /* Length of region to read, A-Axis.               */

   int   zSize, aSize, xStart, xLength, yStart, yLength, zStart, zLength,
         aStart, aLength;

   char  **inList;
   float *array;
{
    int  ia, foldA, error;

    for( ia = aStart - 1; ia < (aStart + aLength - 1); ia++ )
       {
        if (ia < 0)
           foldA = aSize + (ia % aSize);
        else
           foldA = ia % aSize;

        error = rdNMRi3D( array,
                          (inList + foldA*zSize),
                          zSize,
                          xStart, xLength,
                          yStart, xLength,
                          zStart, zLength );

        if (error) return( 1 );

        array += xLength*yLength*zLength;
       }

    return( 0 );
}

/***/
/* Read a 2D file region, given a file unit and a header, folds along Y-Axis if needed.
/* Adjusted to allow pipe-format single-file data.
/***/

int  rdNMRiU64( array,     /* On return, will contain data read.               */
                inUnit,    /* File unit of NMR format file to read.            */
                xStart,    /* Origin 1 start of region to read, X-Axis.        */
                xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE.  */
                yStart,    /* Origin 1 start of region to read, Y-Axis.        */
                yLength,   /* Length of region to read, Y-Axis, 1->SPECNUM.    */
                fdata )    /* Accurate header of file connected to inUnit.     */

   float   *array, fdata[FDATASIZE];
   NMR_INT xStart, xLength, yStart, yLength;

   FILE_UNIT( inUnit );
{
    NMR_INT size, specnum, zSize, aSize, jump, foldY, iy, offset;
    int     pipeFlag, dimCount, quadState, error;

/***/
/* Extract size parameters.
/***/

    specnum  = (NMR_INT)fdata[FDSPECNUM];
    size     = (NMR_INT)fdata[FDSIZE];
    zSize    = (NMR_INT)fdata[FDF3SIZE];
    aSize    = (NMR_INT)fdata[FDF4SIZE];

    dimCount = (int)fdata[FDDIMCOUNT];
    pipeFlag = (int)fdata[FDPIPEFLAG];

    if (specnum < 1)   specnum = 1;
    if (dimCount == 1) specnum = 1;

    if (pipeFlag)
       {
        if (dimCount >= 3 && zSize > 0) specnum *= zSize;
        if (dimCount >= 4 && zSize > 0 && aSize > 0) specnum *= aSize;
       }

    if (fdata[FDQUADFLAG] == 1.0)
       quadState = 1;
    else
       quadState = 2;

    if (xStart < 1 || (xStart + xLength - 1) > size*quadState)
       {
#ifdef NMR64
        FPR( stderr, "NMR I/O Error: X Limits (%ld %ld of %ld).\n", (long)xStart, (long)(xStart + xLength - 1), (long)size*quadState );
#else
        FPR( stderr, "NMR I/O Error: X Limits (%d %d of %d).\n", xStart, xStart + xLength - 1, size*quadState );
#endif
        return( 1 );
       }

/***/
/* If no clipping is required, read the region directly.
/***/

    jump = size*quadState;

    if (yStart > 0 && (yStart + yLength - 1) <= specnum)
       {
        offset = FDATASIZE + (yStart - 1)*jump + xStart - 1;
        error  = mv_read64( inUnit, xLength, offset, yLength, jump, array );
        return( error );
       }
    
/***/
/* If clipping is required:
/*    Read rows one at a time, adjusting for folding.
/***/

    for( iy = yStart - 1; iy < (yStart + yLength - 1); iy++ )
       {
        if (iy < 0)
           foldY = specnum + (iy % specnum);
        else
           foldY = (iy % specnum);

        offset = FDATASIZE + foldY*jump + xStart - 1;

        if ((error = mv_read64( inUnit, xLength, offset, 1, jump, array )))
           {
            FPR( stderr, "NMR I/O Error Reading Folded Y Axis.\n" );
            return( error );
           }

        array += xLength;
       }

    return( error );
}

/***/
/* rdNMRiZU: read a 2D file region, given a file unit and a header;
/*           non-folding version. 
/*
/*           Adjusted to allow pipe-format single-file data.
/***/

int  rdNMRiZU( array,     /* On return, will contain data read.               */
               inUnit,    /* File unit of NMR format file to read.            */
               xStart,    /* Origin 1 start of region to read, X-Axis.        */
               xLength,   /* Length of region to read, X-Axis, 1->QUAD*SIZE.  */
               yStart,    /* Origin 1 start of region to read, Y-Axis.        */
               yLength,   /* Length of region to read, Y-Axis, 1->SPECNUM.    */
               fdata )    /* Accurate header of file connected to inUnit.     */

   float *array, fdata[FDATASIZE];
   int   xStart, xLength, yStart, yLength;

   FILE_UNIT( inUnit );
{
    int size, specnum, zSize, aSize, pipeFlag, dimCount, quadState, offset, error;

/***/
/* Extract size parameters.
/***/

    specnum  = fdata[FDSPECNUM];
    size     = fdata[FDSIZE];
    zSize    = fdata[FDF3SIZE];
    aSize    = fdata[FDF4SIZE];

    dimCount = fdata[FDDIMCOUNT];
    pipeFlag = fdata[FDPIPEFLAG];

    if (specnum < 1)   specnum = 1;
    if (dimCount == 1) specnum = 1;

    if (pipeFlag)
       {
        if (dimCount >= 3 && zSize > 0) specnum *= zSize;
        if (dimCount >= 4 && zSize > 0 && aSize > 0) specnum *= aSize;
       }

    if (fdata[FDQUADFLAG] == 1.0)
       quadState = 1;
    else
       quadState = 2;

    if (xStart < 1 || (xStart + xLength - 1) > size*quadState)
       {
        FPR( stderr, "NMR I/O Error: X Limits.\n" );
        return( 1 );
       }

    offset = FDATASIZE + (yStart - 1)*size*quadState + xStart - 1;
    error  = mv_read( inUnit, xLength, offset, yLength, size*quadState, array );

    return( error );
}
    
/***/
/* Write a 2D file region, given a file name.
/***/

int  wrNMRi64( array,     /* Should contain data to write.                    */
               outName,   /* Name of NMR format file to write.                */
               xStart,    /* Origin 1 start of region to write, X-Axis.       */
               xLength,   /* Length of region to write, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to write, Y-Axis.       */
               yLength )  /* Length of region to write, Y-Axis, 1->SPECNUM.   */

   NMR_INT xStart, xLength, yStart, yLength;
   char    *outName;
   float   *array;
{
    float  fdata[FDATASIZE]; 
    int    error;

    FILE_UNIT( outUnit );

/***/
/* Open the file.
/* Read the header.
/* Write the requested region.
/* Close the file.
/***/

    if ((error = dataOpen( outName, &outUnit, FB_READWRITE )))
       {
        FPR( stderr, "NMR I/O Error Opening File %s\n", outName );
        return( error );
       }

    if ((error = rdFDATAU( outUnit, fdata )))
       {
        FPR( stderr, "NMR I/O Error Reading Header in File %s\n", outName );
        (void) dataClose( outUnit );
        return( error );
       }

    if ((error = wrNMRiU64( array, outUnit, xStart, xLength, yStart, yLength, fdata )))
       {
        FPR( stderr, "NMR I/O Error Writing Data to %s\n", outName );
        (void) dataClose( outUnit );
        return( error );
       }

    (void) dataClose( outUnit );

    return( 0 );
}

/***/
/* wrNMRi3D: write a 3D file region, given a file name list.
/***/

int  wrNMRi3D( array,     /* Should contain data to write.                    */
               outList,   /* List of NMR format files to write.               */
               fileCount, /* Number of files listed in outList.               */
               xStart,    /* Origin 1 start of region to write, X-Axis.       */
               xLength,   /* Length of region to write, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to write, Y-Axis.       */
               yLength,   /* Length of region to write, Y-Axis, 1->SPECNUM.   */
               zStart,    /* Origin 1 start of region to write, Z-Axis.       */
               zLength )  /* Length of region to write, Z-Axis.               */

   int   fileCount, xStart, xLength, yStart, yLength, zStart, zLength;
   char  **outList;
   float *array;
{
    int  iz, error;
    char *outName;

    if (zStart < 1 || (zStart + zLength - 1) > fileCount)
       {
#ifdef NMR64
        FPR( stderr, "NMR I/O Error with Z Bounds %ld %ld\n", (long)zStart, (long)zLength );
#else
        FPR( stderr, "NMR I/O Error with Z Bounds %d %d\n", zStart, zLength );
#endif
        return( 1 );
       }

/***/
/* For each file in the specifed Z-Axis limits:
/*   Write the requested region.
/***/
  
    for( iz = zStart - 1; iz < (zStart + zLength - 1); iz++ )
       {
        outName = *(outList + iz);
         
        if ((error = wrNMRi( array, outName, xStart, xLength, yStart, yLength )))
           {
            FPR( stderr, "NMR I/O Error Writing Data to %s\n", outName );
            return( error );
           }

        array += xLength*yLength;
       }

    return( 0 );
}

/***/
/* wrNMRi4D: write 4D file region, given a file name list.
/***/

int  wrNMRi4D( array,     /* Matrix of data to write.                         */
               outList,   /* List of NMR format files to write.               */
               zSize,     /* Length of Z-Axis dimension of spectrum.          */
               aSize,     /* Length of A-Axis dimension of spectrum.          */
               xStart,    /* Origin 1 start of region to write, X-Axis.       */
               xLength,   /* Length of region to write, X-Axis, 1->QUAD*SIZE. */
               yStart,    /* Origin 1 start of region to write, Y-Axis.       */
               yLength,   /* Length of region to write, Y-Axis, 1->SPECNUM.   */
               zStart,    /* Origin 1 start of region to write, Z-Axis.       */
               zLength,   /* Length of region to write, Z-Axis.               */
               aStart,    /* Origin 1 start of region to write, A-Axis.       */
               aLength )  /* Length of region to, write, A-Axis.              */

   int   zSize, aSize, xStart, xLength, yStart, yLength, zStart, zLength,
         aStart, aLength;

   char  **outList;
   float *array;
{
    int  ia, error;

    if (aStart < 1 || aStart > aSize)
       {
#ifdef NMR64
        FPR( stderr, "NMR I/O Error with A Bounds %ld %ld\n", (long)aStart, (long)aLength );
#else
        FPR( stderr, "NMR I/O Error with A Bounds %d %d\n", aStart, aLength );
#endif
        return( 1 );
       }

    for( ia = aStart - 1; ia < (aStart + aLength - 1); ia++ )
       {
        error = wrNMRi3D( array,
                          (outList + ia*zSize),
                          zSize,
                          xStart, xLength,
                          yStart, xLength,
                          zStart, zLength );

        if (error) return( 1 );

        array += xLength*yLength*zLength;
       }

    return( 0 );
}

/***/
/* Write a 2D file region, given a file unit and a header.
/***/

int  wrNMRiU64( array,     /* Should contain data to write.                     */
                outUnit,   /* File Unit of NMR format file to write.            */
                xStart,    /* Origin 1 start of region to write, X-Axis.        */
                xLength,   /* Length of region to write, X-Axis, 1->QUAD*SIZE.  */
                yStart,    /* Origin 1 start of region to write, Y-Axis.        */
                yLength,   /* Length of region to write, Y-Axis, 1->SPECNUM.    */
                fdata )    /* Accurate header of file connected to outUnit.     */

   float   *array, fdata[FDATASIZE];
   NMR_INT xStart, xLength, yStart, yLength;

   FILE_UNIT( outUnit );
{
    NMR_INT size, specnum, zSize, aSize, jump, offset;
    int     pipeFlag, dimCount, quadState, error;

/***/
/* Extract size parameters.
/***/

    specnum  = (NMR_INT)fdata[FDSPECNUM];
    size     = (NMR_INT)fdata[FDSIZE];
    zSize    = (NMR_INT)fdata[FDF3SIZE];
    aSize    = (NMR_INT)fdata[FDF4SIZE];

    dimCount = (int)fdata[FDDIMCOUNT];
    pipeFlag = (int)fdata[FDPIPEFLAG];

    if (specnum < 1)   specnum = 1;
    if (dimCount == 1) specnum = 1;

    if (fdata[FDQUADFLAG] == 1.0)
       quadState = 1;
    else
       quadState = 2;

    if (pipeFlag)
       {
        if (dimCount >= 3 && zSize > 0) specnum *= zSize;
        if (dimCount >= 4 && zSize > 0 && aSize > 0) specnum *= aSize;
       }

    if (xStart < 1 || (xStart + xLength - 1) > size*quadState)
       {
        FPR( stderr, "NMR I/O Error: X Limits.\n" );
        return( 1 );
       }

/***/
/* Write the region:
/***/

    jump   = size*quadState;
    offset = FDATASIZE + (yStart - 1)*jump + xStart - 1;

    error  = mv_write64( outUnit, xLength, offset, yLength, jump, array );

    return( error );
}

/***/
/* Find and set display scale values of matrix.
/***/

int getNMRMinMax64( matrix, size, specnum, quadSize, fdata )

   float   *matrix, fdata[FDATASIZE];
   NMR_INT size, specnum;
   int     quadSize;
{
    float minVal, maxVal;

    minVal = *matrix;
    maxVal = *matrix;

    while( specnum-- )
       {
        (void) minMax264( matrix, size, &minVal, &maxVal );
        matrix += size*quadSize;
       }

    fdata[FDMIN]       = minVal;
    fdata[FDMAX]       = maxVal;
    fdata[FDDISPMIN]   = minVal;
    fdata[FDDISPMAX]   = maxVal;
    fdata[FDSCALEFLAG] = 1.0;

    return( 0 );
}

/***/
/* scaleNMR: reset display scaling variables of named NMR file.
/***/

int scaleNMR( inName )

   char  *inName;
{
   int   error;
   float thisMin, thisMax;

   if (!inName)    return( 1 );
   if (!inName[0]) return( 1 );

   thisMin = 0.0;
   thisMax = 0.0;

   error = scaleNMR2( inName, 1, 0, 1, &thisMin, &thisMax );

   return( error );
}

/***/
/* scaleNMR2: reset display scaling variables of named NMR file.
/***/

int scaleNMR2( inName, realFlag, imagFlag, saveFlag, minVal, maxVal )

   char  *inName;
   int   realFlag, imagFlag, saveFlag;
   float *minVal, *maxVal;
{
    float *matrix, fdata[FDATASIZE];
    int   byteSwapState, autoSwapState, swapDone, error;
    int   ix1, nx, size, specnum, zSize, aSize, quadState, pipeFlag, cubeFlag, dimCount;

    *minVal = 0.0;
    *maxVal = 0.0;

    if (!fileExists( inName )) return( 1 );

    autoSwapState = getAutoSwapFlag();
    byteSwapState = getByteSwapFlag();

    (void) setByteSwapFlag( 0 );

    error = rdFDATAS( inName, fdata, &swapDone );

    if (error) goto shutdown;

    size     = fdata[FDSIZE];
    specnum  = fdata[FDSPECNUM];
    zSize    = fdata[FDF3SIZE];
    aSize    = fdata[FDF4SIZE];

    dimCount = fdata[FDDIMCOUNT];
    pipeFlag = fdata[FDPIPEFLAG];
    cubeFlag = fdata[FDCUBEFLAG];
 
    if (specnum < 1)   specnum = 1;
    if (dimCount == 1) specnum = 1;

    if (cubeFlag)
       {
        if (dimCount >= 4 && zSize > 0 && aSize > 0) specnum *= zSize; 
       }
    else if (pipeFlag)
       {
        if (dimCount >= 3 && zSize > 0) specnum *= zSize; 
        if (dimCount >= 4 && zSize > 0 && aSize > 0) specnum *= aSize; 
       }

    if (fdata[FDQUADFLAG] == 1.0)
       quadState = 1;
    else
       quadState = 2;

    ix1 = 1;
    nx  = size;

    if (quadState == 2)
       {
        if (realFlag)
           {
            if (imagFlag)
               {
                size *= 2;
                nx    = size;
               }
           }
        else
           {
            ix1 += size;
           }
       }

    if (!(matrix = fltAlloc( "scaleNMR", size*specnum )))
       {
        error = 1;
        goto shutdown;
       }

    if ((error = rdNMRi( matrix, inName, ix1, nx, 1, specnum ))) goto shutdown;

    (void) getNMRMinMax64( matrix, (NMR_INT)size, (NMR_INT)specnum, 1, fdata );

    *minVal = fdata[FDMIN];
    *maxVal = fdata[FDMAX];

    if (saveFlag)
       {
        if (swapDone) (void) swapHdr( fdata );
        if ((error = wrFDATA( inName, fdata ))) goto shutdown;
       }

shutdown:

    (void) setAutoSwapFlag( autoSwapState );
    (void) setByteSwapFlag( byteSwapState );

    (void) deAlloc( "scaleNMR", matrix, sizeof(float)*size*specnum );

    return( 0 ); 
}

/***/
/* scaleNMR3D: reset display scaling of 3D NMR list.
/***/

int scaleNMR3D( inList, zSize )

   char **inList;
   int  zSize;
{
    int i;

    for( i = 0; i < zSize; i++ )
       {
        if (scaleNMR( inList[i] )) return( 1 );
       }

    return( 0 );
}

/***/
/* Bottom.
/***/
