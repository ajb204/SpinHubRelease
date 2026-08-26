/***/
/* smxUtil: facilities for manipulating submatrix data.
/***/

#include <stdio.h>

#include "smxutil.h"

static int smxSize[SMX_MAXDIM];
static int smxEdge[SMX_MAXDIM];
static int smxJump[SMX_MAXDIM];
static int smxX1[SMX_MAXDIM];
static int smxXN[SMX_MAXDIM];

static int matSize[SMX_MAXDIM];
static int matJump[SMX_MAXDIM];
static int matX1[SMX_MAXDIM];
static int matXN[SMX_MAXDIM];

static int smxBlocks[SMX_MAXDIM];
static int blockJump[SMX_MAXDIM];

static int srcLoc[SMX_MAXDIM];
static int destLoc[SMX_MAXDIM];

static int dimCount  = 0;
static int blockSize = 0;
static int wordSize  = 0;

static char *smx = (char *)NULL;
static char *mat = (char *)NULL;

static int smx2matrixR();
static int matrix2smxR();

#define ISWAP( IA, IB ) itemp = IA; IA = IB; IB = itemp

/***/
/* smx2matrix: move selected region of data in submatrix format
/*             to a selected region in a regular matrix.
/*
/* 1. Data limits are origin=1.
/*
/* 2. Values listed as "optional" can be given as null pointers,
/*    and will assume full data region.
/*
/* 3. Submatrix sizes must be multiples of edge sizes.
/***/

int smx2matrix( thisSMX,       /* Submatrix input.                         */
                thisMat,       /* Matrix output.                           */
                thisMatSize,   /* Sizes of matrix data dimensions.         */
                thisMatX1,     /* Data Limit of matrix dest, optional.     */
                thisMatXN,     /* Data Limit of matrix dest, optional.     */
                thisSMXSize,   /* Sizes of submatrix data dimensions.      */ 
                thisSMXX1,     /* Data Limit of submatrix src, optional.   */
                thisSMXXN,     /* Data Limit of submatrix src, optional.   */
                thisEdge,      /* Dimensions of submatrix block.           */
                thisWordSize,  /* Word size in bytes of input and output.  */ 
                thisDimCount ) /* Dimension count.                         */

   int  *thisSMXSize, *thisSMXX1, *thisSMXXN, *thisEdge;
   int  *thisMatSize, *thisMatX1, *thisMatXN;
   int  thisWordSize, thisDimCount;

   char *thisSMX, *thisMat;
{
    int  error;

    error = smxInit( thisMatSize, thisMatX1, thisMatXN,
                     thisSMXSize, thisSMXX1, thisSMXXN,
                     thisEdge, thisWordSize, thisDimCount );

    if (error) return( error );

    smx = thisSMX;
    mat = thisMat;

    smx2matrixR( dimCount );

    return( 0 );
}

/***/
/* matrix2smx: move selected region of data in regular matrix format
/*              to a selected region in a submatrix format.
/*
/* 1. Data limits are origin=1.
/*
/* 2. Values listed as "optional" can be given as null pointers,
/*    and will assume full data region.
/*
/* 3. Submatrix sizes must be multiples of edge sizes.
/***/

int matrix2smx( thisMat,       /* Matrix input.                            */
                thisSMX,       /* Submatrix output.                        */
                thisMatSize,   /* Sizes of matrix data dimensions.         */
                thisMatX1,     /* Data Limit of matrix dest, optional.     */
                thisMatXN,     /* Data Limit of matrix dest, optional.     */
                thisSMXSize,   /* Sizes of submatrix data dimensions.      */
                thisSMXX1,     /* Data Limit of submatrix src, optional.   */
                thisSMXXN,     /* Data Limit of submatrix src, optional.   */
                thisEdge,      /* Dimensions of submatrix block.           */
                thisWordSize,  /* Word size in bytes of input and output.  */
                thisDimCount ) /* Dimension count.                         */

   int  *thisSMXSize, *thisSMXX1, *thisSMXXN, *thisEdge;
   int  *thisMatSize, *thisMatX1, *thisMatXN;
   int  thisWordSize, thisDimCount;

   char *thisSMX, *thisMat;
{
    int error;

    error = smxInit( thisMatSize, thisMatX1, thisMatXN,
                     thisSMXSize, thisSMXX1, thisSMXXN,
                     thisEdge, thisWordSize, thisDimCount );

    if (error) return( error );

    smx = thisSMX;
    mat = thisMat;

    matrix2smxR( dimCount );

    return( 0 );
}

/***/
/* smx2matrixR: convert submatrix data by recursion over dimensions.
/***/

static int smx2matrixR( dim )

   int dim;
{
    char *src, *dest;
    int  i;

    dim--;

    for( srcLoc[dim] = smxX1[dim]; srcLoc[dim] <= smxXN[dim]; srcLoc[dim]++ )
       { 
        destLoc[dim] = matX1[dim] + srcLoc[dim] - smxX1[dim];

        if (dim == 0)
           {
            src  = smx + getSMXLoc( srcLoc );
            dest = mat + getMatLoc( destLoc );

            for( i = 0; i < wordSize; i++ ) *dest++ = *src++;
           }
        else
           {
            smx2matrixR( dim );
           }
       }

    return( 0 );
}

/***/
/* matrix2smxR: convert matrix data by recursion over dimensions.
/***/

static int matrix2smxR( dim )

   int dim;
{
    char *src, *dest;
    int  i;

    dim--;

    for( srcLoc[dim] = matX1[dim]; srcLoc[dim] <= matXN[dim]; srcLoc[dim]++ )
       {
        destLoc[dim] = smxX1[dim] + srcLoc[dim] - matX1[dim];

        if (dim == 0)
           {
            src  = mat + getMatLoc( srcLoc );
            dest = smx + getSMXLoc( destLoc );

            for( i = 0; i < wordSize; i++ ) *dest++ = *src++;
           }
        else
           {
            matrix2smxR( dim );
           }
       }

    return( 0 );
}

int smxInit( thisMatSize, thisMatX1, thisMatXN,
             thisSMXSize, thisSMXX1, thisSMXXN,
             thisEdge, thisWordSize, thisDimCount )

   int *thisSMXSize, *thisSMXX1, *thisSMXXN, *thisEdge;
   int *thisMatSize, *thisMatX1, *thisMatXN;
   int thisWordSize, thisDimCount;
{
    int i, j, itemp;

    smx       = (char *)NULL;
    mat       = (char *)NULL;

    dimCount  = 0;
    blockSize = 0;
    wordSize  = 0;

    for( i = 0; i < SMX_MAXDIM; i++ )
       {
        smxSize[i] = 0;
        smxEdge[i] = 0;
        smxJump[i] = 0;
        smxX1[i]   = 0;
        smxXN[i]   = 0;

        matSize[i] = 0;
        matJump[i] = 0;
        matX1[i]   = 0;
        matXN[i]   = 0;

        blockJump[i] = 0;
        smxBlocks[i] = 0;

        srcLoc[i]    = 0;
        destLoc[i]   = 0;
       }

    if (thisDimCount > SMX_MAXDIM) return( SMX_ERR_MAXDIM ); 
    if (thisDimCount < 1) return( 0 ); 

    for( i = 0; i < thisDimCount; i++ )
       {
        smxSize[i] = thisSMXSize[i];
        matSize[i] = thisMatSize[i];
        smxEdge[i] = thisEdge[i];

        matX1[i]   = thisMatX1 ? thisMatX1[i] : 1;
        smxX1[i]   = thisSMXX1 ? thisSMXX1[i] : 1;

        matXN[i]   = thisMatXN ? thisMatXN[i] : matSize[i];
        smxXN[i]   = thisSMXXN ? thisSMXXN[i] : smxSize[i];

        if (smxSize[i] < 1 || matSize[i] < 1 || smxEdge[i] < 1)
           {
            return( SMX_ERR_NULLSIZE );
           }

        if (matX1[i] > matXN[i]) { ISWAP( matX1[i], matXN[i] ); }
        if (smxX1[i] > smxXN[i]) { ISWAP( smxX1[i], smxXN[i] ); }

        if (matX1[i] < 1 || matXN[i] > matSize[i]) return( SMX_ERR_MATLIM );
        if (smxX1[i] < 1 || smxXN[i] > smxSize[i]) return( SMX_ERR_SMXLIM );

        if (smxSize[i] % smxEdge[i]) return( SMX_ERR_EDGE );

        if (smxXN[i] - smxX1[i] != matXN[i] - matX1[i]) 
          {
           return( SMX_ERR_MISMATCH );
          }
       }

    dimCount  = thisDimCount;
    wordSize  = thisWordSize;
    blockSize = wordSize;

    for( i = 0; i < dimCount; i++ ) 
       {
        smxBlocks[i]  = smxSize[i]/smxEdge[i];
        blockSize    *= smxEdge[i];
       }

    for( i = dimCount - 1; i >= 0; i-- )
       {
        matJump[i]   = wordSize;
        smxJump[i]   = wordSize;
        blockJump[i] = blockSize;

        for( j = i - 1; j >= 0; j-- )
           {
            matJump[i]   *= matSize[j];
            smxJump[i]   *= smxEdge[j];
            blockJump[i] *= smxBlocks[j];
           }
       }

    return( 0 );
}

/***/
/* getSMXLoc: find byte offset of a coord in a regular matrix; requires 
/*            use of smxInit() initialization routine above.
/***/

int getMatLoc( iList )

   int *iList;
{
    int i, loc;

    loc = 0;

    for( i = 0; i < dimCount; i++ ) 
       {
        loc += (iList[i] - 1)*matJump[i]; 
       }

    return( loc );
}

/***/
/* getSMXLoc: find byte offset of a coord in a submatrix; requires
/*            use of smxInit() initialization routine above.
/***/

int getSMXLoc( iList )
   
   int *iList;
{
    int loc, i, nDiv, nMod;

    loc = 0;

    for( i = 0; i < dimCount; i++ )
       {
        nDiv = (iList[i] - 1) / smxEdge[i];
        nMod = (iList[i] - 1) % smxEdge[i];

        loc += nDiv*blockJump[i] + nMod*smxJump[i];
       }

   return( loc );
}
