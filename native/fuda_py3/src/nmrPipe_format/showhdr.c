/*********************************************************************/
/*                           NMRPipe                                 */
/*                     Copyright 1992-2016                           */
/*                        Frank Delaglio                             */
/*********************************************************************/

/***/
/* showHdr: display file parameters by current dimension.
/***/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define FDATA_LOCLIST

#include "cmndargs.h"
#include "testsize.h"
#include "specunit.h"
#include "fdatap.h"
#include "nmrio.h"
#include "inquire.h"
#include "prec.h"

static char axisNames[] = "XYZABC";
static char masterPrefix[NAMELEN+1], prefix[NAMELEN+1];

#define PR (void)printf
 
int main( argc, argv )

   int  argc;
   char **argv;
{
    int   i, j, n, verbose, fixFlag, dumpFlag, testFlag, listFlag, swapDone, error;
    char  *sPtr, inName[NAMELEN+1], templateName[NAMELEN+1], thisPrefix[NAMELEN+1];
    float fdata[FDATASIZE], tmpFDATA[FDATASIZE];

/***/
/* Extract command line arguments.
/*
/* For each file argument:
/*  Read the header.
/*  Fix the head if needed.
/*  Display the header contents.
/***/

    error    = 0;
    swapDone = 0;
    fixFlag  = 1;
    dumpFlag = 0;
    testFlag = 0;
    listFlag = 0;
    n        = 0;

    masterPrefix[0] = '\0';
    prefix[0]       = '\0';
    inName[0]       = '\0';
    templateName[0] = '\0';

    (void) initDataIO();
 
    if (error = getparms( argc, argv, &fixFlag, &dumpFlag, &testFlag, &listFlag, &verbose ))
       {
        (void) printf( "%sShowHDR Error extracting command line arguments.\n", prefix );
        return( error );
       }

    for( i = 1; i < argc; i++ )
       {
        (void) setHdrNull( fdata );

        if (!strcmp( argv[i], "-prefix" )) 
           {
            i++;
            continue;
           }

        if (!isFlag( argv[i] ))  /* No getNthArg needed. */
           {
            sPtr = getNthArg( argc, argv, i );

            (void) strcpy( templateName, sPtr );
            (void) sprintf( inName, sPtr, 1, 1 );

            if (!strcasecmp( masterPrefix, "Auto" ))
               (void) sprintf( prefix, "%s: ", inName );
            else
               (void) strcpy( prefix, masterPrefix );

            if (!fileExists( inName ))
               {
                if (testFlag)
                   (void) printf( "%svalid 0 input %s 1 1 %s ndim 0 sizeList 0 0 0 0 qSize 0", n ? "\n" : "", sPtr, inName );
                else
                   (void) printf( "ShowHdr Error finding input file %s\n\n", inName );
                   
                n++;
                continue;
               }

            if (error = rdFDATAS( inName, fdata, &swapDone ))
               {
                if (testFlag)
                   (void) printf( "%svalid 0 input %s 1 1 %s ndim 0 sizeList 0 0 0 0 qSize 0", n ? "\n" : "", sPtr, inName );
                else
                   (void) printf( "%sShowHdr Error reading file %s\n", prefix, inName );

                n++;
                continue;
               }

            if (testFlag) 
               {
                (void) copyHdr( tmpFDATA, fdata );

                if (HDR_BAD == testHdr( tmpFDATA ))
                   {
                    (void) printf( "%svalid 0 input %s 1 1 %s ndim 0 sizeList 0 0 0 0 qSize 0", n ? "\n" : "", sPtr, inName );

                    n++;
                    continue;
                   }
               }

            if (fixFlag) (void) fixfdata( fdata );

            if (testFlag)
               {
                if (n) (void) printf( "\n" );
                (void) showTestInfo( templateName, inName, fdata );
               }
            else if (listFlag)
               {
                (void) showListInfo( templateName, inName, fdata );
               }
            else
               {
                (void) showHdr( templateName, inName, fdata, swapDone, verbose );
                (void) dumpHdr( inName, fdata, dumpFlag );
               }

            n++;
           }
       }

    (void) fflush( stdout );

    return( 0 );
}

/***/
/* showListInfo: list all names associated with an NMR template.
/***/

int showListInfo( templateName, inName, fdata )

   char  *templateName, *inName;
   float fdata[FDATASIZE];
{
    int  n, dimCount, fmtCount, fileCount, pipeFlag, cubeFlag, iz, ia, xSize, ySize, zSize, aSize;
    char thisName[NAMELEN+1];

    dimCount  = getParmI( fdata, FDDIMCOUNT,  NULL_DIM );
    pipeFlag  = getParmI( fdata, FDPIPEFLAG,  NULL_DIM );
    cubeFlag  = getParmI( fdata, FDCUBEFLAG,  NULL_DIM );
    fileCount = getParmI( fdata, FDFILECOUNT, NULL_DIM );

    fmtCount  = getFmtCount( templateName );

    if (dimCount == 3)
       {
        if (fmtCount == 0 && pipeFlag == 0) dimCount = 2;
       }
    else if (dimCount == 4)
       {
        if (fmtCount == 0 && pipeFlag == 0)
           dimCount = 2;
        else if (fmtCount == 0 && cubeFlag != 0)
           dimCount = 3;
       }

    xSize = dimCount >= 1 ? getParmI( fdata, NDSIZE, CUR_XDIM ) : 0;
    ySize = dimCount >= 2 ? getParmI( fdata, NDSIZE, CUR_YDIM ) : 1;
    zSize = dimCount >= 3 ? getParmI( fdata, NDSIZE, CUR_ZDIM ) : 1;
    aSize = dimCount >= 4 ? getParmI( fdata, NDSIZE, CUR_ADIM ) : 1;

    if (fmtCount == 2)
       {
        for( ia = 1; ia <= aSize; ia++ )
           {
            for( iz = 1; iz <= zSize; iz++ )
               {
                (void) sprintf( thisName, templateName, ia, iz );
                (void) printf( "%s\n", thisName );
               }
           }
       }
    else if (fmtCount == 1)
       {
        if (dimCount == 4)
           {
            for( ia = 1; ia <= aSize; ia++ )
               {
                (void) sprintf( thisName, templateName, ia );
                (void) printf( "%s\n", thisName );
               }
           }
        else
           {
            for( iz = 1; iz <= zSize; iz++ )
               {
                (void) sprintf( thisName, templateName, iz );
                (void) printf( "%s\n", thisName );
               }
           }
       }
    else
       {
        (void) printf( "%s\n", inName );
       }

    return( 0 );
}

/***/
/* showTestInfo: show basic info to characterize an NMRPipe-format input file.
/***/

int showTestInfo( templateName, inName, fdata )

   char  *templateName, *inName;
   float fdata[FDATASIZE]; 
{
    int n, dimCount, fmtCount, fileCount, validFlag, pipeFlag, cubeFlag, iz, ia, xSize, ySize, zSize, aSize, qSize, zFT, aFT;

    validFlag = 1;

    n = strlen( templateName );

    if (n > 2)
       {
        if (templateName[0] == '.' && templateName[1] == '/') templateName += 2;
       }

    n = strlen( inName );

    if (n > 2)
       {
        if (inName[0] == '.' && inName[1] == '/') inName += 2;
       }

    dimCount  = getParmI( fdata, FDDIMCOUNT,  NULL_DIM );
    pipeFlag  = getParmI( fdata, FDPIPEFLAG,  NULL_DIM );
    cubeFlag  = getParmI( fdata, FDCUBEFLAG,  NULL_DIM );
    fileCount = getParmI( fdata, FDFILECOUNT, NULL_DIM );

    fmtCount  = getFmtCount( templateName );

    if (dimCount == 3)
       {
        if (fmtCount == 0 && pipeFlag == 0) dimCount = 2;
       }
    else if (dimCount == 4)
       {
        if (fmtCount == 0 && pipeFlag == 0) 
           dimCount = 2;
        else if (fmtCount == 0 && cubeFlag != 0)
           dimCount = 3;
       }
 
    xSize     = dimCount >= 1 ? getParmI( fdata, NDSIZE, CUR_XDIM ) : 0;
    ySize     = dimCount >= 2 ? getParmI( fdata, NDSIZE, CUR_YDIM ) : 1;
    zSize     = dimCount >= 3 ? getParmI( fdata, NDSIZE, CUR_ZDIM ) : 1;
    aSize     = dimCount >= 4 ? getParmI( fdata, NDSIZE, CUR_ADIM ) : 1;

    qSize     = getQuad( fdata, FDQUADFLAG, NULL_DIM ); 

    zFT       = dimCount >= 3 ? getParmI( fdata, NDFTFLAG, CUR_ZDIM ) : 0;
    aFT       = dimCount >= 4 ? getParmI( fdata, NDFTFLAG, CUR_ADIM ) : 0;
    iz        = 1;
    ia        = 1;

    if (dimCount >= 3 && zFT != 0)
       {
        iz = zSize/2;
        if (iz < 1) iz = 1;
       }

    if (dimCount >= 4 && aFT != 0)
       {
        ia = aSize/2;
        if (ia < 1) ia = 1;
       }

    (void) printf( "valid %d input %s %d %d %s ndim %d sizeList %d %d %d %d qSize %d fileCount %d", 
                   validFlag, templateName, iz, ia, inName, dimCount, xSize, ySize, zSize, aSize, qSize, fileCount );

    return( 0 );
}

/***/
/* dumpHdr: dump all header values as floats.
/***/

int dumpHdr( inName, fdata, dumpFlag )

   float fdata[FDATASIZE];
   char  *inName;
   int   dumpFlag;
{
    char *sPtr;
    int  i;

    if (!dumpFlag) return( 0 );

    (void) printf( "%sFDATA Values for %s:\n", prefix, inName );

    sPtr = (char *)NULL;

    for( i = 0; i < FDATASIZE; i++ )
       {
        (void) getNameByVal( &sPtr, fdataLocList, (float)i );
        (void) printf( "%s   %3d. % 18.6f %s\n", prefix, i, fdata[i], sPtr ? sPtr : " " );
       }

    (void) printf( "%s\n", prefix );

    return( 0 );
}

/***/
/* showHdr: extract and display header information.
/***/

int showHdr( templateName, inName, fdata, swapDone, verbose )

   float fdata[FDATASIZE];
   int   swapDone, verbose;
   char  *templateName, *inName;
{
    int   xSize, ySize, zSize, aSize, planes, quadState, itemp, dim,
          firstPlane, lastPlane, pipeMode, cubeFlag, dimCount, origDimCount, fileCount,
          fmtCount, threadCount, threadID, nusDim, hour, min, sec, month, 
          day, year, planeType, scans, n;

    NMR_INT sliceCount64;

    struct FileSize fInfo;
    float  rtemp, rtemp2;
    char   *stemp, ctemp[NAMELEN+1];
     
/***/
/* General:
/***/

    dimCount  = getParmI( fdata, FDDIMCOUNT,  0 );
    nusDim    = getParmI( fdata, FDNUSDIM,    0 );
    pipeMode  = getParmI( fdata, FDPIPEFLAG,  0 );
    cubeFlag  = getParmI( fdata, FDCUBEFLAG,  0 );
    planeType = getParmI( fdata, FD2DPHASE,   0 );
    fileCount = getParmI( fdata, FDFILECOUNT, 0 );

    fmtCount     = getFmtCount( templateName );
    origDimCount = dimCount;

    (void) printf( "%s%s: %s DIM: %d", prefix, fmtCount ? "TEMPLATE" : "FILE", templateName, dimCount );

    if (nusDim) (void) printf( " NUSDIM: %d", nusDim );

    (void) printf( " QUAD: %s",  getParmStr( fdata, FDQUADFLAG, 0 ));

    if (dimCount == 1)
       {
        (void) printf( " 2DMODE: %s", planeType == 0 ? "None" : getParmStr( fdata, FD2DPHASE, 0 ));
       }
    else
       {
        (void) printf( " 2DMODE: %s", getParmStr( fdata, FD2DPHASE, 0 ));
       }

    (void) printf( " %s", getParmStr( fdata, FDTRANSPOSED, 0 ));

    (void) printf( "\n" );

    (void) getFileBytes( inName, &fInfo );

    if (fInfo.status == FILE_STATUS_INT_OK)
       (void) printf( "%sBYTES: %s", prefix, getNMRIntStr( fInfo.iTotalBytes ));
    else if (fInfo.status == FILE_STATUS_INT_TOO_SMALL)
       (void) printf( "%sBYTES: %.0lf", prefix, fInfo.dTotalBytes );
    else
       (void) printf( "%sBYTES: Unknown", prefix );

    (void) getNMRBytes( fdata, NMR_SINGLE_FILE, &fInfo );

    if (fInfo.status == FILE_STATUS_INT_OK)
       (void) printf( " PRED: %s", getNMRIntStr( fInfo.iTotalBytes ));
    else if (fInfo.status == FILE_STATUS_INT_TOO_SMALL)
       (void) printf( " PRED: %.0lf", fInfo.dTotalBytes );
    else
       (void) printf( " PRED: Unknown" );

    if (fmtCount)
       {
        (void) getNMRBytes( fdata, NMR_FILE_SERIES, &fInfo );

        if (fInfo.status == FILE_STATUS_INT_OK)
           (void) printf( " TOTAL: %s", getNMRIntStr( fInfo.iTotalBytes ));
        else if (fInfo.status == FILE_STATUS_INT_TOO_SMALL)
           (void) printf( " TOTAL: %.0lf", fInfo.dTotalBytes );
        else
           (void) printf( " TOTAL: Unknown" );
       }
    
    (void) printf( " MIN: %g",   getParm( fdata, FDMIN, 0 ));
    (void) printf( " MAX: %g",   getParm( fdata, FDMAX, 0 ));
    (void) printf( " VALID: %d", getParmI( fdata, FDSCALEFLAG, 0 ));

    (void) printf( "\n" );
    (void) printf( "%sORDER:", prefix );

    for( dim = 0; dim < dimCount; dim++ )
       {
        (void) printf( " %d", getParmI(fdata,FDDIMORDER+dim,0));
       }

    (void) printf( " PIPE: %d", pipeMode );
    (void) printf( " CUBE: %d", cubeFlag );

    xSize  = getParmI( fdata, NDSIZE, CUR_XDIM );

    ySize  = getParmI( fdata, NDSIZE, CUR_YDIM );
    ySize  = ySize < 1 ? 1 : ySize;

    zSize  = getParmI( fdata, NDSIZE, CUR_ZDIM );
    zSize  = zSize < 1 ? 1 : zSize;

    aSize  = getParmI( fdata, NDSIZE, CUR_ADIM );
    aSize  = aSize < 1 ? 1 : aSize;

    planes = zSize*aSize;

    if (1 == getParmI( fdata, FDQUADFLAG, 0 ))
       quadState = 1;
    else
       quadState = 2;

    (void) printf( " FILES: %d", fileCount );

    if (dimCount == 1)
       {
        (void) printf( " %dx%d", xSize, quadState );
       }
    else if (dimCount == 2)
       {
        (void) printf( " %dx%dx%d", xSize, ySize, quadState );
       }
    else if (dimCount == 3)
       {
        if (pipeMode || fmtCount == 1)
           (void) printf( " %dx%dx%dx%d", xSize, ySize, zSize, quadState );
        else
           (void) printf( " %dx%dx%d", xSize, ySize, quadState );
       }
    else if (dimCount == 4)
       {
        if (pipeMode || fmtCount)
           (void) printf( " %dx%dx%dx%dx%d", xSize, ySize, zSize, aSize, quadState );
        else if (cubeFlag)
           (void) printf( " %dx%dx%dx%d", xSize, ySize, zSize, quadState );
        else
           (void) printf( " %dx%dx%d", xSize, ySize, quadState );
       } 

    if (dimCount == 3)
       {
        if (pipeMode)
           (void) printf ( " 3D Stream" );
        else if (fmtCount == 1)
           (void) printf ( " Plane Series" );
        else
           (void) printf ( " 2D Plane" );
       }
    else if (dimCount == 4)
       {
        if (pipeMode)
           (void) printf ( " 4D Stream" );
        else if (fmtCount == 2)
           (void) printf ( " Plane Series" );
        else if (fmtCount == 1)
           (void) printf ( " Cube Series" );
        else if (cubeFlag)
           (void) printf ( " 3D Cube" );
        else
           (void) printf ( " 2D Plane" );
       }

    if (0.0 != (rtemp = getParm( fdata, FDTAU, 0 )))
       {
        (void) printf( " TAU: %g", rtemp );
       }

    if (swapDone) (void) printf( " Swapped" );

    (void) printf( "\n" );

/***/
/* Parallel Processing:
/***/

    firstPlane   = getParmI( fdata, FDFIRSTPLANE,  0 );
    lastPlane    = getParmI( fdata, FDLASTPLANE,   0 );
    threadID     = getParmI( fdata, FDTHREADID,    0 );
    threadCount  = getParmI( fdata, FDTHREADCOUNT, 0 );

    sliceCount64 = getFDATASliceCount64( fdata );

    if (firstPlane || lastPlane)
       {
#ifdef NMR64
        (void) printf( "%sFIRST PLANE: %d LAST PLANE: %d SLICE COUNT: %ld\n",
                       prefix, firstPlane, lastPlane, sliceCount64 );
#else
        (void) printf( "%sFIRST PLANE: %d LAST PLANE: %d SLICE COUNT: %d\n",
                       prefix, firstPlane, lastPlane, sliceCount64 );
#endif
       }

   if (threadCount || threadID)
       {
        (void) printf( "%sTHREAD COUNT: %d THREAD ID: %d\n", prefix, threadCount, threadID );
       }

/***/
/* Date and time of conversion:
/***/

    month = getParmI( fdata, FDMONTH, 0 );
    day   = getParmI( fdata, FDDAY,   0 );
    year  = getParmI( fdata, FDYEAR,  0 );
    hour  = getParmI( fdata, FDHOURS, 0 );
    min   = getParmI( fdata, FDMINS,  0 );
    sec   = getParmI( fdata, FDSECS,  0 );

    if (verbose)
       {
        n = 0;

        month = getParmI( fdata, FDMONTH, 0 );
        day   = getParmI( fdata, FDDAY,   0 );
        year  = getParmI( fdata, FDYEAR,  0 );
        hour  = getParmI( fdata, FDHOURS, 0 );
        min   = getParmI( fdata, FDMINS,  0 );
        sec   = getParmI( fdata, FDSECS,  0 );

        if (month && day && year)
           {
            if (!n) (void) printf( "%s", prefix );

            (void) printf( "CONVERTED: %d/%d/%d %d:%02d %s",
                           month, day, year,
                           hour > 12 ? hour - 12 : hour,
                           min, hour > 12 ? "PM" : "AM" );

            n++;
           }

        scans = getParmI( fdata, FDSCANS, 0 );

        if (scans)
           {
            (void) printf( "%s", n ? " " : prefix );
            (void) printf( "SCANS: %d", scans );
            n++;
           }

        rtemp = getParm( fdata, FDTEMPERATURE, 0 );

        if (rtemp != 0.0)
           {
            (void) printf( "%s", n ? " " : prefix );
            (void) printf( "T: %.1f", rtemp );
            n++;
           }

        (void) strcpy( ctemp, getParmStr( fdata, FDTITLE, NULL_DIM ));

        if (strlen( ctemp ))
           {
            (void) printf( "%s", n ? " " : prefix );
            (void) printf( "TITLE: %s\n", ctemp );
            n++;
           }

        (void) strcpy( ctemp, getParmStr( fdata, FDCOMMENT, NULL_DIM ));

        if (strlen( ctemp ))
           {
            (void) printf( "%s", n ? " " : prefix );
            (void) printf( "COMMENT: %s\n", ctemp );
            n++;
           }
       }

/***/
/* Show Axis Labels For Each Column:
/***/

    if (nusDim > dimCount) dimCount = nusDim;

    (void) printf( "%s\n", prefix );
    (void) printf( "%s%10s", prefix, " " );

    for( dim = 0; dim < dimCount; dim++ )
       {
        (void) printf( "     %c%-7s ", axisNames[dim], "-Axis" );
       }

    (void) printf( "\n" );

/***/
/* Actual sizes:
/***/
 
    (void) printf( "%s\n", prefix );
    (void) printf( "%s%-11s", prefix, "DATA SIZE:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDSIZE, dim );
        (void) printf( " %12d", itemp );
       }

/***/
/* Apodization sizes:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "APOD SIZE:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPOD, dim );
        (void) printf( " %12d", itemp );
       }

/***/
/* Spectral Width:
/***/

    if (planeType == FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "WIDTH:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDSW, dim );

            if (dim > 2)
               {
                if (!strcasecmp( "IR", getParmStr( fdata, NDLABEL, dim )))
                   {
                    rtemp /= C_CM_SEC;
                   }
               }

            if (fabs( (double)rtemp ) > 9.0e5)
               (void) printf( " %12e", rtemp );
            else
               (void) printf( " %12.3f", rtemp );
           }
       }
    else
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "SW Hz:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDSW, dim );
    
            if (fabs( (double)rtemp ) > 9.0e5)
               (void) printf( " %12e", rtemp );
            else
               (void) printf( " %12f", rtemp );
           }
       }

    if (verbose && planeType != FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "SW PPM:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp  = getParm( fdata, NDSW,  dim );
            rtemp2 = getParm( fdata, NDOBS, dim );
            rtemp  = rtemp2 == 0.0 ? 0.0 : rtemp/rtemp2;

            (void) printf( "    % 9.3f", rtemp );
           }

        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "Hz/POINT:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp  = getParm( fdata, NDSW, dim );
            rtemp2 = getParm( fdata, NDSIZE, dim );
            rtemp  = rtemp2 == 0.0 ? 0.0 : rtemp/rtemp2;

            (void) printf( " %12.3f", rtemp );
           }

        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "AQTIME SEC:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp  = getParm( fdata, NDAPOD, dim );
            rtemp2 = getParm( fdata, NDSW,   dim );
            rtemp  = rtemp2 == 0.0 ? 0.0 : rtemp/rtemp2;

            (void) printf( " %12f", rtemp );
           }
       }

/***/
/* Observe Frequency:
/***/

    if (planeType != FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "OBS MHz:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDOBS, dim );
            (void) printf( " %12f", rtemp );
           }
       }

/***/
/* Axis Origin:
/***/

    if (planeType != FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "ORIG Hz:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDORIG, dim );
    
            if (fabs( (double)rtemp ) > 9.0e5)
               (void) printf( " %12e", rtemp );
            else
               (void) printf( " %12f", rtemp );
           }
       }

/***/
/* FT Flag:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "DOMAIN:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        stemp = getParmStr( fdata, NDFTFLAG, dim );
        (void) printf( " %12s", stemp );
       }

/***/
/* Quad Flag:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "MODE:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        stemp = getParmStr( fdata, NDQUADFLAG, dim );
        (void) printf( " %12s", stemp );
       }

/****/
/* Axis labels:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "NAME:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        stemp = getParmStr( fdata, NDLABEL, dim );
        (void) printf( " %12s", stemp );
       }

/***/
/* Begin of verbose mode listing:
/***/

    if (!verbose)
       {
        (void) printf( "\n%s\n", prefix );
        return( 0 );
       }

/***/
/* Apodization Parameters:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "APOD NAME:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        stemp = getParmStr( fdata, NDAPODCODE, dim );
        (void) printf( " %12s", stemp );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "APOD Q1:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDAPODQ1,   dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "APOD Q2:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDAPODQ2,   dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "APOD Q3:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDAPODQ3,   dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }
    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "LB:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDLB,       dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "GB:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDGB,       dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "GOFF:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = getParm( fdata, NDGOFF,     dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }


    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "C1:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDAPODCODE, dim );
        rtemp = 1.0 + getParm( fdata, NDC1, dim );

        if (itemp)
           (void) printf( "    % 9.3f", rtemp );
        else
           (void) printf( " %12s", "None" );
       }

/***/
/* Processing Sizes:
/*   Original Time Domain Size.
/*   Zero Fill Amount
/*   Size at Fourier Transform.
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "TD SIZE:" );
 
    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDTDSIZE, dim );

        if (itemp == 0)
           (void) printf( " %12s", "None" );
        else
           (void) printf( " %12d", itemp );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "ZF SIZE:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDZF, dim );

        if (itemp > 0)
           (void) printf( " %11dx", itemp );
        else if (itemp == 0)
           (void) printf( " %12s", "None" );
        else
           (void) printf( " %12d", -itemp );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "FT SIZE:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDFTSIZE, dim );

        if (itemp == 0)
           (void) printf( " %12s", "None" );
        else
           (void) printf( " %12d", itemp );
       }

/***/
/* Phase Parameters:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "P0:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        rtemp = getParm( fdata, NDP0, dim );
        (void) printf( "    % 9.3f", rtemp );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "P1:" );

    for( dim = 1; dim <= dimCount; dim++ )
       {
        rtemp = getParm( fdata, NDP1, dim );
        (void) printf( "    % 9.3f", rtemp );
       }

/***/
/* Extraction Regions:
/***/

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "EXT X1:" );
 
    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDX1, dim );
 
        if (itemp > 0)
           (void) printf( " %12d", itemp );
        else
           (void) printf( " %12s", "None" );
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "EXT XN:" );
 
    for( dim = 1; dim <= dimCount; dim++ )
       {
        itemp = getParmI( fdata, NDXN, dim );
 
        if (itemp > 0)
           (void) printf( " %12d", itemp );
        else
           (void) printf( " %12s", "None" );
       }
 
/***/
/* Frequency Domain Limits:
/***/

    if (planeType == FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "FIRST:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            if (dim > 2)
               {
                if (!strcasecmp( "IR", getParmStr( fdata, NDLABEL, dim )))
                   rtemp = iPnt2spec( fdata, dim, 1, LAB_WN );
                else
                   rtemp = iPnt2spec( fdata, dim, 1, LAB_PPM );
               }
            else
               {
                rtemp = iPnt2spec( fdata, dim, 1, LAB_MM );
               }

            (void) printf( "    % 9.3f", rtemp );
           }

        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "LAST:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            itemp = getParmI( fdata, NDSIZE, dim );

            if (dim > 2)
               {
                if (!strcasecmp( "IR", getParmStr( fdata, NDLABEL, dim )))
                   rtemp = iPnt2spec( fdata, dim, itemp, LAB_WN );
                else
                   rtemp = iPnt2spec( fdata, dim, itemp, LAB_PPM );
               }
            else
               {
                rtemp = iPnt2spec( fdata, dim, itemp, LAB_MM );
               }

            (void) printf( "    % 9.3f", rtemp );
           }
       }
    else
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "PPM FIRST:" );
 
        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = iPnt2spec( fdata, dim, 1, LAB_PPM );
            (void) printf( "    % 9.3f", rtemp );
           }

        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "PPM LAST:" );
 
        for( dim = 1; dim <= dimCount; dim++ )
           {
            itemp = getParmI( fdata, NDSIZE, dim );
            rtemp = iPnt2spec( fdata, dim, itemp, LAB_PPM );
            (void) printf( "    % 9.3f", rtemp );
           }
       }

    if (planeType != FD_IMAGE)
       {
        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "PPM OFFSET:" );

        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDOFFPPM, dim );
            (void) printf( "    % 9.3f", rtemp );
           }

        (void) printf( "\n" );
        (void) printf( "%s%-11s", prefix, "PPM CAR:" );
 
        for( dim = 1; dim <= dimCount; dim++ )
           {
            rtemp = getParm( fdata, NDCAR, dim );
            (void) printf( "    % 9.3f", rtemp );
           }
       }

    (void) printf( "\n" );
    (void) printf( "%s%-11s", prefix, "CENTER:" );
 
    for( dim = 1; dim <= dimCount; dim++ )
       {
        rtemp = getParm( fdata, NDCENTER, dim );
        itemp = rtemp == ZERO_EQUIV ? 0 : (int)rtemp;
 
        (void) printf( " %12d", itemp );
       }
 
/***/
/* End:
/***/

    (void) printf( "\n%s\n", prefix );

    return( 0 );
}

/***/
/* getparms: extract command line arguments for showHdr.
/***/

int getparms( argc, argv, fixFlag, dumpFlag, testFlag, listFlag, verbose )

   int  argc, *fixFlag, *dumpFlag, *testFlag, *listFlag, *verbose;
   char **argv;
{
    *dumpFlag = 0;
    *fixFlag  = 1;
    *testFlag = 0;
    *listFlag = 0;
    *verbose  = 0;

    masterPrefix[0] = '\0';
    prefix[0]       = '\0';

    if (argc == 1 || flagLoc( argc, argv, "-help" ))
       {
        PR( "ShowHdr: Display NMR File Header Information:\n" );
        PR( "%s [-nofix] [-verb] [fileList]\n", argv[0] );
        PR( " -prefix  preText  Prefix for Output Lines or Keyword: Auto None.\n" );
        PR( " -nofix            Suppress Adjustment of Input Header.\n" );
        PR( " -dump             Dump All Header Values as Floats.\n" );
        PR( " -test             Test for Valid Data, Report Input Info.\n" );
        PR( " -list             List File Names in Template Series.\n" );
        PR( " -verb             Print Verbose Parameter Listing.\n" );
        return( 0 );
       }

    (void) strArgD( argc, argv, "-prefix", masterPrefix );

    if (flagLoc( argc, argv, "-nofix" )) *fixFlag  = 0;
    if (flagLoc( argc, argv, "-dump" ))  *dumpFlag = 1;
    if (flagLoc( argc, argv, "-test" ))  *testFlag = 1;
    if (flagLoc( argc, argv, "-list" ))  *listFlag = 1;
    if (flagLoc( argc, argv, "-verb" ))  *verbose  = 1;

    if (!strcasecmp( masterPrefix, "None" )) masterPrefix[0] = '\0';

    (void) strcpy( prefix, masterPrefix );

    return( 0 );
}

/***/
/* Bottom.
/***/
