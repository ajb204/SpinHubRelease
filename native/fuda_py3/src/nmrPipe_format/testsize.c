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
/* testSize: utilities for testing sizes of input and output files.
/***/

#include <stdio.h>
#include <limits.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/file.h>
#include <sys/stat.h>

#define FPR (void)fprintf

#define USE_VFS

#ifdef WINNT
#include <sys/mount.h>
#include <sys/stat.h>
#undef USE_VFS
#endif

#ifdef MAC_OSX
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/mount.h>
#undef USE_VFS
#endif

#ifdef SOLARIS
#include <sys/statfs.h>
#define STATFS statvfs
#endif

#ifdef SGI
#include <sys/statfs.h>
#undef USE_VFS
#endif

#ifdef CRAY
#include <sys/statfs.h>
#undef USE_VFS
#endif

#ifdef ALPHA
#include <sys/mount.h>
#include <sys/ustat.h>
#undef USE_VFS
#endif

#ifdef IBM
#include <sys/statfs.h>
#endif

#ifdef WINXP
#include <sys/statvfs.h>
#define STATFS statvfs
#else
#ifdef USE_VFS 
#include <sys/vfs.h>
#endif
#endif

#ifndef STATFS
#define STATFS statfs
#endif

#include "prec.h"
#include "fdatap.h"
#include "testsize.h"


int getFSCheck(); //AJB 14/09/22
int safe_stat(); //AJB 14/09/22

/***/
/* testNMRin: tests that true size of input file matches header.
/***/

int testNMRin( path, fdata, multiFile, verbose )

   char  *path;
   float fdata[FDATASIZE];
   int   multiFile, verbose;
{
    struct FileSize pInfo, tInfo;

/***/
/* Extract the predicted and true file sizes.
/* Stop if either extraction fails.
/***/

    (void) getNMRBytes( fdata, multiFile, &pInfo );
    (void) getFileBytes( path, &tInfo );

    if (pInfo.status == FILE_STATUS_FAIL || tInfo.status == FILE_STATUS_FAIL || pInfo.dTotalBytes < 0.0 || tInfo.dTotalBytes < 0.0)
       {
        return( FILE_SIZE_FAIL );
       }

/***/
/* If the true and predicted sizes can be expressed as ints:
/*   Compare the int sizes to see if they are consistant.
/* Otherwise:
/*   Compare the double-precision float sizes instead.
/***/

    if (pInfo.status == FILE_STATUS_INT_OK && tInfo.status == FILE_STATUS_INT_OK)
       {
        if (pInfo.iTotalBytes == tInfo.iTotalBytes)
           {
            return( FILE_SIZE_OK );
           }
        else if (tInfo.iTotalBytes < pInfo.iTotalBytes)
           {
            if (verbose)
               {
                FPR( stderr, "\n" );
                FPR( stderr, "File Size Check Error:" );
                FPR( stderr, " Insufficient Data in %s\n", path );
                FPR( stderr, "Predicted Size: %s bytes.\n", getNMRIntStr( pInfo.iTotalBytes ));
                FPR( stderr, "Available Size: %s bytes.\n", getNMRIntStr( tInfo.iTotalBytes ));
                FPR( stderr, "\n" );
               }

            return( FILE_SIZE_TOO_SMALL );
           }
        else
           {
            FPR( stderr, "\n" );
            FPR( stderr, "File Size Check Warning:" );
            FPR( stderr, " Extra Data in %s\n", path );
            FPR( stderr, "Predicted Size: %s bytes.\n", getNMRIntStr( pInfo.iTotalBytes ));
            FPR( stderr, "Available Size: %s bytes.\n", getNMRIntStr( tInfo.iTotalBytes ));

            return( FILE_SIZE_TOO_LARGE );
           }
       }
    else
       {
        if (pInfo.dTotalBytes == tInfo.dTotalBytes)
           {
            return( FILE_SIZE_OK );
           }
        else if (tInfo.dTotalBytes < pInfo.dTotalBytes)
           {
            if (verbose)
               {
                FPR( stderr, "\n" );
                FPR( stderr, "File Size Check Error:" );
                FPR( stderr, " Insufficient Data in %s\n", path );
                FPR( stderr, "Predicted Size: %.0lf bytes.\n", pInfo.dTotalBytes );
                FPR( stderr, "Available Size: %.0lf bytes.\n", tInfo.dTotalBytes );
                FPR( stderr, "\n" );
               }

            return( FILE_SIZE_TOO_SMALL );
           }
        else
           {
            if (verbose)
               {
                FPR( stderr, "\n" );
                FPR( stderr, "File Size Check Warning:" );
                FPR( stderr, " Extra Data in %s\n", path );
                FPR( stderr, "Predicted Size: %.0lf bytes.\n", pInfo.dTotalBytes );
                FPR( stderr, "Available Size: %.0lf bytes.\n", tInfo.dTotalBytes );
               }

            return( FILE_SIZE_TOO_LARGE );
           }
       }

    return( FILE_SIZE_OK );
}

/***/
/* testNMRout: tests if predicted length of output exceeds space limits.
/***/

int testNMRout( path, fdata, multiFile, verbose )

   char  *path;
   float fdata[FDATASIZE];
   int   multiFile, verbose;
{
    struct FileSize pInfo, aInfo;
    int    noSpaceFlag;

/***/
/* Extract the predicted and available file sizes.
/* Stop if either extraction fails.
/***/

    if (!getFSCheck()) return( FILE_SIZE_OK );

    (void) getNMRBytes( fdata, multiFile, &pInfo );
    (void) getFreeBytes( path, &aInfo );

    if (pInfo.status == FILE_STATUS_FAIL || aInfo.status == FILE_STATUS_FAIL || pInfo.dTotalBytes < 0.0 || aInfo.dTotalBytes < 0.0)
       {
        FPR( stderr, "\n" );
        FPR( stderr, "Error Checking File Space for %s\n", path );
        FPR( stderr, "To suppress this test, use environment variable settings:\n" );
        FPR( stderr, "\n" );
        FPR( stderr, "   setenv NMR_FSCHECK 0\n" );
        FPR( stderr, "\n" );

        return( FILE_SIZE_FAIL );
       }

/***/
/* If the predicted and available sizes can be expressed as ints:
/*   Compare the int sizes to see if they are consistant.
/* Otherwise:
/*   Compare the double-precision float sizes instead.
/***/

    if (pInfo.status == FILE_STATUS_INT_OK && aInfo.status == FILE_STATUS_INT_OK)
       {
        if (pInfo.iTotalBytes > aInfo.iTotalBytes)
           {
            if (verbose)
               {
                FPR( stderr, "\n" );
                FPR( stderr, "File Size Check Error:" );
                FPR( stderr, " Insufficient Space for %s\n", path );
                FPR( stderr, "Predicted Size: %s bytes.\n", getNMRIntStr( pInfo.iTotalBytes ));
                FPR( stderr, "Available Size: %s bytes.\n", getNMRIntStr( aInfo.iTotalBytes ));
                FPR( stderr, "\n" );
                FPR( stderr, "To suppress this test, use environment variable settings:\n" );
                FPR( stderr, "\n" );
                FPR( stderr, "   setenv NMR_FSCHECK 0\n" );
                FPR( stderr, "\n" );
               }

            return( FILE_SIZE_TOO_LARGE );
           }
       }
    else
       {
        if (pInfo.dTotalBytes > aInfo.dTotalBytes)
           {
            if (verbose)
               {
                FPR( stderr, "\n" );
                FPR( stderr, "File Size Check Error:" );
                FPR( stderr, " Insufficient Space for %s\n", path );
                FPR( stderr, "Predicted Size: %.0lf bytes.\n", pInfo.dTotalBytes );
                FPR( stderr, "Available Size: %.0lf bytes.\n", aInfo.dTotalBytes );
                FPR( stderr, "\n" );
                FPR( stderr, "To suppress this test, use environment variable settings:\n" );
                FPR( stderr, "\n" );
                FPR( stderr, "   setenv NMR_FSCHECK 0\n" );
                FPR( stderr, "\n" );
               }

            return( FILE_SIZE_TOO_LARGE );
           }
       }

    return( FILE_SIZE_OK );
}

/***/
/* getNMRBytes: returns expected size of NMRPipe data based on header.
/***/

double getNMRBytes( fdata, multiFile, fInfo )

   struct FileSize *fInfo;
   float  fdata[FDATASIZE];
   int    multiFile;
{
    NMR_INT xSize, ySize, zSize, aSize, pipeFlag, cubeFlag, fileCount, dimCount, quadState;
    
    (void) fileSizeInit( fInfo );

    xSize  = (NMR_INT) getParmI( fdata, NDSIZE, CUR_XDIM );
    ySize  = (NMR_INT) getParmI( fdata, NDSIZE, CUR_YDIM );
    zSize  = (NMR_INT) getParmI( fdata, NDSIZE, CUR_ZDIM );
    aSize  = (NMR_INT) getParmI( fdata, NDSIZE, CUR_ADIM );

    fileCount  = (NMR_INT) getParmI( fdata, FDFILECOUNT, 0 );
    dimCount   = (NMR_INT) getParmI( fdata, FDDIMCOUNT,  0 );

    pipeFlag   = getParmI( fdata, FDPIPEFLAG, 0 );
    cubeFlag   = getParmI( fdata, FDCUBEFLAG, 0 );

    if (dimCount < 4) aSize = 1;
    if (dimCount < 3) zSize = 1;
    if (dimCount < 2) ySize = 1;

    if (xSize < 1) xSize = 1;
    if (ySize < 1) ySize = 1;
    if (zSize < 1) zSize = 1;
    if (aSize < 1) aSize = 1;

    if (!multiFile || fileCount < 1) fileCount = 1;

    if (1 == getParmI( fdata, FDQUADFLAG, 0 ))
       quadState = 1;
    else
       quadState = 2;


    if (pipeFlag)
       {
        fInfo->dTotalBytes = (double)sizeof(float)*(FDATASIZE + (double)xSize*ySize*zSize*aSize*quadState);
        fInfo->iTotalBytes = (NMR_INT)sizeof(float)*(FDATASIZE + (NMR_INT)xSize*ySize*zSize*aSize*quadState);
       }
    else
       {
        if (dimCount == 4)
           {
            if (cubeFlag)
               {
                fInfo->dTotalBytes = (double)fileCount*(sizeof(float)*(FDATASIZE + (double)xSize*ySize*zSize*quadState));
                fInfo->iTotalBytes = (NMR_INT)fileCount*(sizeof(float)*(FDATASIZE + xSize*ySize*zSize*quadState));
               }
            else
               {
                fInfo->dTotalBytes = (double)fileCount*(sizeof(float)*(FDATASIZE + (double)xSize*ySize*quadState));
                fInfo->iTotalBytes = (NMR_INT)fileCount*(sizeof(float)*(FDATASIZE + xSize*ySize*quadState));
               }
           }
        else 
           {
            fInfo->dTotalBytes = (double)fileCount*(sizeof(float)*(FDATASIZE + (double)xSize*ySize*quadState));
            fInfo->iTotalBytes = (NMR_INT)fileCount*(sizeof(float)*(FDATASIZE + (NMR_INT)xSize*ySize*quadState));
           }
       }

    if (fInfo->dTotalBytes > NMR_INT_MAX || fInfo->dTotalBytes < 0.0)
       {
        fInfo->blockCount = fInfo->dTotalBytes/fInfo->blockSize;
        fInfo->extraBytes = fInfo->dTotalBytes - fInfo->blockSize*fInfo->blockCount;
        fInfo->status     = FILE_STATUS_INT_TOO_SMALL;
       }
    else
       {
        fInfo->blockCount = fInfo->iTotalBytes/fInfo->blockSize;
        fInfo->extraBytes = fInfo->iTotalBytes - fInfo->blockSize*fInfo->blockCount;
       }

    return( fInfo->dTotalBytes );
}

/***/
/* getFreeBytes: returns space available on filesystem associated with path.
/***/

double getFreeBytes( path, fInfo )

   struct FileSize *fInfo; 
   char   *path;
{
    struct STATFS buff;

    errno = 0;

    (void) fileSizeInit( fInfo );

#define USE_STATFS

#ifdef SGI
    if (statfs( path, &buff, sizeof( struct statfs ), 0 ))
       {
        if (errno) (void) perror( "File Status Notice" );

        fInfo->status = FILE_STATUS_FAIL;
        return( (double) -1.0 );
       }
    else
       {
        fInfo->blockSize  = buff.f_bsize;
        fInfo->blockCount = buff.f_bfree;
       }

#undef USE_STATFS
#endif

#ifdef CRAY
    if (statfs( path, &buff, sizeof( struct statfs ), 0 ))
       {
        if (errno) (void) perror( "File Status Notice" );

        fInfo->status = FILE_STATUS_FAIL;
        return( (double) -1.0 );
       }
    else
       {
        fInfo->blockSize  = buff.f_bsize;
        fInfo->blockCount = buff.f_bfree;
       }

#undef USE_STATFS
#endif

#ifdef ALPHA
    if (statfs( path, &buff, sizeof( buff ) ))
       {
        if (errno) (void) perror( "File Status Notice" );

        fInfo->status = FILE_STATUS_FAIL;
        return( (double) -1.0 );
       }
    else
       {
        fInfo->blockSize  = buff.f_bsize;
        fInfo->blockCount = buff.f_bavail;
       }

#undef USE_STATFS
#endif

#ifdef USE_STATFS
    if (STATFS( path, &buff ))
       {
        if (errno) (void) perror( "File Status Notice" );

        fInfo->status = FILE_STATUS_FAIL;
        return( (double) -1.0 );
       }
    else
       {
        fInfo->blockSize  = buff.f_bsize;
        fInfo->blockCount = buff.f_bavail;
       }
#endif

    fInfo->dTotalBytes = (double)fInfo->blockSize*fInfo->blockCount;
    fInfo->iTotalBytes = (NMR_INT)fInfo->blockSize*fInfo->blockCount;

    if (fInfo->dTotalBytes > NMR_INT_MAX || fInfo->dTotalBytes < 0.0)
       {
        fInfo->status = FILE_STATUS_INT_TOO_SMALL;
       }

    return( fInfo->dTotalBytes );
}

/***/
/* getFileBytes: returns space used by file associated with path.
/*               Returns negative number if test fails.
/***/

double getFileBytes( path, fInfo )

   struct FileSize *fInfo;
   char   *path;
{
    struct stat buff;

    errno = 0;

    (void) fileSizeInit( fInfo );

    if (safe_stat( path, &buff ))
       {
        if (errno) (void) perror( "File Status Notice" );

        fInfo->status = FILE_STATUS_FAIL;
        return( (double) -1.0 );
       }
    else
       {
        fInfo->dTotalBytes = (double)buff.st_size;
        fInfo->iTotalBytes = buff.st_size;
       }

    if (fInfo->dTotalBytes > NMR_INT_MAX)
       {
        fInfo->blockCount = fInfo->dTotalBytes/fInfo->blockSize;
        fInfo->extraBytes = fInfo->dTotalBytes - fInfo->blockSize*fInfo->blockCount;
        fInfo->status     = FILE_STATUS_INT_TOO_SMALL;
       }
    else
       {
        fInfo->blockCount = buff.st_size/fInfo->blockSize;
        fInfo->extraBytes = buff.st_size - fInfo->blockSize*fInfo->blockCount;
       }

    return( fInfo->dTotalBytes );
}

/***/
/* fileSizeInit: initialize a FileSize structure.
/***/

int fileSizeInit( fInfo )

   struct FileSize *fInfo;
{
    fInfo->status      = FILE_STATUS_INT_OK;
    fInfo->blockSize   = DEF_BLOCK_SIZE;
    fInfo->blockCount  = 0;
    fInfo->extraBytes  = 0;
    fInfo->iTotalBytes = 0;
    fInfo->dTotalBytes = 0.0;

    return( 0 );
}

/***/
/* Return a string representing an NMR_INT
/***/

char *GetNMRIntStr( i )

   NMR_INT i;
{
   static char str[NAMELEN+1];
   
#ifdef NMR64
   (void) sprintf( str, "%ld", i );
#else
   (void) sprintf( str, "%d", i );
#endif

   return ( str );
}

/***/
/* Bottom.
/***/
