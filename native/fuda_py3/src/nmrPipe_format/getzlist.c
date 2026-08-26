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
/* getZList: utilities to read or build a systematic file list.
/***/

#include <stdio.h>
#include <string.h>

#include "cmndargs.h"
#include "memory.h"
#include "rdtext.h"
#include "fdatap.h"
#include "prec.h"

#define MAXZPLANE 256 
#define MAXAPLANE 128 

int fileExists(); //AJB 14/09/22

/***/
/* getZList: read or build a Z or AZ format file list.
/***/

int getZList( listName, nameIsTemplate, nameCount, zCount, aCount, listPtr )

   int  nameIsTemplate, nameCount, zCount, aCount; 
   char *listName, ***listPtr;
{
    int  iz, ia, lineCount, useAZFmt, error;
    char ctemp[NAMELEN+1], **sPtr;

    error = 0;

/***/
/* If the name is a template:
/*  Allocate space for the lines of text.
/*  For each plane:
/*   Create a name for the plane.
/*   Add it to the list.
/*
/* If the name is a file list:
/*  Read the list of plane names from the file.
/***/

    if (nameIsTemplate)
       {
        *listPtr = (char **) voidAlloc( "get", sizeof(char *)*(1 + nameCount) );
  
        if (!*listPtr)
           {
            (void) fprintf( stderr, "GetZList Error allocating memory.\n" );
            return( 1 );
           }

        sPtr            = *listPtr;
        sPtr[nameCount] = (char *) NULL;
        useAZFmt        = isAZFmt( listName );
 
        for( ia = 0; ia < aCount; ia++ )
           { 
            for( iz = 0; iz < zCount; iz++ )
               {
                if (useAZFmt)
                   (void) sprintf( ctemp, listName, ia + 1, iz + 1 );
                else
                   (void) sprintf( ctemp, listName, iz + ia*zCount + 1 );

                if (!(*sPtr = charAlloc( "get", strlen(ctemp)+1 ))) 
                   {
#ifndef S_SPLINT_S
                    return( 1 );
#endif
                   }

                (void) strcpy( *sPtr, ctemp );

                sPtr++;
               }
           }
       }
    else
       {
        *listPtr = rdText( listName, &lineCount, &error );

        if (!*listPtr)
           {
            (void) fprintf( stderr, "GetZList Error allocating memory.\n" );
            return( 1 );
           }

        if (lineCount != nameCount)
           {
            (void) fprintf( stderr, "GetZList Error: List Size Mismatch.\n" );
            return( 1 );
           }
       }

#ifdef S_SPLINT_S
   (void) free( (void *)sPtr );
#endif
 
    return( error );
}

/***/
/* getZListN: read or build a Z or AZ format file list structure with plane/cube offsets.
/***/

//int getZListN( char *listName,int  zCount, int aCount, int dimCount, struct ZListInfo ***listPtr );

int getZListN( listName, zCount, aCount, dimCount, listPtr )

   char   *listName;
   int    zCount, aCount; 
int dimCount;
struct ZListInfo ***listPtr;
{
    int    i, n, iz, ia, fmtCount;
    struct ZListInfo **zPtr, *itemPtr;
    char   ctemp[NAMELEN+1], *sPtr;

    if (zCount < 1) zCount = 1;
    if (aCount < 1) aCount = 1;

    if (dimCount < 4) aCount = 1;

    n = 1 + zCount*aCount;
   
    if (!(*listPtr = (struct ZListInfo **)voidAlloc( "get", sizeof(struct ZListInfo *)*n )))
       {
        (void) fprintf( stderr, "GetZListN Error Allocating List.\n" );
        return( 1 );
       }

    zPtr = *listPtr;

    for( i = 0; i < n; i++ ) 
       {
        zPtr[i] = (struct ZListInfo *)NULL;
       }

    for( i = 0; i < n; i++ ) 
       {
        if (!(itemPtr = (struct ZListInfo *)voidAlloc( "get", sizeof(struct ZListInfo))))
           {
            (void) fprintf( stderr, "GetZListN Error Allocating List Item.\n" );
            return( 1 );
           }

        zPtr[i] = itemPtr;

        itemPtr->name = (char *)NULL;
        itemPtr->zOff = 0;
        itemPtr->aOff = 0;
       }

    fmtCount = getFmtCount( listName );
    zPtr     = *listPtr;
    i        = 0;
 
    for( ia = 0; ia < aCount; ia++ )
       { 
        for( iz = 0; iz < zCount; iz++ )
           {
            itemPtr = zPtr[i++];

            if (fmtCount == 2)
               (void) sprintf( ctemp, listName, ia + 1, iz + 1 );
            else if (fmtCount == 1)
               (void) sprintf( ctemp, listName, dimCount == 4 ? ia + 1 : iz + 1 );
            else
               (void) strcpy( ctemp, listName );

            if (!(sPtr = strDup( ctemp )))
               {
                (void) fprintf( stderr, "GetZListN Error Allocating Filename.\n" );
                return( 1 );
               }

            itemPtr->name = sPtr;
            itemPtr->zOff = iz;
            itemPtr->aOff = ia;
           }
       }

    return( 0 );
}

/***/
/* freeZListN: free a Z or AZ format file list.
/***/

int freeZListN( zCount, aCount, dimCount, zList )

   int    zCount, aCount;
int dimCount;
   struct ZListInfo **zList;
{
    int    i, n, iz, ia, error;
    struct ZListInfo *itemPtr;

    if (!zList) return( 0 );

    if (zCount < 1) zCount = 1;
    if (aCount < 1) aCount = 1;

    if (dimCount < 4) aCount = 1;

    n = 1 + zCount*aCount;

    for( i = 0; i < n; i++ )
       {
        itemPtr = zList[i];
        if (!itemPtr) break;

        if (itemPtr->name) (void) strFree( itemPtr->name );
        (void) deAlloc( "get", itemPtr, sizeof(struct ZListInfo) );
       }

    (void) deAlloc( "get", zList, sizeof(struct ZListInfo *)*n );

    return( 0 );
}

/***/
/* getZCount: finds file count in 2D series.
/***/

int getZCount( listName, nameIsTemplate, nameCount, zCount, aCount, parFlag )

   int  nameIsTemplate, *nameCount, *zCount, *aCount, parFlag;
   char *listName;
{
    int   iz, ia, zSize, aSize, useAZFmt;
    char  inName[NAMELEN+1];
    float fdata[FDATASIZE];

/***/
/* In Parallel processing mode:
/*  A series of files must be tested for existance.
/* Otherwise:
/*  Only the first file in a series needs to exist.
/***/

    if (parFlag)
       {
        zSize = MAXZPLANE;
        aSize = MAXAPLANE;
       }
    else
       {
        zSize = 1;
        aSize = 1;
       }

/***/
/* If the file name is a template:
/*  Try to find some file in the series, and read its header.
/*  Find the Z and A sizes.
/*
/* If the file name is a text file list:
/*  Count the lines of text in the file.
/***/
 
    if (nameIsTemplate)
       {
        useAZFmt = isAZFmt( listName );

        for( ia = 1; ia <= aSize; ia++ )
           { 
            for( iz = 1; iz <= zSize; iz++ )
               {
                if (useAZFmt)
                   (void) sprintf( inName, listName, ia, iz );
                else
                   (void) sprintf( inName, listName, iz  + (ia-1)*zSize );

                if (fileExists( inName ))
                   {
                    if (rdFDATA( inName, fdata ))
                       {
                        (void) fprintf( stderr,
                                        "GetZCnt Error, file %s\n",
                                        inName );

                        return( 1 );
                       }

                    (void) fixfdata( fdata );

                    *zCount = getParm( fdata, NDSIZE, CUR_ZDIM );
                    *aCount = getParm( fdata, NDSIZE, CUR_ADIM );

                    if (*zCount < 1) *zCount = 1;
                    if (*aCount < 1) *aCount = 1;

                    *nameCount = (*zCount)*(*aCount);

                    return( 0 );
                   }
               }
           }

        (void) fprintf( stderr, "GetZCnt Error: no input files found.\n" );
        return( 1 );
       }
    else
       {
        if (cntText( listName, nameCount ))
           {
            (void) fprintf( stderr, "GetZCnt Error in file %s\n", listName );
            return( 1 );
           }

        *zCount = *nameCount;
        *aCount = 1;
       }

    return( 0 );
}

/***/
/* showZList: displays name list, as for debugging.
/***/

int showZList( nameList, nameCount )

   char **nameList;
   int  nameCount;
{
    int i;

    (void) fprintf( stderr, "\n" );

    for( i = 1; i <= nameCount; i++ )
       {
        (void) fprintf( stderr, "%5d. %s\n", i, *nameList );
        nameList++;
       }

    (void) fprintf( stderr, "\n" );

    return( 0 );
}
