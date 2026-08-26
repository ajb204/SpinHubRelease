/******************************************************************************/
/*                                                                            */
/*                   ---- NIH NMR Software System ----                        */
/*                         Copyright 1992 - 1996                              */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/

/***/
/* stralloc: modules for string manipulation which use dynamic allocation,
/*           but only when needed.
/***/

#include <string.h>
#include <stdio.h>

#include "stralloc.h"
#include "memory.h"
#include "token.h"
#include "prec.h"

#define FPR (void)fprintf

static char *ctemp       = (char *)NULL;
static int  ctempBuffLen = 0; 

static char ctempBuff[NAMELEN+1];

/***/
/* strcatAlloc: dynamic version of strcat() library function:
/*
/* strcat( dest, src );
/*
/*  Becomes:
/*
/* strcatAlloc( &dest, src, cBuff, cBuffLen, &buffLen );
/*
/*  Where:
/*
/*   dest:     is now a pointer which can be changed;
/*   cBuff:    is a static or automatic char array to use as an initial buffer;
/*   cBuffLen: is the size of the cBuff array;
/*   buffLen:  is the current size of dest.
/*
/*  Initially:
/*
/*    char *dest, cBuff[NAMELEN+1];
/*    int  cBuffLen, buffLen;
/*
/*    dest     = cBuff;
/*    cBuff[0] = NULL;
/*    cBuffLen = NAMELEN+1;
/*    buffLen  = cBuffLen;
/***/

char *strcatAlloc( destPtr, src, cBuff, cBuffLen, maxLength )

   char **destPtr, *src, *cBuff;
   int  cBuffLen, *maxLength;
{
    int  length, newLength;
    char *dest;

/***/
/* Stop if there is no pointer for return result.
/* Initialize pointer if needed.
/***/

    if (!destPtr)
       {
        FPR( stderr, "StrCatAlloc Error: Null pointer.\n" );
        return( (char *)NULL );
       }

    if (!*destPtr)
       {
        *maxLength = cBuffLen;
        *destPtr   = cBuff;
       }

    dest = *destPtr;

    if (!src) return( dest );

/***/
/* If current buffer is too small:
/*  Increment buffer size by a fixed length (NAMELEN) if this is enough.
/*  Otherwise, increment buffer to required length plus some extra.
/***/

    length = istrlen( src ) + istrlen( dest ) + 2;

    if (length >= *maxLength)
       {
        if (length >= *maxLength + NAMELEN)
           newLength = length + NAMELEN + 1;
        else
           newLength = *maxLength + NAMELEN + 1;

        if (dest == cBuff)
           {
            if (!(dest = charAlloc( "str", newLength )))
               {
                FPR( stderr, "StrCatAlloc Error: alloc failed.\n" );
                return( (char *)NULL );
               }

            (void) strcpy( dest, cBuff );
           }
        else
           {
            if (!(dest = charReAlloc( "str", dest, newLength, *maxLength )))
               {
                FPR( stderr, "StrCatAlloc Error: realloc failed.\n" );
                return( (char *)NULL );
               }
           }

        *destPtr   = dest;
        *maxLength = newLength;
       }

/***/
/* Perform the catenation, return the result:
/***/

    (void) strcat( dest, src );

    return( dest );
}

/***/
/* strpcatAlloc: as strcatAlloc above, but catenates src before
/*               contents of dest, rather than after.
/***/

char *strpcatAlloc( destPtr, src, cBuff, cBuffLen, maxLength )

   char **destPtr, *src, *cBuff;
   int  cBuffLen, *maxLength;
{
    char c, *dest;
    int  srcLoc, destLoc;

    if (!(dest = strcatAlloc( destPtr, src, cBuff, cBuffLen, maxLength )))
       {
        return( (char *)NULL );
       }
 
    destLoc = istrlen( dest ) - 1;
    srcLoc  = destLoc - istrlen( src );

    while( srcLoc >= 0 )
       {
        dest[destLoc--] = dest[srcLoc--];
       }

    destLoc = 0;

    while( *src )
       {
        dest[destLoc++] = *src++;
       }

    return( dest );
}

/***/
/* strcpyAlloc: dynamic version of strcpy, as strcat version above.
/***/

char *strcpyAlloc( destPtr, src, cBuff, cBuffLen, maxLength )

   char **destPtr, *src, *cBuff;
   int  cBuffLen, *maxLength;
{
    char *dest;

/***/
/* Test for null return pointer, initialize it if needed:
/***/

    if (!destPtr)
       {
        FPR( stderr, "StrCpyAlloc Error: Null pointer.\n" );
        return( (char *)NULL );
       }

    if (!*destPtr)
       {
        *maxLength = cBuffLen;
        *destPtr   = cBuff;
       }

/***/
/* Perform the copy by catenating the src string to a null string:
/***/

    dest    = *destPtr;
    dest[0] = '\0';

    (void) strcatAlloc( &dest, src, cBuff, cBuffLen, maxLength );

    if (!dest) return( (char *)NULL ) ; 

    *destPtr = dest;

    return( dest );
}

/***/
/* fgetsAlloc: dynamic version of fgets string reading function.
/***/

char *fgetsAlloc( strPtr, maxLength, inUnit, 
                  inBuff, inBuffLen,
                  cBuff, cBuffLen )

   int  cBuffLen, *maxLength, inBuffLen;
   char **strPtr, *inBuff, *cBuff;
   FILE *inUnit;
{
    char *str, *result;

/***/
/* Test for null return pointer, initialize it if needed:
/***/

    if (!strPtr)
       {
        FPR( stderr, "FGetSAlloc Error: Null pointer.\n" );
        return( (char *)NULL );
       }

    if (inBuffLen < 3)
       {
        FPR( stderr, "FGetSAlloc Error: Buffer too small.\n" );
        return( (char *)NULL );
       }

    if (!*strPtr)
       {
        *maxLength = cBuffLen;
        *strPtr    = cBuff;
       }

    str    = *strPtr;
    str[0] = '\0';

/***/
/* While an end-of-line has not yet been encountered:
/*  Read at most a buffer full of text.
/*  Stop on end-of-file.
/*  Catenate the buffer to the result with dynamic allocation if needed.
/***/

    result = (char *)NULL;

    do
       {
        if (!fgets( inBuff, inBuffLen-1, inUnit )) /* fgetsAlloc not needed. */
           {
            return( result );
           }

        (void) strcatAlloc( &str, inBuff, cBuff, cBuffLen, maxLength );

        result = str;
       }
    while( !strchr( inBuff, '\n' ));

    *strPtr = result;
 
    return( result );
}

/***/
/* strAllocInit: set initial values for dynamic string functions.
/***/

int strAllocInit( cPtr, cBuff, cBuffLen, maxLength )

   int  cBuffLen, *maxLength;
   char **cPtr, *cBuff;
{
    if (!cPtr) return( 0 );

    *maxLength = cBuffLen;
    *cPtr      = cBuff;
    cBuff[0]   = '\0';
   
    return( 0 );
}

/***/
/* strAllocFree: deallocate string buffer if needed.
/***/

int strAllocFree( cPtr, cBuff, cBuffLen, maxLength )

   int  cBuffLen, *maxLength;
   char **cPtr, *cBuff;
{
    if (!cPtr) return( 0 );

    if (*cPtr != cBuff)
       {
        (void) deAlloc( "str", *cPtr, sizeof(char)*(*maxLength) );
       }

    *maxLength = cBuffLen;
    *cPtr      = cBuff;
    
    return( 0 );
}

/***/
/* ctempAlloc: provide a work buffer at least n chars long.
/***/

char *ctempAlloc( n )

   int n;
{
    if (!ctemp)
       {
        ctemp        = ctempBuff;
        ctempBuffLen = NAMELEN+1;
       }

    if (n > ctempBuffLen)
       {
        FPR( stderr, "Expanding to %d\n", n );

        (void) ctempFree();

        ctempBuffLen = n;

        if (!(ctemp = charAlloc( "str", ctempBuffLen+1 )))
           {
            FPR( stderr, "CTempAlloc memory allocation error.\n" );

            ctemp        = ctempBuff;
            ctempBuffLen = NAMELEN + 1;
           } 
       }

    ctemp[0] = '\0';

    return( ctemp );
}

/***/
/* ctempFree: deallocate work buffer.
/***/

int ctempFree()
{
   if (ctemp != ctempBuff)
      {
       (void) deAlloc( "str", ctemp, sizeof(char)*(ctempBuffLen+1) ); 
      }

   ctemp        = (char *)NULL;
   ctempBuffLen = 0;
   ctempBuff[0] = '\0';

   return( 0 );
}
 
/***/
/* Bottom.
/***/
