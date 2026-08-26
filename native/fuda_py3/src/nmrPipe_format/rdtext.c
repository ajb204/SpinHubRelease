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

#include <stdio.h>
#include <string.h>

#include "stralloc.h"
#include "memory.h"
#include "token.h"
#include "prec.h"
#include "rdtext.h"

#define MEM_INC 256

/***/
/* rdText: read the contents of text-file named inName.
/*         Returns the buffer allocated for the text data.
/*
/*         Use freeText to free the returned buffer.
/***/

char **rdText( inName, lineCount, error )

  char *inName;
  int  *lineCount, *error;
{
    char  cBuff[NAMELEN+1], inBuff[NAMELEN+1];
    char  **buffer, **bPtr, *cPtr, *line;
    int   length, inLen, cBuffLen, maxLen;
    FILE  *inUnit;

/***/
/* Initialize:
/***/

    *error     = 0;
    *lineCount = 0;

    buffer    = (char **) NULL;
    inUnit    = (FILE *) NULL;

    line      = cBuff;
    inLen     = NAMELEN+1;
    cBuffLen  = NAMELEN+1;
    maxLen    = NAMELEN+1;

    cBuff[0]  = '\0';
    inBuff[0] = '\0';
 
    if (!(inUnit = fopen( inName, "r" )))
       {
        *error = -1;
        return( buffer );
       }

/***/
/* Count lines in file.
/* Allocate memory for string pointers.
/***/

    while( fgetsAlloc( &line, &maxLen, inUnit, inBuff, inLen, cBuff, cBuffLen ))
       {
        (*lineCount)++;
       }

    buffer = (char **) voidAlloc( "rdText", sizeof(char *)*(1 + *lineCount) );

    if (!buffer)
       {
        *error = -1;
        return( buffer );
       }

/***/
/* Rewind file and read lines.
/* For each line, remove trailing newline, allocate memory and save.
/***/

    bPtr = buffer;

    (void) rewind( inUnit );

    while( fgetsAlloc( &line, &maxLen, inUnit, inBuff, inLen, cBuff, cBuffLen ))
       {
        length = istrlen( line );

        if (length > 0) line[length-1] = '\0';

        if (!(cPtr = strDup( line )))
           {
            *error = -1;

            (void) strAllocFree( &line, cBuff, cBuffLen, &maxLen );
            (void) fclose( inUnit );

            return( buffer );
           }

        *(bPtr++) = cPtr; 
       }

/***/
/* Add trailing null pointer to list of string pointers.
/***/

    *bPtr = (char *) NULL;

    (void) strAllocFree( &line, cBuff, cBuffLen, &maxLen );
    (void) fclose( inUnit );

    return( buffer );
}

/***/
/* rdText2: read lines of text from non-seeking text source at inUnit.
/***/

char **rdText2( inUnit, lineCount, endOfText, verboseFlag, error )

  FILE *inUnit;
  char *endOfText;
  int  *lineCount, verboseFlag, *error;
{
    int  lineLength, buffLength, oldLength, maxLen, cBuffLen, inLen;
    char **buffer, *line, cBuff[NAMELEN+1], inBuff[NAMELEN+1];

/***/
/* Initialize:
/***/

    *error     = 0;
    *lineCount = 0;

    buffer     = (char **) NULL;
    buffLength = 0;

    line       = cBuff;
    inLen      = NAMELEN+1;
    cBuffLen   = NAMELEN+1;
    maxLen     = NAMELEN+1;

    cBuff[0]   = '\0';
    inBuff[0]  = '\0';

/***/
/* For each line:
/*   Display line if needed.
/*   Remove trailing newline.
/*   Expand string list as needed.
/*   Allocate and copy new string.
/*   Stop if string is endOfText marker.
/***/

    while( fgetsAlloc( &line, &maxLen, inUnit, inBuff, inLen, cBuff, cBuffLen ))
       {
        if (verboseFlag) (void) fprintf( stderr, "%s", line );

        lineLength = istrlen( line );

        if (lineLength > 0) line[lineLength-1] = '\0';

        if (*lineCount + 2 > buffLength)
           {
            oldLength   = buffLength;
            buffLength += MEM_INC;

            buffer = (char **) voidReAlloc( "rdText", buffer, sizeof(char *)*buffLength, sizeof(char *)*oldLength );
     
            if (!buffer)
               {
                *error = -1;
                return( buffer );
               }
           }

#ifndef S_SPLINT_S
        if (!(buffer[*lineCount] = strDup( line ))) 
           {
            *error = -1;
            return( buffer );
           }
#endif
        (*lineCount)++;

        if (endOfText && !strcasecmp( line, endOfText )) break;
       }

/***/
/* Add trailing null pointer to list of string pointers.
/* Contract string list to smallest size needed.
/***/

    if (buffer)
       {
        buffer[*lineCount] = (char *) NULL;
        buffer = (char **) voidReAlloc( "rdText", buffer, sizeof(char *)*(*lineCount + 1), sizeof(char *)*buffLength );
       }

    (void) strAllocFree( &line, cBuff, cBuffLen, &maxLen );

    return( buffer );
}

/***/
/* rdTextC: read a text file into a single character array.
/***/

char *rdTextC( inName, lineCount, length, error )

  char *inName;
  int  *lineCount, *length, *error;
{
    char  cBuff[NAMELEN+1], inBuff[NAMELEN+1];
    char  *buffer, *sPtr, *line;

    int   i, inLen, cBuffLen, maxLen;
    FILE  *inUnit;

/***/
/* Initialize:
/***/

    *error     = 0;
    *lineCount = 0;
    *length    = 0;

    buffer    = (char *) NULL;
    inUnit    = (FILE *) NULL;

    line      = cBuff;
    inLen     = NAMELEN+1;
    cBuffLen  = NAMELEN+1;
    maxLen    = NAMELEN+1;

    cBuff[0]  = '\0';
    inBuff[0] = '\0';
 
    if (!(inUnit = fopen( inName, "r" )))
       {
        *error = -1;
        return( buffer );
       }

/***/
/* Count lines in file and tally length of file text.
/* Allocate the text buffer.
/***/

    while( fgetsAlloc( &line, &maxLen, inUnit, inBuff, inLen, cBuff, cBuffLen ))
       {
        *length += istrlen( line ) + 1;
        (*lineCount)++;
       }

    (*length)++;

    if (!(buffer = charAlloc( "rdText", *length )))
       {
        *error = -1;
        return( buffer );
       }

    for( i = 0; i < *length; i++ ) buffer[i] = '\0'; 

/***/
/* Rewind the file and read lines into the buffer.
/***/

    sPtr = buffer;

    (void) rewind( inUnit );

    while( fgetsAlloc( &line, &maxLen, inUnit, inBuff, inLen, cBuff, cBuffLen ))
       {
        (void) strcpy( sPtr, line );
        sPtr += istrlen( line ) + 1; 
       }

    (void) strAllocFree( &line, cBuff, cBuffLen, &maxLen );
    (void) fclose( inUnit );

    return( buffer );
}

/***/
/***/
/* cntText: count lines in a text file.
/***/

int cntText( inName, lineCount )

  char *inName;
  int  *lineCount;

{
    char  line[NAMELEN+1];
    FILE  *inUnit;

/***/
/* Initialize:
/***/

    *lineCount = 0;
    inUnit     = (FILE *) NULL;

    if (!(inUnit = fopen( inName, "r" )))
       {
        return( 1 );
       }

    while( fgets( line, NAMELEN, inUnit ))   /* fgetsAlloc not needed. */
       {
        if (strchr( line, '\n' )) (*lineCount)++;
       }

    (void) fclose( inUnit );

    return( 0 );
}

/***/
/* buildList: create space for a text list with "length" strings.
/***/

char **buildList( length )

   int length;
{
    char **buffer, **cPtr;

    if (!(buffer = (char **) voidAlloc( "rdText", sizeof(char *)*(1+length) )))
       {
        return( buffer );
       }

    cPtr = buffer;

    while( length-- )
       {
        if (!(*cPtr = charAlloc( "rdText", NAMELEN+1 )))
           {
#ifndef S_SPLINT_S
            return( cPtr );
#endif
           }

        cPtr++;
       }

    *cPtr = (char *) NULL;

    return( buffer );
}

/***/
/* freeText: free text buffer as allocated by rdText.
/***/

void freeText( buffer )

   char **buffer;
{
    char **strPtr;
    int  count;

    if (!buffer) return;

    strPtr = buffer;
    count  = 0;

    while( *strPtr )
       {
        (void) deAlloc( "strDup", *strPtr, istrlen(*strPtr) + 1 );
        strPtr++;
        count++;
       }

    (void) deAlloc( "rdText", buffer, sizeof(char *)*(count + 1) );
}

/***/
/* Bottom.
/***/
