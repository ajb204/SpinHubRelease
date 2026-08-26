/***/
/* text2argv: converts lines of text to an argv-style list of pointers.
/*
/* Nota Bene:
/*   1. Continuation lines marked with "\" are allowed.
/*   2. Continuation characters must have a leading space " \".
/*   2. A given token cannot be split by continuation characters.
/*   3. Argv list is overwritten with each call.
/*   4. Argv pointers point to tokens in lineList, not copies.
/*   5. On return, text in lineList has whitespace replaced by NULLs.
/***/

#include <string.h>
#include <ctype.h>

#define TOK_INC 512

#include "memory.h"

static char **argvBuff   = (char **) NULL;
static int  argvBuffSize = 0;

int text2argv( lineList,  /* List of lines to analyze.                       */
               argcPtr,   /* On return, number of tokens in argv list.       */
               argvPtr,   /* On return, null terminated argv list of tokens. */
               nextLine ) /* On return, points to next line in list.         */

   char **lineList, ***nextLine, ***argvPtr;
   int  *argcPtr;
{
    char *sPtr, lastChar;
    int  argc, error;

    error = 0;
    argc  = 0;

/***/
/* Allocate the local argv list buffer:
/***/

    if (!argvBuffSize)
       {
        argvBuffSize += TOK_INC;
        argvBuff  = (char **) voidAlloc( "argv", sizeof(char *)*argvBuffSize );

        if (!argvBuff)
           {
            argvBuffSize = 0;
            error        = 1;

            goto shutdown;
           }
       }

/***/
/* For each line:
/*   Check for sufficient space in argv array.
/*
/*   Scan for tokens in the line:
/*     Replace spaces with nulls.
/*     Skip end-of-line backslash character.
/*     Increment argc.
/*     Skip past token body.
/*
/* Stop if last line was not a continuation.
/*
/* Null-terminate the list.
/***/

    while( *lineList )
       {
        if (argc >= argvBuffSize - 1)
           {
            error = 2;
            goto shutdown;
           }

        sPtr = *lineList;

        while( *sPtr )
           {
            if (isspace( *sPtr ))
               {
                lastChar = *sPtr;
                *sPtr++  = '\0';
               }
            else if (*sPtr == '\\')
               {
                if (sPtr[1] == '\n' || sPtr[1] == '\0')
                   {
                    lastChar = *sPtr;
                    break;
                   }
                else
                   {
                    lastChar = '\0';
                    sPtr++;
                   }
               }
            else
               {
                if (*sPtr == '\'')
                   {
                    argvBuff[argc++] = ++sPtr;
                    while( *sPtr && *sPtr != '\'') sPtr++;
                    if (*sPtr) *sPtr = ' ';
                   }
                else
                   {
                    argvBuff[argc++] = sPtr;
                    while( *sPtr && !isspace( *sPtr )) sPtr++;
                   }

                lastChar = *sPtr;
               }
           }

        lineList++;

        if (lastChar != '\\') break;
       }

    argvBuff[argc] = (char *) NULL;

shutdown:

    *nextLine = lineList;
    *argvPtr  = argvBuff;
    *argcPtr  = argc;

    return( error );
}

/***/
/* text2argvFree: deallocate local argv buffer.
/***/

int text2argvFree()
{
    (void) deAlloc( "argv", argvBuff, sizeof(char *)*argvBuffSize );

    argvBuff     = (char **) NULL;
    argvBuffSize = 0;

    return( 0 );
}

/***/
/* str2argv: create an argv-style list from text in a single string.
/***/

int str2argvC( baseCmnd, cmnd, argcPtr, argvPtr )

   char *baseCmnd, *cmnd;

   int  *argcPtr;
   char ***argvPtr;
{
    char **thisArgv;
    int  thisArgc, i, n;

    *argcPtr = 0;
    *argvPtr = (char **)NULL;

    if (!cmnd) return( 0 );

    thisArgc = 0;
    thisArgv = (char **)NULL;

    thisArgc = 1 + cntTokenS( cmnd );
    
    if (!(thisArgv = (char **)voidAlloc( "str2argv", sizeof(char *)*(thisArgc + 1) )))
       {
        return( 1 );
       }

    for( i = 0; i < thisArgc + 1; i++ ) thisArgv[i] = (char *)NULL;

    *argcPtr = thisArgc;
    *argvPtr = thisArgv;

    for( i = 0; i < thisArgc; i++ )
       {
        if (i)
           thisArgv[i] = strDup( getTokenS( i, cmnd ));
        else
           thisArgv[i] = strDup( baseCmnd );

        if (!thisArgv[i]) return( 1 );
       }

    return( 0 );
}

int str2argv( cmnd, argcPtr, argvPtr )

   char *cmnd, ***argvPtr;
   int  *argcPtr;
{
    char **thisArgv;
    int  thisArgc, i, n;

    *argcPtr = 0;
    *argvPtr = (char **)NULL;

    if (!cmnd) return( 0 );

    thisArgc = 0;
    thisArgv = (char **)NULL;

    thisArgc = cntTokenS( cmnd );

    if (!(thisArgv = (char **)voidAlloc( "str2argv", sizeof(char *)*(thisArgc + 1) )))
       {
        return( 1 );
       }

    for( i = 0; i < thisArgc + 1; i++ ) thisArgv[i] = (char *)NULL;

    *argcPtr = thisArgc;
 
    for( i = 0; i < thisArgc; i++ )
       {
        if (!(thisArgv[i] = strDup( getTokenS( i+1, cmnd ))))
           {
            return( 1 );
           }
       }

    return( 0 );
}

int str2argvFree( thisArgc, thisArgv )

   int  thisArgc;
   char **thisArgv;
{
    int i;

    if (thisArgc < 0 || !thisArgv) return( 0 );

    for( i = 0; i < thisArgc; i++ )
       {
        (void) strFree( thisArgv[i] );
        thisArgv[i] = (char *)NULL;
       }

    (void) deAlloc( "str2argv", thisArgv, sizeof(char *)*(thisArgc + 1) );

    return( 0 );
}

