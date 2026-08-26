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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "token.h"
#include "namelist.h"
#include "cmndargs.h"
#include "atof.h"
#include "prec.h"

#define FPR (void) fprintf

/***/
/* intStrArgD: interpret command line via Name/Value table.
/*             Returns negative number on error; see "namelist.h" for
/*             summary of return codes.
/***/

int intStrArgD( argc, argv, flagName, value, name2Val )

   int    argc;
   int    *value;
   char   **argv, *flagName;
   struct NameVal *name2Val;

{
    int status;
    float rtemp;

    rtemp  = *value;
    status = fltStrArgD( argc, argv, flagName, &rtemp, name2Val );
    *value = rtemp;

    return( status );
}

/***/
/* fltStrArgD: interpret command line via Name/Value table;
/*             Returns negative number on error; see "namelist.h" for
/*             summary of return codes.
/***/

int fltStrArgD( argc, argv, flagName, value, name2Val )

   int    argc;
   float  *value;
   char   **argv, *flagName;
   struct NameVal *name2Val;

{
    int  status, loc;
    char *aPtr;

    status = KEYWORD_NONE;

    if ((loc = flagLoc( argc, argv, flagName )))
       {
        if (loc >= argc - 1)
           {
            FPR( stderr, "Keyword Error: bad or missing value for" );
            FPR( stderr, " '%s'\n", flagName );
            return( KEYWORD_BAD );
           }

        aPtr   = getNthArg( argc, argv, loc+1 );
        status = getValByName( flagName, aPtr, name2Val, value );
       }

    return( status );
}

/***/
/* strStrArgD: interpret command line via Name/Value table;
/*             Returns negative number on error; see "namelist.h" for
/*             summary of return codes.
/***/

int strStrArgD( argc, argv, flagName, value, name2Val )

   int    argc;
   char   **value;
   char   **argv, *flagName;
   struct NameVal *name2Val;
{
    int  status, i, loc;
    char *aPtr;

    status = KEYWORD_NONE;

    if ((loc = flagLoc( argc, argv, flagName )))
       {
        if (loc >= argc - 1)
           {
            FPR( stderr, "Keyword Error: bad or missing value for" );
            FPR( stderr, " '%s'\n", flagName );
            return( KEYWORD_BAD );
           }

        aPtr = getNthArg( argc, argv, loc+1 );

        for( i = 0; name2Val[i].name; i++ )
           {
            if (!strcasecmp( name2Val[i].name, aPtr ))
               {
                *value = name2Val[i].name;
                return( loc );
               }
           }
       }

    return( status );
}

/***/
/* getNameByVal: find a pointer to a name in the given namelist which
/*               corresponds to a given value.  If more than one value
/*               matches, the last match is returned.
/*
/*               Returns the number of values which matched.
/***/

int getNameByVal( thisNamePtr, name2Val, value )

   struct NameVal *name2Val;
   char   **thisNamePtr;
   float  value;
{
    int i, count;

#ifndef S_SPLINT_S
    *thisNamePtr = (char *)NULL;
#endif

    count        = 0;

    for( i = 0; name2Val[i].name; i++ )
       {
        if (value == name2Val[i].val)
           {
            *thisNamePtr = name2Val[i].name;
            count++;
           }
       }

    return( count );
}

char *getNameStrByVal( name2Val, value )

   struct NameVal *name2Val;
   float  value;
{
    static char s[] = "???";
    int i;

    for( i = 0; name2Val[i].name; i++ )
       {
        if (value == name2Val[i].val) return( name2Val[i].name );
       }

    return( s );
}

/***/
/* getValByName: find the value in a namelist corresponding to a given name.
/***/

int getValByName( flagName, thisArg, name2Val, value )

   char   *flagName, *thisArg;
   struct NameVal *name2Val;
   float  *value;
{
    int i;

    if (!thisArg) return( KEYWORD_BAD );

    if (isalpha( thisArg[0] ))
       {
        for( i = 0; name2Val[i].name; i++ )
           {
            if (!strcasecmp( name2Val[i].name, thisArg ))
               {
                *value = name2Val[i].val;
                return( KEYWORD_OK );
               }
           }

        FPR( stderr, "Keyword Error: bad value for" );
        FPR( stderr, " '%s'; Valid Options:\n\n", flagName );
        FPR( stderr, "%-14s %s\n\n", "Code Name", "Numerical Value" ); 

        for( i = 0; name2Val[i].name; i++ )
           {
            FPR( stderr,
                 " %-14s % 2.0f\n",
                 name2Val[i].name, name2Val[i].val );
           }

        return( KEYWORD_BAD );
       }
    else
       {
        *value = ATOF( thisArg );

        return( KEYWORD_NUM );
       }
}

/***/
/* Bottom.
/***/
