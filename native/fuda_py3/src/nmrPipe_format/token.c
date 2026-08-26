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
/* token: tools for tokenizing blank-delimited tokens.
/*        Defined in "token.h", along with global
/*        variable "badToken"; badToken will be set
/*        to zero after a successful call to any
/*        of these procedures; or non-zero otheriwse.
/*        Non-zero "badToken" indicates that the
/*        returned value is undefined.
/*
/* In the procedures listed, "string" is the character
/* string to be tokenized, "n" is token number to extract,
/* starting at 1.
/*
/* The following extract the Nth token from string as
/* a string, float, and integer, respectively. The string versions
/* store their return value in a static space which will be over-written
/* by other calls to these token utilities.
/*
/*        char *getToken( n, string ) 
/*        char *getTokenS( n, string ) (for quoted strings)
/*        float getTokenR( n, string )
/*        int   getTokenI( n, string )
/*
/* The following also extract the Nth token
/*
/* The following extract the leftmost and rightmost tokens
/* in string, compressing multiple spaces, and starting from 
/* the Nth token:
/*
/*        char *leftToken( n, string )
/*        char *rightToken( n, string )
/*
/* The following returns the number of tokens in string:
/*
/*        int   cntToken( string )
/*        int   cntTokenS( string ) (for quoted strings).
/*
/* This function returns the token number of key if it is found
/* in string:
/*
/*        int   findToken( string, key )
/*
/*
/* This function returns the number of characters in the longest token:
/*
/*        int   getTokenMaxLength( string )
/*
/* NOTA BENE: char * token functions above will return
/*            pointers to character arrays that will be
/*            overwritten in successive calls; return
/*            strings should be copied as needed.
/***/

/***/
/* For experimentation, compile this stand-alone test:
/***/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "stralloc.h"
#include "memory.h"
#include "minmax.h"
#include "token.h"
#include "prec.h"

static char *lValue = (char *)NULL;
static char *rValue = (char *)NULL;
static char *tValue = (char *)NULL;
static char *cValue = (char *)NULL;

static char lBuff[NAMELEN+1];
static char rBuff[NAMELEN+1];
static char tBuff[NAMELEN+1];
static char cBuff[NAMELEN+1];

static int  lBuffLen = 0, rBuffLen = 0, tBuffLen = 0, cBuffLen = 0;

#define FPR (void)fprintf
#define PR  (void)printf

#ifdef STANDALONE

int main( argc, argv )

   int  argc;
   char **argv;
{
    char string[NAMELEN+1], *token;
    int  i, loc1, locN, count, n, ival;
    float rval;

    loc1 = -1;
    locN = -1;

    PR( "Line: " );
    (void) gets( string );
    PR( "Token: " );
    (void) scanf( "%d", &n );

    PR( "Your line: %s\n", string );
    count = istrlen( string );
    PR( "Locations: ", string );
    for( i = 0; i < count; i++ ) PR( "%d", i % 10 );
    PR( "\n" );
    
    count = cntToken( string ); 
    PR( "The line has %d regular token(s).\n", count );

    count = cntTokenS( string );
    PR( "The line has %d quoted token(s).\n", count );

    token = getToken( n, string );
    PR( "Token(%d) Status(%d): %s\n", n, badToken, token ); 

    token = getTokenS( n, string );
    PR( "Quoted Token(%d) Status(%d): %s\n", n, badToken, token ); 

    n     = getTokenLoc( n, string, &loc1, &locN );
    PR( "Quoted Token(%d) Status(%d) Loc %d to %d\n", n, badToken, loc1, locN ); 
    ival  = getTokenI( n, string );
    PR( "TokenI(%d) Status(%d): %d\n", n, badToken, ival ); 

    rval  = getTokenR( n, string );
    PR( "TokenR(%d) Status(%d): %f\n", n, badToken, rval ); 
    
    token  = leftToken( n, string );
    PR( "LeftToken(%d) Status(%d): %s\n", n, badToken, token ); 

    token  = rightToken( n, string );
    PR( "RightToken(%d) Status(%d): %s\n", n, badToken, token ); 

    return( 0 );
}
#endif

int cntToken( string )

    char *string;
{
    int wasSpace, count;

    count = 0; wasSpace = 1; badToken = 0;

    if (!string) return( 0 );

    while( *string )
       {
        if (wasSpace)
           {
            if (!isspace( (int) *string ))
               {
                count++;
                wasSpace = 0;
               } 
           }
        else
           {
            if (isspace( (int) *string ))
               {
                wasSpace = 1;   
               }
           }

        string++;
       }

    return( count );
}

int cntTokenS( string )

    char *string;
{
    int  wasSpace, count;
    char qc;

    count = 0; wasSpace = 1; badToken = 0;

    if (!string) return( 0 );

    while( *string )
       {
        if (wasSpace)
           {
            if (!isspace( (int) *string ))
               {
                count++;
                wasSpace = 0;

                if (*string == '"' || *string == '\'')
                   {
                    qc = *string++;

                    while( *string && *string != qc) string++;
                   }
               }
           }
        else
           {
            if (isspace( (int) *string ))
               {
                wasSpace = 1;
               }
           }

        string++;
       }

    return( count );
}

int cntTokenC( string, tokChars )

    char *string, *tokChars;
{
    int  wasToken, count;
    char qc, *sPtr;

    badToken  = 0;
    count     = 0;
    wasToken  = 1;

    sPtr = string;

    if (!string || !tokChars)
       {
        badToken = 1;
        return( count );
       }

    while( *string )
       {
        if (wasToken)
           {
            count++;

            while( *string )
               {
                if (strchr( tokChars, (int) *string ))
                   {
                    wasToken = 1;
                    break;
                   }
                else
                   {
                    wasToken = 0;

                    if (*string == '"' || *string == '\'')
                       {
                        qc = *string++;
                        while( *string && *string != qc ) string++;
                       }
                   }

                string++;
               }
           }
        else
           {
            if (strchr( tokChars, (int) *string ))
               {
                wasToken = 1;
               }
           }

        if (!*string) break;

        string++;
       }

    return( count );
}

float getTokenR( n, string )

    int   n;
    char *string;
{
    float value;

    value    = 0.0;
    badToken = 0;

    if (!string)
       {
        badToken = 1;
        return( value );
       }

    if (!sscanf( getToken( n, string ), "%f", &value ))
       {
        badToken = 1;
       }

    return( value );
}
      
int getTokenI( n, string )

    char *string;
    int   n;
{
    int value;

    value    = 0;
    badToken = 0;

    if (!string)
       {
        badToken = 1;
        return( value );
       }

    if (!sscanf( getToken( n, string ), "%d", &value ))
       {
        badToken = 1;
       }

    return( value );
}

NMR_INT getTokenN( n, string )

    char *string;
    int   n;
{
    NMR_INT value;

    value    = 0;
    badToken = 0;

    if (!string || n < 1)
       {
        badToken = 1;
        return( value );
       }

#ifdef NMR64
    if (!sscanf( getToken( n, string ), "%ld", &value ))
       {
        badToken = 1;
       }
#else
    if (!sscanf( getToken( n, string ), "%d", &value ))
       {
        badToken = 1;
       }
#endif

    return( value );
}

char *leftToken( n, string )

    int  n;
    char *string;
{
    int oldN, i, count;

    badToken = 0;

    if (!lValue)
       {
        lBuffLen = NAMELEN + 1;
        lValue   = lBuff;
        lBuff[0] = '\0';
       }

    lValue[0] = '\0';

    if (n < 1)
       {
        return( lValue );
       }

    if (!string)
       {
        badToken = 1;
        return( lValue );
       }

    count = cntToken( string );
    oldN  = n;

    n = n > count ? count : n;
    n = n < 1 ? 1 : n;

    for( i = 1; i <= n; i++ )
       {
        (void) strcatAlloc( &lValue,
                            getToken( i, string ),
                            lBuff, NAMELEN+1,
                            &lBuffLen ); 

        if (i != n)
           {
            (void) strcatAlloc( &lValue, " ", lBuff, NAMELEN+1, &lBuffLen );
           }
       }

    if (oldN > count) badToken = 1;

    return( lValue );
}

char *rightToken( n, string )

    char *string;
    int  n;
{
    int i, oldN, count;

    badToken = 0;

    if (!rValue)
       {
        rBuffLen = NAMELEN + 1;
        rValue   = rBuff;
        rBuff[0] = '\0';
       }

    rValue[0] = '\0';

    if (!string)
       {
        badToken = 1;
        return( rValue );
       }

    count  = cntToken( string );
    oldN   = n;

    if (n > count)
       {
        return( rValue );
       }
 
    n = n > count ? count : n;
    n = n < 1 ? 1 : n;

    for( i = n; i <= count; i++ )
       {
        (void) strcatAlloc( &rValue,
                            getToken( i, string ),
                            rBuff, NAMELEN+1,
                            &rBuffLen );

        if (i != count)
           {
            (void) strcatAlloc( &rValue, " ", rBuff, NAMELEN+1, &rBuffLen );
           }
       }

    if (oldN > count) badToken = 1;

    return( rValue );
}

char *getToken( n, string )

    int   n;
    char *string;
{
    int    i, count, skip, nn;
    char   *sPtr;

    if (!tValue)
       {
        tBuffLen = NAMELEN + 1;
        tValue   = tBuff;
       }

    badToken  = 0;
    tValue[0] = '\0';

    if (!string || n < 1)
       {
        badToken = 1;
        return( tValue );
       }
    
    count = cntToken( string );

    if ((n < 1) || (n > count))
       {
        badToken = 1;
        return( tValue );
       }

    nn = istrlen( string ) + 1;

    if (nn > tBuffLen)
       {
        if (tValue != tBuff) 
           {
            (void) deAlloc( "token", tValue, sizeof(char)*tBuffLen );
           }

        tBuffLen = NAMELEN + 1;
        tValue   = tBuff;

        if ((sPtr = charAlloc( "token", nn )))
           {
            tBuffLen = nn;
            tValue   = sPtr;
           }
        else
           {
            FPR( stderr, "GetToken: Memory Allocation Failure.\n" );
           }

        tValue[0] = '\0';
       }

    sPtr = string;

    for( i = 0; i < n; i++ )
       {
        (void) sscanf( sPtr, "%s%n", tValue, &skip );
        sPtr += skip;
       } 

    return( tValue );
}

int findToken( str, key )

   char *str, *key;
{
    int n, count;

    count = cntToken( str );
    if (count < 1) return( 0 );

    for( n = 1; n <= count; n++ )
       {
        if (!strcmp( key, getToken( n, str ))) return( n );
       }

    return( 0 );
}

int getTokenMaxLength( str )

    char *str;
{
    int i, n, thisLength, maxLength;

    n         = cntToken( str );
    maxLength = 0;

    for( i = 1; i <= n; i++ )
       {
        thisLength = strlen( getToken( i, str ));
        if (thisLength > maxLength) maxLength = thisLength;
       }

    return( maxLength ); 
}

int getTokenLoc( n, strOrig, loc1, locN )

    char *strOrig;
    int  n, *loc1, *locN;
{
    int    loc, wasSpace, count, nn;
    char   *string, *sPtr, qc;

    badToken  = 0;
    count     = 0; 
    wasSpace  = 1;
    *loc1     = -1;
    *locN     = -1;

    if (!strOrig)
       {
        badToken = 1;
        return( 0 );
       }

    string = strOrig;

    while( *string )
       {
        if (wasSpace)
           {
            if (!isspace( (int)*string ))
               {
                count++;
                wasSpace = 0;

                if (*string == '"' || *string == '\'')
                   {
                    qc = *string++;

                    if (n == count)
                       {
                        while( *string && *string != qc)
                           {
                            loc = string - strOrig;

                            if (*loc1 < 0) *loc1 = loc;

                            *locN   = loc;
                            string++;
                           }
                        
                        return( 1 + *locN - *loc1 );
                       }
                    else
                       {
                        while( *string && *string != qc) string++;
                       }
                   }
                else
                   {
                    if (n == count)
                       {
                        while( *string && !isspace((int)*string) )
                           {
                            loc = string - strOrig;

                            if (*loc1 < 0) *loc1 = loc;

                            *locN = loc;
                            string++;
                           }

                        return( 1 + *locN - *loc1 );
                       }
                   }
               }
           }
        else
           {
            if (isspace( (int) *string ))
               {
                wasSpace = 1;
               }
           }

        string++;
       }

    badToken = 1;

    return( 0 );
}

char *getTokenS( n, string )

    char *string;
    int  n;
{
    int    wasSpace, count, nn;
    char   *sPtr, qc;

    if (!tValue)
       {
        tBuffLen = NAMELEN + 1;
        tValue   = tBuff;
       }

    badToken  = 0;
    count     = 0; 
    wasSpace  = 1;

    tValue[0] = '\0';

    if (!string || n < 1)
       {
        badToken = 1;
        return( tValue );
       }

    nn = istrlen( string ) + 1;

    if (nn > tBuffLen)
       {
        if (tValue != tBuff)
           {
            (void) deAlloc( "token", tValue, sizeof(char)*tBuffLen );
           }

        tBuffLen = NAMELEN + 1;
        tValue   = tBuff;

        if ((sPtr = charAlloc( "token", nn )))
           {
            tBuffLen = nn;
            tValue   = sPtr;
           }
        else
           {
            FPR( stderr, "GetToken: Memory Allocation Failure.\n" );
           }

        tValue[0] = '\0';
       }

    sPtr = tValue;

    while( *string )
       {
        if (wasSpace)
           {
            if (!isspace( (int)*string ))
               {
                count++;
                wasSpace = 0;

                if (*string == '"' || *string == '\'')
                   {
                    qc = *string++;

                    if (n == count)
                       {
                        while( *string && *string != qc)
                           {
                            *sPtr++ = *string++;
                           }
                        
                        *sPtr++ = '\0';

                        return( tValue );
                       }
                    else
                       {
                        while( *string && *string != qc) string++;
                       }
                   }
                else
                   {
                    if (n == count)
                       {
                        while( *string && !isspace((int)*string) )
                           {
                            *sPtr++ = *string++;
                           }

                        *sPtr++ = '\0';

                        return( tValue );
                       }
                   }
               }
           }
        else
           {
            if (isspace( (int) *string ))
               {
                wasSpace = 1;
               }
           }

        string++;
       }

    badToken = 1;

    return( tValue );
}

char *getTokenC( n, string, tokChars )

    char *string, *tokChars;
    int  n;
{
    int  nn, wasToken, count;
    char qc, *sPtr;

    if (!cValue)
       {
        cBuffLen = NAMELEN + 1;
        cValue   = cBuff;
       }

    badToken  = 0;
    count     = 0;
    wasToken  = 1;

    cValue[0] = '\0';

    if (!string || !tokChars || n < 1)
       {
        badToken = 1;
        return( cValue );
       }

    nn = istrlen( string ) + 1;

    if (nn > cBuffLen)
       {
        if (cValue != cBuff)
           {
            (void) deAlloc( "token", cValue, sizeof(char)*cBuffLen );
           }

        cBuffLen = NAMELEN + 1;
        cValue   = cBuff;

        if ((sPtr = charAlloc( "token", nn )))
           {
            cBuffLen = nn;
            cValue   = sPtr;
           }
        else
           {
            FPR( stderr, "GetToken: Memory Allocation Failure.\n" );
           }

        cValue[0] = '\0';
       }

    sPtr = cValue;

    while( *string )
       {
        if (wasToken)
           {
            count++;

            if (n == count)
               {
                while( *string && !strchr( tokChars, (int) *string ))
                    {
                     if (*string == '"' || *string == '\'')
                        {
                         qc = *string++;
                         while( *string && *string != qc ) *sPtr++ = *string++;
                        }
                     else
                        { 
                         *sPtr++ = *string;
                        }

                     string++;
                    }

                 *sPtr++ = '\0';

                 return( cValue );
                }

            while( *string )
               {
                if (strchr( tokChars, (int) *string ))
                   {
                    wasToken = 1;
                    break;
                   }
                else
                   {
                    wasToken = 0;

                    if (*string == '"' || *string == '\'')
                       {
                        qc = *string++;
                        while( *string && *string != qc ) string++;
                       }
                   }

                string++;
               }
           }
        else
           {
            if (strchr( tokChars, (int) *string ))
               {
                wasToken = 1;
               }
           }

        if (!*string) break;

        string++;
       }

    badToken = 1;

    return( cValue );
}

int strTrimLoc( s, loc1, locN )

   char *s;
   int  *loc1, *locN;
{
    int  i, n;
    char *sPtr;

    *loc1 = -1;
    *locN = -1;

    if (!s) return( 0 );

    n = istrlen( s );
    if (!n) return( 0 );

    *loc1 = 0;
    *locN = n - 1;

    for( sPtr = s; *sPtr; sPtr++ )
       {
        if (!isspace( *sPtr )) break;
        (*loc1)++;
       } 

    for( sPtr = s + n - 1; *sPtr; sPtr-- )
       {
        if (!isspace( *sPtr )) break;
        (*locN)--;
       }

    return( 1 + *locN - *loc1 );
}

int isInteger( s )

   char *s;
{
   int  i, loc1, locN, nC, nD;

   if (!s) return( 0 );

   loc1 = -1;
   locN = -1;

   nC = strTrimLoc( s, &loc1, &locN );
   if (nC < 1) return( 0 ); 

   nD = 0;

   for( i = loc1; i <= locN; i++ ) 
      {
       if (s[i] == '-' || s[i] == '+')
          {
           if (i > loc1) return( 0 );
          }
       else if (s[i] >= '0' && s[i] <= '9')
          {
           nD++; 
          }
       else
          {
           return( 0 ); 
          }
      }

   if (nD) return( 1 );

   return( 0 );
}

int isFloat( s )

   char *s;
{
   int  i, loc1, locN, locE, nC, nD, nE, nS, nP, nX;

   if (!s) return( 0 );

   loc1 = -1;
   locN = -1;
   locE = -1;

   nC = strTrimLoc( s, &loc1, &locN );
   if (nC < 1) return( 0 );

   nD = 0;
   nE = 0;
   nS = 0;
   nP = 0;
   nX = 0;

   for( i = loc1; i <= locN; i++ )
      {
       if (s[i] == 'e' || s[i] == 'E')
          {
           if (nE || !nD) return( 0 );
           locE = i;
           nE++;
          }
       else if (s[i] == '-' || s[i] == '+')
          {
           if (nS > 1) return( 0 );
           nS++;

           if (i == loc1) continue;
           if (!nE) return( 0 );
           if (i != locE + 1)  return( 0 );
          }
       else if (s[i] == '.')
          {
           if (nP || nE) return( 0 );
           nP++;
          }
       else if (s[i] >= '0' && s[i] <= '9')
          {
           if (nE)
              nX++;
           else
              nD++;
          }
       else
          {
           return( 0 );
          }
      }

   if (nD < 1) return( 0 );
   if (nE && !nX) return( 0 );

   return( 1 );
}

int freeToken()
{
    (void) strAllocFree( &lValue, lBuff, NAMELEN+1, &lBuffLen );
    (void) strAllocFree( &rValue, rBuff, NAMELEN+1, &rBuffLen );
    (void) strAllocFree( &tValue, tBuff, NAMELEN+1, &tBuffLen );

    lValue = (char *)NULL; lBuffLen = 0;
    rValue = (char *)NULL; rBuffLen = 0;
    tValue = (char *)NULL; tBuffLen = 0;

    return( 0 );
}

int istrlen( s )

   char *s;
{
   if (!s) return( 0 );
   return( (int)strlen( s ) );
}

int str2hash( str, i1, iN )

   char   *str;
   int    i1, iN;
{
    int      ihash, itmp, length;
    unsigned int h, g;
    char     *p;

    if (i1 > iN)
       {
        itmp = i1;
        i1   = iN;
        iN   = itmp;
       }

    length = 1 + iN - i1;
 
    if (!length) return( 0 );

    h = 0;

    for( p = str; *p; p++ )
       {
        h = (h << 4) + *p;
        g = h & 0xf0000000;

        if (g)
           {
            h = h ^ (g >> 24);
            h = h ^ g;
           }
       }

    h %= length;

    ihash = i1 + (h < 0 ? -h : h);

    return( ihash );
}

/***/
/* Bottom.
/***/
