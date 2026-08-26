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

#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include <stdio.h>
#include <errno.h>

#include "inquire.h"
#include "prec.h"

#ifdef CONVEX
#define S_ISSOCK(X) (_S_ISSOCK((X)))
#endif

#ifdef WIN95
#define S_ISSOCK(S) ((S) & 0)
#define S_ISDIR(S)  ((S) & _S_IFDIR)
#define S_ISFIFO(S) ((S) & _S_IFIFO)
#define S_ISREG(S)  ((S) & _S_IFREG)
#define S_ISCHR(S)  ((S) & _S_IFCHR)
#endif


int safe_stat(); //AJB 14/09/22
int safe_fstat(); //AJB 14/09/22

static char *lastChar();

/***/
/* fileExists: return 1 if file testName can be opened and read.
/*             return 0 otherwise.
/***/

int fileExists( testName )

   char *testName;

{
    FILE *testUnit;
    char  buffer[1];

    if (!(testUnit = fopen( testName, "r" )))
       return( 0 );

    if (sizeof(buffer) != fread( buffer, sizeof(buffer), 1, testUnit ))
       return( 0 );

    (void) fclose( testUnit );

    return( 1 );
} 

/***/
/* fileExists0: return 1 if file testName exists, even if it is zero length.
/*              return 0 otherwise.
/***/

int fileExists0( testName )

   char *testName;
{
    struct stat statBuff;
    int    status;

    errno = 0;

    status = safe_stat( testName, &statBuff );

    if (errno == ENOENT || status != 0)
       {
        return( 0 );
       }

    return( 1 );
}

/***/
/* getDirName: extract directory from file path. 
/***/

int getDirName( thisPath, dirName )

   char *thisPath;
   char *dirName;
{
    char *sPtr;
    int  status;

    status = 0;

    (void) strcpy( dirName, thisPath );

    if ((sPtr = lastChar( dirName, '/' )))
       {
        *sPtr  = '\0';
        status = 1;
       } 
    else
       {
        (void) strcpy( dirName, "." );
       }

    return( status );
}

/***/
/* getFileName: extract file name from full path.
/***/

int getFileName( thisPath, fileName )

   char *thisPath;
   char *fileName;
{
    (void) strcpy( fileName, fileTail( thisPath ));
    return( 0 );
}

/***/
/* dirExists: return 1 if path to file testName exists.
/*            return 0 otherwise.
/*            
/*            On return, dirName contains the path itself.
/***/

int dirExists( testName, dirName )

   char *testName, *dirName;
{
    struct stat statBuff;
    int    status;
    char   *sPtr;

    if (testName && *testName == '!') return( 1 );

    (void) strcpy( dirName, testName );

    errno = 0;

    if ((sPtr = lastChar( dirName, '/' )))
       {
        *sPtr = '\0';

        status = safe_stat( dirName, &statBuff );

        if (errno == ENOENT || status != 0 || !S_ISDIR(statBuff.st_mode))
           {
            return( 0 );
           }
       }
    else
       {
        (void) strcpy( dirName, "./" );
       }

    return( 1 );
}

/***/
/* isTTY: returns 1 if descriptor fd seems to be connected to a TTY.
/*        returns 0 otherwise.
/***/

int isTTY( fd )

   int fd;
{
    struct stat statBuff;

    if (safe_fstat( fd, &statBuff )) return( 0 );

    if (S_ISSOCK(statBuff.st_mode) ||
         S_ISFIFO(statBuff.st_mode) ||
          S_ISDIR(statBuff.st_mode)  ||
           S_ISREG(statBuff.st_mode)) return( 0 );

    if (S_ISCHR(statBuff.st_mode)) return( 1 );

    return( 0 );
}

char *fileTail( fName )

   char *fName;
{
    static char tailName[NAMELEN+1];
    char   *sPtr;

    tailName[0] = '\0';

    if (strlen( fName ) > NAMELEN) return( (char *)NULL );

    if ((sPtr = lastChar( fName, '/' )))
       (void) strcpy( tailName, sPtr+1 );
    else
       (void) strcpy( tailName, fName );

    return( tailName );
}

char *fileRoot( fName )

   char *fName;
{
    static char rootName[NAMELEN+1];
    char   *sPtr;

    rootName[0] = '\0';

    if (strlen( fName ) > NAMELEN) return( (char *)NULL );

    (void) strcpy( rootName, fName );

    if ((sPtr = lastChar( rootName, '.' ))) *sPtr = '\0';

    return( rootName );
}

char *fileExt( fName )

   char *fName;
{
    static char extName[NAMELEN+1];
    char   *sPtr;

    extName[0] = '\0';

    if (strlen( fName ) > NAMELEN) return( (char *)NULL );

    if ((sPtr = lastChar( fName, '.' )))
       (void) strcpy( extName, sPtr+1 );
    else
       *extName = '\0';

    return( extName );
}

char *fileDir( fName )

   char *fName;
{
    static char dirName[NAMELEN+1];
    char   *sPtr;

    dirName[0] = '\0';

    if (strlen( fName ) > NAMELEN) return( (char *)NULL );

    (void) strcpy( dirName, fName );

    if ((sPtr = lastChar( dirName, '/' )))
       *sPtr = '\0';
    else
       (void) strcpy( dirName, "." );

    return( dirName );
}

/***/
/* lastChar: return ptr to last occurrence of c in string s.
/***/

static char *lastChar( s, c )

   char *s, c;
{
    int i;

    if (!s) return( (char *)NULL );

    i = (int)strlen( s );

    while( i-- )
       {
        if (s[i] == c) return( s+i );
       }

    return( (char *)NULL );
}

