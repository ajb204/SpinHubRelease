/***/
/*  NMRPipe Copyright 1992-2018 Frank Delaglio
/***/

/***/
/* DataIO: routines for binary I/O via C or system calls.
/*
/* error = DataRead( fileUnit, buffer, byteCount );
/*         DataWrite( fileUnit, buffer, byteCount );
/*         DataReadB( fileUnit, buffer, byteCount, maxTries, timeOut );
/*         DataWriteB( fileUnit, buffer, byteCount, maxTries, timeOut );
/*         DataRecvB( fileUnit, buffer, byteCount, maxTries, timeOut );
/*         DataSendB( fileUnit, buffer, byteCount, maxTries, timeOut );
/*         DataPos( fileUnit, bytePosition );
/*         DataPos2( fileUnit, bytePosition, buff, buffSize );
/*
/* error = dataOpen( fileName, fileUnit, openMode );
/*         dataOpen2( fileName, fileUnit, openStr );
/*         dataClose( fileUnit );
/*         dataCreate( fileName );
/*         dirCreate( dirName );
/**/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef WIN95
#include <sys/stat.h>
#include <errno.h>
#include <winsock.h>
#include <direct.h>
#include <fcntl.h>

#define popen(S1,S2) _popen(S1,S2)
typedef unsigned int mode_t;
#else
#include <sys/errno.h>
#include <dirent.h>
#endif

#ifdef MAC_OSX
#include <sys/stat.h>
#endif

#include <string.h>

#include "dataio2.h"
#include "dataio.h"
#include "inquire.h"
#include "testsize.h"

#ifdef NMRTIME
#include "nmrserver.h"
#include "nmrpipe.h"
#endif

#include "nmrtime.h"
#include "memory.h"

#ifdef FB_SYS_IO

#ifdef WIN95
#include <time.h>
#undef O_BINARY
#define O_BINARY _O_BINARY | _O_RANDOM
#else
#include <sys/time.h>
#include <sys/file.h>
#undef  O_BINARY
#define O_BINARY 0
#endif

#include <fcntl.h>
#include <errno.h>

#ifndef WIN95
#include <unistd.h>
#endif

#endif /* FB_SYS_IO */

#include "memory.h"
#include "lintfix.h"
#include "prec.h"

#ifdef WINXP
#include <fcntl.h>
#endif

#if defined (WIN95) || defined (WINNT) || defined (MAC_OSX)
static int selectFlag = 0;
#else
static int selectFlag = 0;
#endif

static NMR_INT rdCount     = 0;
static NMR_INT rdRetries   = 0;
static NMR_INT rdWaitCount = 0;

static NMR_INT wrCount     = 0;
static NMR_INT wrRetries   = 0;
static NMR_INT wrWaitCount = 0;

static int defTimeOut = DATAIO_TIMEOUT;

static int autoSwapFlag   = 1;
static int byteSwapFlag   = 0;
static int blockingWrFlag = 0;
static int fsCheckFlag    = 0;
static int stdBufInFlag   = 0;
static int stdBufOutFlag  = 0;
static int stdBufInSize   = 0;
static int stdBufOutSize  = 0;
static int initDone       = 0;
static int ioDebugFlag    = 0;

static char *stdBufIn     = (char *)NULL;
static char *stdBufOut    = (char *)NULL;

static int  checkIOError  = 1;

#define FPR (void)fprintf


int safe_stat(); //AJB 14/09/22
ssize_t recv();//AJB 14/09/22
ssize_t send();//AJB 14/09/22
int dataWriteCB();//AJB 14/09/22

/***/
/* initDataIO: initialize I/O settings via environment variables.
/***/

int initDataIO()
{
    int  status;
    char *sPtr;

    if (initDone) return( 0 );

    status = 0;

    if ((sPtr = getenv( "NMR_IO_SELECT" )))
       {
        selectFlag = atoi( sPtr );
        status += 1;
       }

    if ((sPtr = getenv( "NMR_IO_TIMEOUT" )))
       {
        defTimeOut = atoi( sPtr );
        status += 2;
       }

    if ((sPtr = getenv( "NMR_AUTOSWAP" )))
       {
        if (!strcasecmp( sPtr, "VERBOSE" )) autoSwapFlag = 2;
        if (!strcasecmp( sPtr, "YES" ))     autoSwapFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    autoSwapFlag = 1;
        if (!strcasecmp( sPtr, "1" ))       autoSwapFlag = 1;

        if (!strcasecmp( sPtr, "NO" ))      autoSwapFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   autoSwapFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       autoSwapFlag = 0;

        status += 4;
       }

    if ((sPtr = getenv( "NMR_FSCHECK" )))
       {
        if (!strcasecmp( sPtr, "YES" ))     fsCheckFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    fsCheckFlag = 1;
        if (!strcasecmp( sPtr, "1" ))       fsCheckFlag = 1;

        if (!strcasecmp( sPtr, "NO" ))      fsCheckFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   fsCheckFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       fsCheckFlag = 0;

        status += 8;
       }

    if ((sPtr = getenv( "NMR_BLOCKING_WR" )))
       {
        if (!strcasecmp( sPtr, "YES" ))     blockingWrFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    blockingWrFlag = 1; 
        if (!strcasecmp( sPtr, "1" ))       blockingWrFlag = 1; 

        if (!strcasecmp( sPtr, "NO" ))      blockingWrFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   blockingWrFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       blockingWrFlag = 0;

        status += 16;
       }

    if ((sPtr = getenv( "NMR_DEBUG_IO" )))
       {
        if (!strcasecmp( sPtr, "YES" ))     ioDebugFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    ioDebugFlag = 1; 
        if (!strcasecmp( sPtr, "1" ))       ioDebugFlag = 1;

        if (!strcasecmp( sPtr, "NO" ))      ioDebugFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   ioDebugFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       ioDebugFlag = 0;

        status += 16;
       }

/*
   if ((sPtr = getenv( "NMR_STDBUF_IN" )))
       {
        if (!strcasecmp( sPtr, "YES" ))     stdBufInFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    stdBufInFlag = 1;
        if (!strcasecmp( sPtr, "1" ))       stdBufInFlag = 1;

        if (!strcasecmp( sPtr, "NO" ))      stdBufInFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   stdBufInFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       stdBufInFlag = 0;

        status += 32;
       }

   if ((sPtr = getenv( "NMR_STDBUF_OUT" )))
       {
        if (!strcasecmp( sPtr, "YES" ))     stdBufOutFlag = 1;
        if (!strcasecmp( sPtr, "TRUE" ))    stdBufOutFlag = 1;
        if (!strcasecmp( sPtr, "1" ))       stdBufOutFlag = 1;

        if (!strcasecmp( sPtr, "NO" ))      stdBufOutFlag = 0;
        if (!strcasecmp( sPtr, "FALSE" ))   stdBufOutFlag = 0;
        if (!strcasecmp( sPtr, "0" ))       stdBufOutFlag = 0;

        status += 64;
       }
*/

    initDone = 1;

    return( status ); 
}

/***/
/* Report I/O debug info.
/***/

int showDebugIO( argc, argv )

   int  argc;
   char **argv;
{
    int i, pid;

    if (!ioDebugFlag) return( 0 );

    pid = getpid();

    FPR( stderr, "\n" );
    FPR( stderr, "PID %d Command:      ", pid );

    for( i = 0; i < argc; i++ ) FPR( stderr, " %s", argv[i] );

    FPR( stderr, "\n" );
    FPR( stderr, "PID %d Reads:         %ld\n", pid, (long)rdCount );
    FPR( stderr, "PID %d Read Retries:  %ld\n", pid, (long)rdRetries );
    FPR( stderr, "PID %d Writes:        %ld\n", pid, (long)wrCount );
    FPR( stderr, "PID %d Write Retries: %ld\n", pid, (long)wrRetries );
    FPR( stderr, "\n\n" );

    return( 0 );
}

int getIODebugFlag()
{
    return( ioDebugFlag );
}

/***/
/* stdBuf utilities.
/***/

int getStdBufIn()
{
   return( stdBufInFlag );
}

int getStdBufOut()
{
   return( stdBufOutFlag );
}

int setStdBufIn( n )

   int n;
{
    if (n < 0) return( 0 );

    if (n > 0)
       {
        if (!(stdBufIn = charAlloc( "stdBufIn", sizeof(float)*n )))
           {
            FPR( stderr, "Error Allocating stdin Buffer.\n" );
            return( 1 );
           }

        stdBufInSize = n;
       }

    FPR( stderr, "stdin Buffer: %d\n", n );

    (void) setbuffer( stdin, stdBufIn, sizeof(float)*n );

    return( 0 );
}

int setStdBufOut( n )

   int n;
{
    if (n < 0) return( 0 );

    if (n > 0)
       {
        if (!(stdBufOut = charAlloc( "stdBufOut", sizeof(float)*n )))
           {
            FPR( stderr, "Error Allocating stdout Buffer.\n" );
            return( 1 );
           }

        stdBufOutSize = n;
       }

    FPR( stderr, "stdout Buffer: %d\n", n );

    (void) setbuffer( stdout, stdBufOut, sizeof(float)*n );

    return( 0 );
}

int stdBufFree()
{
    (void) deAlloc( "stdBufIn",  stdBufIn,  sizeof(float)*stdBufInSize );
    (void) deAlloc( "stdBufOut", stdBufOut, sizeof(float)*stdBufOutSize );

    stdBufIn  = (char *)NULL;
    stdBufOut = (char *)NULL;

    stdBufInSize  = 0;
    stdBufOutSize = 0;

    return( 0 );
}

/***/
/* getDefTimeout: returns current default IO timeout.
/***/

int getDefTimeout()
{
   (void) initDataIO();
   return( defTimeOut );
}

int setDefTimeout( i )

   int i;
{
    (void) initDataIO();

    defTimeOut = i;

    return( i );
}

/***/
/* getFSCheck: returns current default setting for checking available file space.
/***/

int getFSCheck()
{
   (void) initDataIO();
   return( fsCheckFlag );
}

int setFSCheck( i )

   int i;
{
   (void) initDataIO();

   fsCheckFlag = i;

   return( i );
}

/***/
/* getIOSelect: returns current default setting for use of select() for blocking I/O.
/***/

int getIOSelect()
{
   (void) initDataIO();
   return( selectFlag );
}

int setIOSelect( i )

   int i;
{
   (void) initDataIO();

   selectFlag = i;

   return( i );
}

/***/
/* dataOpen: open a binary file named "fileName".
/***/

int dataOpen( fileName, fileUnit, openMode )

   char *fileName;
   int   openMode;

   FILE_UNIT( *fileUnit );

{
    struct stat buff;
    char   *modeName;
    int    error;
    mode_t mode;

/***/
/* Initialize error status:
/***/

    error = 0;
    mode  = 0666;

    (void) initDataIO();

/***/
/* Test for I/O to one-way pipe:
/***/

   if (fileName && *fileName == '!')
      {
       error = dataPOpen( fileName+1, fileUnit, openMode );
       return( error );
      }

/***/
/* Test for existance of file and create it if needed:
/***/

    if (openMode == FB_UNKNOWN)
       {
        if (safe_stat( fileName, &buff ) && ENOENT == errno)
           {
            if ((error = dataCreate( fileName )))
               {
                return( error );
               }
           }

        openMode = FB_READWRITE;
       }

/***/
/* Open file using system call:
/***/

#ifdef FB_SYS_IO

    errno = 0;

    switch( openMode )
       {
        case FB_READ:
           *fileUnit = open( fileName, O_RDONLY | O_BINARY );
           modeName  = "Read";
           break;

        case FB_WRITE:
           *fileUnit = open( fileName, O_WRONLY | O_BINARY );
           modeName  = "Write";
           break;

        case FB_READWRITE:
           *fileUnit = open( fileName, O_RDWR | O_BINARY  );
           modeName  = "Read/Write";
           break;

        case FB_CREATE:

#ifdef WIN95
           if ((*fileUnit = creat( fileName, mode )) >= 0)
              {
               (void) close( *fileUnit );
              }

           *fileUnit = open( fileName, O_WRONLY | O_BINARY );
#else
           if ((*fileUnit = creat( fileName, mode )) < 0)
              {
               *fileUnit = open( fileName, O_WRONLY | O_BINARY );
              }
#endif

           modeName  = "Create";
           break;

        case FB_CREATE_KEEP:
           *fileUnit = open( fileName, O_RDWR | O_BINARY | O_CREAT, (mode_t)0600 );
           modeName  = "CreateKeep";
           break;

        default:
           *fileUnit = -1;
           modeName  = "???";
           break;
       }
     
    if (*fileUnit < 0)
       {
        (void) fprintf( stderr, "DATAIO SYSOPEN Error %d:\n", errno );
        (void) fprintf( stderr, "File: %s Access: %s\n", fileName, modeName );

        if (errno) (void) perror( "SYSOPEN Message" );

        error = 1;
        *fileUnit = UNIT_CAST( 0 );  
       }

/***/
/* Open file via C routines:
/***/

#else

    switch( openMode )
       {
        case FB_READ:
           *fileUnit = fopen( inName, "r" );
           break;
        case FB_WRITE:
           *fileUnit = fopen( inName, "w" );
           break;
        case FB_READWRITE:
           *fileUnit = fopen( inName, "r+" );
           break;
        case FB_CREATE:
           *fileUnit = fopen( inName, "w+" );
           break;
        default:
           *fileUnit = UNIT_CAST( 0 );
           break;
       }

    if (!*fileUnit)
       {
        (void) fprintf( stderr, "DATAIO C-OPEN Error:\n" );
        (void) fprintf( stderr, "File %s Mode %d`n", openMode );
        error = 1;
       }

#endif

/***/
/* Return open status either way:
/***/

    return( error );
}

/***/
/* dataPOpen: open a one-way pipe according to mode string.
/***/

int dataPOpen( cmndName, fileUnit, openMode )

   char *cmndName;
   int   openMode;

   FILE_UNIT( *fileUnit );

{
    FILE  *pUnit;
    int   error;

/***/
/* Initialize error status:
/***/

    error     = 0;
    *fileUnit = -1;

    (void) initDataIO();

    switch( openMode )
       {
        case FB_READ:
           if ((pUnit = popen( cmndName, "r" )))
              {
               *fileUnit = fileno( pUnit );
              }

           break;

        case FB_WRITE:
        case FB_CREATE:

           if ((pUnit = popen( cmndName, "w" )))
              {
               *fileUnit = fileno( pUnit );
              }

           break;

        default:
           break;
       }

    if (*fileUnit < 0)
       {
        (void) fprintf( stderr, "DATAIO POPEN Error %d:\n", errno );
        (void) fprintf( stderr, "Cmnd: %s Access: %d\n", cmndName, openMode );

        if (errno) (void) perror( "SYSOPEN Message" );

        error = 1;
       }

    return( error );
}

/***/
/* dataOpen2: open a binary file named "fileName" according to mode string.
/***/

int dataOpen2( fileName, fileUnit, openStr )

   char *fileName, *openStr;
   int  *fileUnit;

{
    int    pipeMode, flags;
    int    error;
    mode_t mode;

    *fileUnit = -1;

    pipeMode  = FB_READ; 
    mode      = 0666;
    flags     = 0;
    error     = 0;

    (void) initDataIO();

    if (strchr( openStr, 'r' ) || strchr( openStr, 'R' )) flags |= O_RDONLY;
    if (strchr( openStr, 'w' ) || strchr( openStr, 'W' )) flags |= O_WRONLY;

    if (flags == (O_RDONLY | O_WRONLY)) flags = O_RDWR;

#ifndef WIN95
    if (strchr( openStr, 'd' ) || strchr( openStr, 'D' )) flags |= O_NDELAY;
    if (strchr( openStr, 'y' ) || strchr( openStr, 'Y' )) flags |= O_NOCTTY;
#ifndef WINXP
#ifdef MAC_OSX
    if (strchr( openStr, 's' ) || strchr( openStr, 'S' )) flags |= O_FSYNC;
#else
    if (strchr( openStr, 's' ) || strchr( openStr, 'S' )) flags |= O_SYNC;
#endif
#endif
#endif
    if (strchr( openStr, 'a' ) || strchr( openStr, 'A' )) flags |= O_APPEND;
    if (strchr( openStr, 'c' ) || strchr( openStr, 'C' )) flags |= O_CREAT;
    if (strchr( openStr, 't' ) || strchr( openStr, 'T' )) flags |= O_TRUNC;
    if (strchr( openStr, 'e' ) || strchr( openStr, 'E' )) flags |= O_EXCL;

    if (!strcasecmp( openStr, "w" ) || !strcasecmp( openStr, "W" ))
       {
        flags = (O_RDWR | O_CREAT | O_TRUNC);
       }

    mode |= O_BINARY;

    if (flags | O_WRONLY) pipeMode = FB_WRITE;
    if (flags | O_APPEND) pipeMode = FB_WRITE;
    if (flags | O_CREAT)  pipeMode = FB_WRITE;

    if (fileName && *fileName == '!')
       {
        error = dataPOpen( fileName+1, fileUnit, pipeMode );
       }
    else
       {
        *fileUnit = open( fileName, flags, mode );
        if (*fileUnit < 0) (void) perror( "DataOpen2 File Error" );
       }

    return( error );
}

/***/
/* DataRead: read byteCount bytes from file inUnit.
/***/

int DataRead( inUnit, array, byteCount )

   FILE_UNIT( inUnit );

   REAL4   *array;
   NMR_INT byteCount;

{
    NMR_INT  count;
    int      error;

/***/
/* Initialize error status:
/***/

    error = 0;

    (void) initDataIO();

/***/
/* Read the data:
/***/

#ifdef FB_SYS_IO
    count = (NMR_INT)read( inUnit, (IOPTR *)array, (IOSIZE)byteCount );
#else
    count = (NMR_INT)fread( (IOPTR *)array, (IOSIZE)1, (IOSIZE)byteCount, inUnit );
#endif

   if (ioDebugFlag) rdCount++;

/***/
/* Interpret the status:
/***/

    if (count != byteCount)
       {
        (void) fprintf( stderr, "DATAIO File Read Error %d.\n", errno );
        (void) fprintf( stderr, "Request: %ld Actual: %ld.\n", (long)byteCount, (long)count );
        (void) perror( "SYSREAD Message" );
        error = 1;
       }

   if (byteSwapFlag) (void) byteSwapV( array, byteCount/sizeof(REAL4) );
 
   return( error );
}

/***/
/* Get and set byteswap modes.
/***/

int getAutoSwapFlag()
{
   return( autoSwapFlag );
}

int getByteSwapFlag()
{
   return( byteSwapFlag );
}

int setAutoSwapFlag( aFlag )

   int aFlag;
{
   autoSwapFlag = aFlag;
   return( aFlag );
}

int setByteSwapFlag( sFlag )

   int sFlag;
{
    byteSwapFlag = sFlag;
    return( sFlag );
}

/***/
/* DataReadB: Read byteCount bytes from file inUnit in pseudo-blocking mode.
/***/

int DataReadB( inUnit, rArray, bytesRequested, maxTries, timeOut )

   FILE_UNIT( inUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   REAL4   *rArray;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;
    char    *array;

    array     = (char *)rArray;
    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &readfds );
            FD_SET( inUnit, &readfds );

            count = (NMR_INT)select( DTABLESIZE, &readfds, FDNULL, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)read( inUnit, (IOPTR *)array, (IOSIZE)bytesLeft );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        count = (NMR_INT)fread( (IOPTR *)array, (IOSIZE)1, (IOSIZE)bytesLeft, inUnit );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( readRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) rdRetries++;
            rdCount++;
           }
       }

    if (byteSwapFlag) (void) byteSwapV( rArray, bytesRequested/sizeof(REAL4) );

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataRecvB: read byteCount bytes from socket inUnit in pseudo-blocking mode.
/***/

int DataRecvB( inUnit, rArray, bytesRequested, maxTries, timeOut )

   FILE_UNIT( inUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   REAL4   *rArray;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;
    char    *array;

    array     = (char *)rArray;
    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &readfds );
            FD_SET( inUnit, &readfds );

            count = (NMR_INT)select( DTABLESIZE, &readfds, FDNULL, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)recv( inUnit, (IOPTR *)array, (IOSIZE)bytesLeft, 0 );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        return( -666 );
#endif

        bytesLeft -= count;

        if (bytesLeft)
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( sReadRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) rdRetries++;
            rdCount++;
           }
       }

    if (byteSwapFlag) (void) byteSwapV( rArray, bytesRequested/sizeof(REAL4) );

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataRecvCB: read byteCount chars from socket inUnit in pseudo-blocking mode.
/***/

int DataRecvCB( inUnit, array, bytesRequested, maxTries, timeOut )

   FILE_UNIT( inUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   char    *array;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;

    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &readfds );
            FD_SET( inUnit, &readfds );

            count = (NMR_INT)select( DTABLESIZE, &readfds, FDNULL, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)recv( inUnit, (IOPTR *)array, (IOSIZE)bytesLeft, 0 );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }

#else
        return( -666 );
#endif

        bytesLeft -= count;

        if (bytesLeft)
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( sReadRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) rdRetries++;
            rdCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* dataReadC: read byteCount bytes from file inUnit into a character array.
/***/

int dataReadC( inUnit, array, byteCount )

   FILE_UNIT( inUnit );

   char     *array;
   NMR_INT  byteCount;
{
    NMR_INT count;
    int     error;

/***/
/* Initialize error status:
/***/

    error = 0;
    errno = 0;

    (void) initDataIO();

/***/
/* Read the data:
/***/

#ifdef FB_SYS_IO
    count = (NMR_INT)read( inUnit, (IOPTR *)array, (IOSIZE)byteCount );
#else
    count = (NMR_INT)fread( array, 1, byteCount, inUnit );
#endif

    if (ioDebugFlag) rdCount++;

/***/
/* Interpret the status:
/***/

    if (count != byteCount)
       {
        (void) fprintf( stderr, "DATAIO File Read Error %d.\n", errno );
        (void) fprintf( stderr, "Request: %ld Actual: %ld.\n", (long)byteCount, (long)count );

        if (errno) (void) perror( "SYSTEM Message" );

        error = 1;
       }

   return( error );
}

/***/
/* dataWriteC: write byteCount bytes from file outUnit into a character array.
/***/

int dataWriteC( outUnit, array, byteCount )

   FILE_UNIT( outUnit );

   char    *array;
   NMR_INT byteCount;
{
    NMR_INT count;
    int     error;

/***/
/* Initialize error status:
/***/

    error = 0;
    errno = 0;

    (void) initDataIO();

    if (blockingWrFlag)
       {
        error = dataWriteCB( outUnit, array, byteCount, -1, 0 );
        return( error );
       }

/***/
/* Write the data:
/***/

#ifdef FB_SYS_IO
    count = (NMR_INT)write( outUnit, (IOPTR *)array, (IOSIZE)byteCount );
#else
    count = (NMR_INT)fwrite( array, 1, byteCount, outUnit );
#endif

    if (ioDebugFlag) wrCount++;

/***/
/* Interpret the status:
/***/

    if (count != byteCount)
       {
        (void) fprintf( stderr, "DATAIO File Write Error %d.\n", errno );
        (void) fprintf( stderr, "Request: %ld Actual: %ld.\n", (long)byteCount, (long)count );

        if (errno) (void) perror( "SYSTEM Message" );

        error = 1;
       }

   return( error );
}

/***/
/* dataReadCB: Read byteCount bytes from file inUnit in pseudo-blocking mode.
/***/

int dataReadCB( inUnit, array, bytesRequested, maxTries, timeOut )

   FILE_UNIT( inUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   char    *array;

{
    NMR_INT count, bytesLeft;
    NMR_INT tries;

    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &readfds );
            FD_SET( inUnit, &readfds );

            count = (NMR_INT)select( DTABLESIZE, &readfds, FDNULL, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)read( inUnit, (IOPTR *)array, (IOSIZE)bytesLeft );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        count = (NMR_INT)fread( array, 1, bytesLeft, inUnit );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( readRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) rdRetries++;
            rdCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataWrite: write byteCount bytes to file outUnit.
/***/

int DataWrite( outUnit, array, byteCount )

   FILE_UNIT( outUnit );

   REAL4   *array;
   NMR_INT byteCount;
{
    NMR_INT count, bytesLeft;
    int     tries, error;

/***/
/* Initialize error status:
/***/

    error = 0;
    errno = 0;

    (void) initDataIO();

    if (blockingWrFlag)
       {
        error = dataWriteB( outUnit, array, byteCount, -1, 0 );
        return( error );
       }

/***/
/* Write the data:
/***/

#ifdef FB_SYS_IO
    count = (NMR_INT)write( outUnit, (IOPTR *)array, (IOSIZE)byteCount );
#else
    count = (NMR_INT)fwrite( array, 1, byteCount, outUnit );
#endif

    if (ioDebugFlag) wrCount++;

/***/
/* Interpret the status:
/***/

    if (count != byteCount)
       {
        (void) fprintf( stderr, "DATAIO File Write Error %d.\n", errno );
        (void) fprintf( stderr, "Request: %ld Actual: %ld.\n", (long)byteCount, (long)count );

        if (errno) (void) perror( "SYSWRITE Message" );

        error = 1;
       }

   return( error );
}

/***/
/* DataWriteB: write byteCount bytes to file outUnit in pseudo-blocking mode.
/***/

int DataWriteB( outUnit, rArray, bytesRequested, maxTries, timeOut )

   FILE_UNIT( outUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   REAL4   *rArray;

{
    NMR_INT count, bytesLeft;
    NMR_INT tries;
    char   *array;

    array     = (char *)rArray;
    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &writefds );
            FD_SET( outUnit, &writefds );

            count = (NMR_INT)select( DTABLESIZE, FDNULL, &writefds, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)write( outUnit, (IOPTR *)array, (IOSIZE)bytesLeft );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        count = (NMR_INT)fwrite( array, 1, bytesLeft, outUnit );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( writeRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) wrRetries++;
            wrCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* dataWriteCB: write byteCount chars to file outUnit in pseudo-blocking mode.
/***/

int dataWriteCB( outUnit, array, bytesRequested, maxTries, timeOut )

   FILE_UNIT( outUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   char    *array;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;

    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &writefds );
            FD_SET( outUnit, &writefds );

            count = (NMR_INT)select( DTABLESIZE, FDNULL, &writefds, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)write( outUnit, (IOPTR *)array, (IOSIZE)bytesLeft );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }

#else
        count = (NMR_INT)fwrite( array, 1, bytesLeft, outUnit );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( writeRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) wrRetries++;
            wrCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataSendB: write byteCount bytes to socket outUnit in pseudo-blocking mode.
/***/

int DataSendB( outUnit, rArray, bytesRequested, maxTries, timeOut )

   FILE_UNIT( outUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   REAL4   *rArray;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;
    char    *array;

    array     = (char *)rArray;
    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &writefds );
            FD_SET( outUnit, &writefds );

            count = (NMR_INT)select( DTABLESIZE, FDNULL, &writefds, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)send( outUnit, (IOPTR *)array, (IOSIZE)bytesLeft, 0 );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        return( -666 );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( sWriteRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) wrRetries++;
            wrCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataSendCB: write byteCount chars to socket outUnit in pseudo-blocking mode.
/***/

int DataSendCB( outUnit, array, bytesRequested, maxTries, timeOut )

   FILE_UNIT( outUnit );

   NMR_INT bytesRequested;
   int     maxTries, timeOut;
   char    *array;
{
    NMR_INT count, bytesLeft;
    NMR_INT tries;

    bytesLeft = bytesRequested;
    tries     = 0;

    (void) initDataIO();

#ifdef DATAIO_USE_MAX_TRIES
    while( bytesLeft && (maxTries < 0 || tries++ < maxTries) )
#else
    while( bytesLeft )
#endif
       {
#ifdef FB_SYS_IO
        if (selectFlag)
           {
            FD_ZERO2( &writefds );
            FD_SET( outUnit, &writefds );

            count = (NMR_INT)select( DTABLESIZE, FDNULL, &writefds, FDNULL, TONULL );
           }

        errno = 0;

        count = (NMR_INT)send( outUnit, (IOPTR *)array, (IOSIZE)bytesLeft, 0 );

        if (checkIOError && (count < 0 || (count == 0 && errno != 0)))
           {
            if (errno == EAGAIN || errno == EWOULDBLOCK)
               count = 0;
            else
               return( -1 );
           }
#else
        return( -666 );
#endif

        bytesLeft -= count;

#ifdef DATAIO_USE_MAX_TRIES
        if (bytesLeft && (maxTries < 0 || tries < maxTries))
#else
        if (bytesLeft)
#endif
           {
            array += count;
            USLEEP( timeOut );
            TIMER_INCR( sWriteRetries );
           }

        if (ioDebugFlag)
           {
            if (bytesLeft > 0) wrRetries++;
            wrCount++;
           }
       }

    if (bytesLeft) return( -1 );

    return( 0 );
}

/***/
/* DataTruncU: truncate or expand file to the given size in bytes.
/***/

int DataTruncU( thisUnit, targetBytes )

   NMR_INT targetBytes;
   UNIT    thisUnit;
{
   int error;

   error = ftruncate( thisUnit, (off_t)targetBytes );

   return( error );
}

/***/
/* DataTrunc: truncate or expand file to the given size in bytes.
/***/

int DataTrunc( thisName, targetBytes, origBytes, vFlag )

   NMR_INT targetBytes, *origBytes;
   char    *thisName;
   int     vFlag;
{
    struct  FileSize trueSizeInfo;
    NMR_INT trueBytes;

    if (origBytes) *origBytes = 0;

    if (!fileExists( thisName ))
       {
        if (vFlag) FPR( stderr, "Truncate Error Finding File %s\n", thisName );
        return( 1 );
       }

    (void) getFileBytes( thisName, &trueSizeInfo );

    trueBytes  = trueSizeInfo.iTotalBytes;
    if (origBytes) *origBytes = trueBytes;

    if (trueBytes == targetBytes)
       {
#ifdef NMR64
        FPR( stderr, 
             "%s File Unchanged. Bytes Before: %ld After: %ld\n",
             thisName, trueBytes, targetBytes ); 
#else
        FPR( stderr,
             "%s File Unchanged. Bytes Before: %d After: %d\n",
             thisName, trueBytes, targetBytes );
#endif

        return( 0 );
       }

    if (vFlag)
       {
#ifdef NMR64
        FPR( stderr, 
             "%s File %s Bytes Before: %ld After: %ld\n",
             trueBytes > targetBytes ? "Truncating." : "Expanding.",
             thisName, trueBytes, targetBytes ); 
#else
        FPR( stderr, 
             "%s File %s Bytes Before: %d After: %d\n",
             trueBytes > targetBytes ? "Truncating." : "Expanding.",
             thisName, trueBytes, targetBytes );
#endif
       }

    if (truncate( thisName, (off_t)targetBytes )) return( 1 );

    return( 0 );
}
 
/***/
/* DataPos: position fileUnit to byte "position".
/***/

int DataPos( fileUnit, position )

   FILE_UNIT( fileUnit );
   NMR_INT position;
{
    int error;
    
#ifdef FB_SYS_IO
    off_t istatus;
#endif

    error = 0;
    errno = 0;

#ifdef FB_SYS_IO
    istatus = lseek( fileUnit, (off_t) position, SEEK_SET );
    if (-1 == istatus) error = 1;
#else
    error = fseek( fileUnit, (long) position, SEEK_SET );
#endif

    if (error)
       {
        (void) fprintf( stderr, "DATAIO File Seek error %d.\n", errno );

        if (errno) (void) perror( "SYSTEM Message" );
       }
      
    return( error ); 
}

/***/
/* DataPos2: use read to move position of fileUnit forward by byteCount.
/***/

int DataPos2( fileUnit, byteCount, buffer, buffSize )

   FILE_UNIT( fileUnit );
   NMR_INT byteCount, buffSize;
   VOID    *buffer;
{
    NMR_INT bytesIn;

    while( byteCount > 0 )
       {
        bytesIn = byteCount > buffSize ? buffSize : byteCount;

        if (dataReadB( fileUnit, (float *)buffer, bytesIn, -1, 0 ))
           {
            (void) fprintf( stderr, "DATAPOS2 Error positioning stream.\n" );
            return( 1 );
           }

        byteCount -= bytesIn;
       }

    return( 0 );
}

/***/
/* ByteSwapV: verbose version, perform 4-byte swap of elements in vec.
/***/

int byteSwapV( vec, length )

   float *vec;
   NMR_INT   length;
{
    union fsw { float f; char s[4]; } in, out;

    if (autoSwapFlag > 1)
       {
#ifdef NMR64
        FPR( stderr, " Byte Swap Data:   %ld\n", length );
#else
        FPR( stderr, " Byte Swap Data:   %d\n", length );
#endif
       }

    while( length-- )
       {
        in.f = *vec;

        out.s[0] = in.s[3];
        out.s[1] = in.s[2];
        out.s[2] = in.s[1];
        out.s[3] = in.s[0];

        *vec++ = out.f;
       }

    return( 0 );
}

/***/
/* dataClose: close file or socket at fileUnit.
/***/

int dataClose( fileUnit )

   FILE_UNIT( fileUnit );
{
    int error;

    error = 0;
    errno = 0;

    if (!fileUnit) return( 0 );

#ifdef FB_SYS_IO
    if (fileUnit > 0)
       error = close( fileUnit );
#else
    if (fileUnit)
       error = fclose( fileUnit );
#endif

    if (error)
       {
        (void) fprintf( stderr, "DATAIO File Close Error %d.\n", errno );

        if (errno) (void) perror( "SYSTEM Message" );
       }

    return( error );
}

/***/
/* dataCreate: create a binary file named "fileName".
/***/

int dataCreate( fileName )

   char *fileName;
{
    FILE_UNIT( fileUnit );
    int error;

    error = dataOpen( fileName, &fileUnit, FB_CREATE );

    if (error)
       (void) fprintf( stderr, "DATAIO File Create Error %d.\n", error );

    (void) dataClose( fileUnit );

    return( error );
}

/***/
/* dirCreate: create directory dirName.
/***/

int dirCreate( dirName )

   char *dirName;
{
    int  error;

#ifdef WIN95
    error = mkdir( dirName );
#else
    error = mkdir( dirName, 0755 );
#endif

    if (error)
       {
        (void) fprintf( stderr,
                        "DATAIO Error Creating Directory %s\n", dirName );

        (void) perror( "MKDIR Message" );
       }

    return( error );
}

/***/
/* delFile: delete a file.
/***/

int delFile( fileName )

   char *fileName;
{
    int error;

    error = 0;

#ifdef WIN95
    if (fileExists( fileName ))
       {
        (void) unlink( fileName );
       }
#else
    if (fileExists0( fileName ))
       {
        if ((error = unlink( fileName )))
           {
            (void) perror( "DATAIO Error Deleting File" );
           }
       }
#endif

    return( error );
}

#ifndef WIN95

/***/
/* delDirContents: delete the regular files in the given directory.
/***/

int delDirContents( dirName )

   char *dirName;
{
    char   fullName[NAMELEN+1], ctemp[NAMELEN+1];
    DIR    *dirUnit;
    int    i;

    struct stat   statBuff;
    struct dirent *dirInfo;

#ifdef SOLARIS
    char buff[sizeof(struct dirent) + _POSIX_PATH_MAX + NAMELEN + 1];
    struct dirent *readdir_r(); /* Not used. */
#endif

    errno   = 0;
    dirUnit = (DIR *)NULL;
    dirInfo = (struct dirent *)NULL;

    if (!(dirUnit = opendir( dirName )))
       {
        (void) fprintf( stderr, "DATAIO Error opening dir %s\n", dirName );

        if (errno) (void) perror( "System Message" );
        return( 1 );
       }

    while( (dirInfo = READDIR( dirUnit )) )
       {
        for( i = 0; i < DIR_NAMELEN; i++ )
           {
            ctemp[i] = dirInfo->d_name[i];
           }

        ctemp[DIR_NAMELEN] = '\0';

        (void) sprintf( fullName, "%s/%s", dirName, ctemp );

        if (!safe_stat( fullName, &statBuff ))
           {
            if (S_ISREG(statBuff.st_mode)) (void) delFile( fullName );
           }
       }

    (void) closedir( dirUnit );

    return( 0 );
}

#endif

/***/
/* safe_stat: ignores file size overflow.
/***/

int safe_stat( path, buff )

   char   *path;
   struct stat *buff;
{
#ifdef LINUX
    long   this_st_dev, this_st_ino, this_st_mode, this_st_nlink, this_st_uid, this_st_gid, this_st_rdev, 
           this_st_size, this_st_blksize, this_st_blocks, this_st_atime, this_st_mtime, this_st_ctime;
#endif

    int status;

    errno = 0;

#ifdef LINUX
    if ((status = stat( path, buff )))
       {
        if (getIODebugFlag()) FPR( stderr, "safe_stat Error. path: %s size: %ld status %d errno %d\n", path, (long)buff->st_size, status, errno );

#ifndef NMR64
        if (errno == EOVERFLOW)
           {
            if (getIODebugFlag()) FPR( stderr, "safe_stat: fixing overflow for path %s\n", path );

            status = safe_stat_vals( path, 
                                     &this_st_dev, &this_st_ino, &this_st_mode, &this_st_nlink, &this_st_uid, &this_st_gid, &this_st_rdev,
                                     &this_st_size, &this_st_blksize, &this_st_blocks, &this_st_atime, &this_st_mtime, &this_st_ctime );

            buff->st_dev     = this_st_dev;
            buff->st_ino     = this_st_ino;
            buff->st_mode    = this_st_mode;
            buff->st_nlink   = this_st_nlink;
            buff->st_uid     = this_st_uid;
            buff->st_gid     = this_st_gid;
            buff->st_rdev    = this_st_rdev;
            buff->st_size    = this_st_size;
            buff->st_blksize = this_st_blksize;
            buff->st_blocks  = this_st_blocks;
            buff->st_atime   = this_st_atime;
            buff->st_mtime   = this_st_mtime;
            buff->st_ctime   = this_st_ctime;

            if (status)
               {
                FPR( stderr, "safe_stat_vals failed to fix overflow for path %s errno %d\n", path, errno );
                return( status );
               }
          }
#endif
      }
#else
    status = stat( path, buff );
#endif
 
   return( status ); 
}

/***/
/* safe_fstat: ignores file size overflow.
/***/

int safe_fstat( fd, buff )

   int    fd;
   struct stat *buff;
{
#ifdef LINUX
    long   this_st_dev, this_st_ino, this_st_mode, this_st_nlink, this_st_uid, this_st_gid, this_st_rdev, 
           this_st_size, this_st_blksize, this_st_blocks, this_st_atime, this_st_mtime, this_st_ctime;
#endif

    int status;

    errno = 0;

#ifdef LINUX
    if ((status = fstat( fd, buff )))
       {
        if (getIODebugFlag()) FPR( stderr, "safe_fstat Error. fd: %d size: %ld status %d errno %d\n", fd, (long)buff->st_size, status, errno );

#ifndef NMR64
        if (errno == EOVERFLOW)
           {
            if (getIODebugFlag()) FPR( stderr, "safe_fstat: fixing overflow for fd %d\n", fd );

            status = safe_fstat_vals( fd, 
                                      &this_st_dev, &this_st_ino, &this_st_mode, &this_st_nlink, &this_st_uid, &this_st_gid, &this_st_rdev,
                                      &this_st_size, &this_st_blksize, &this_st_blocks, &this_st_atime, &this_st_mtime, &this_st_ctime );

            buff->st_dev     = this_st_dev;
            buff->st_ino     = this_st_ino;
            buff->st_mode    = this_st_mode;
            buff->st_nlink   = this_st_nlink;
            buff->st_uid     = this_st_uid;
            buff->st_gid     = this_st_gid;
            buff->st_rdev    = this_st_rdev;
            buff->st_size    = this_st_size;
            buff->st_blksize = this_st_blksize;
            buff->st_blocks  = this_st_blocks;
            buff->st_atime   = this_st_atime;
            buff->st_mtime   = this_st_mtime;
            buff->st_ctime   = this_st_ctime;

            if (status)
               {
                FPR( stderr, "safe_fstat_vals failed to fix overflow for fd %d errno %d\n", fd, errno );
                return( status );
               }
          }
#endif
      }
#else
    status = fstat( fd, buff );
#endif
  
   return( status ); 
}

#ifdef NMR64
int dummy_dataio_safe_stat_vals()
{
   return( 0 );
}
#else
#define _FILE_OFFSET_BITS 64

/***/
/* safe_stat_vals: transfer stat64() results for 32-bit use.
/***/

int safe_stat_vals( path, this_st_dev, this_st_ino, this_st_mode, this_st_nlink, this_st_uid, this_st_gid, this_st_rdev,
                    this_st_size, this_st_blksize, this_st_blocks, this_st_atime, this_st_mtime, this_st_ctime )

   char   *path;
   long   *this_st_dev, *this_st_ino, *this_st_mode, *this_st_nlink, *this_st_uid, *this_st_gid, *this_st_rdev, 
          *this_st_size, *this_st_blksize, *this_st_blocks, *this_st_atime, *this_st_mtime, *this_st_ctime;
{
   struct stat buff;
   int    status;

   errno = 0;

   status = stat( path, &buff );

   if (status)
      {
       if (errno == EOVERFLOW)
          {
           errno  = 0;
           status = 0;
          }
      }

   *this_st_dev     = (long)buff.st_dev;
   *this_st_ino     = (long)buff.st_ino;
   *this_st_mode    = (long)buff.st_mode;
   *this_st_nlink   = (long)buff.st_nlink;
   *this_st_uid     = (long)buff.st_uid;
   *this_st_gid     = (long)buff.st_gid;
   *this_st_rdev    = (long)buff.st_rdev;
   *this_st_size    = (long)buff.st_size;
   *this_st_blksize = (long)buff.st_blksize;
   *this_st_blocks  = (long)buff.st_blocks;
   *this_st_atime   = (long)buff.st_atime;
   *this_st_mtime   = (long)buff.st_mtime;
   *this_st_ctime   = (long)buff.st_ctime;
    
   return( status ); 
}

/***/
/* safe_fstat_vals: transfer fstat64() values for 32-bit use.
/***/

int safe_fstat_vals( fd, this_st_dev, this_st_ino, this_st_mode, this_st_nlink, this_st_uid, this_st_gid, this_st_rdev,
                     this_st_size, this_st_blksize, this_st_blocks, this_st_atime, this_st_mtime, this_st_ctime )

   int    *fd;
   long   *this_st_dev, *this_st_ino, *this_st_mode, *this_st_nlink, *this_st_uid, *this_st_gid, *this_st_rdev,
          *this_st_size, *this_st_blksize, *this_st_blocks, *this_st_atime, *this_st_mtime, *this_st_ctime;

{
   struct stat buff;
   int    status;

   errno = 0;

   status = fstat( fd, &buff );

   if (status) 
      { 
       if (errno == EOVERFLOW) 
          { 
           errno  = 0;
           status = 0;
          }
      }

   *this_st_dev     = (long)buff.st_dev;
   *this_st_ino     = (long)buff.st_ino;
   *this_st_mode    = (long)buff.st_mode;
   *this_st_nlink   = (long)buff.st_nlink;
   *this_st_uid     = (long)buff.st_uid;
   *this_st_gid     = (long)buff.st_gid;
   *this_st_rdev    = (long)buff.st_rdev;
   *this_st_size    = (long)buff.st_size;
   *this_st_blksize = (long)buff.st_blksize;
   *this_st_blocks  = (long)buff.st_blocks;
   *this_st_atime   = (long)buff.st_atime;
   *this_st_mtime   = (long)buff.st_mtime;
   *this_st_ctime   = (long)buff.st_ctime;

   return( status );
}
#endif

/***/
/* Bottom
/***/
