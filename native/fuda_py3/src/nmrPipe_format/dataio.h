/***/
/* Include file for C and system I/O routines which
/* can be selected via conditional compiling.
/*
/*   FB_SYS_IO = 1 System I/O (default)
/*   FB_SYS_IO = 0 C I/O
/***/

#ifndef FB_SYS_IO
#define FB_SYS_IO 1
#endif

#if defined (LINUX) || defined(MAC_OSX) || defined(WINXP)
#define DATAIO_TIMEOUT 0 
#else
#define DATAIO_TIMEOUT 512
#endif

#include <unistd.h>

#include "memory.h"
#include "prec.h"

int showDebugIO();

/***/
/* Defines for manipulating directory files:
/***/

#ifdef SGI
#define DIR_NAMELEN (strlen(dirInfo->d_name))
#endif

#ifdef CRAY
#define DIR_NAMELEN (strlen(dirInfo->d_name))
#endif

#if defined(LINUX)
#define DIR_NAMELEN (dirInfo->d_reclen)
#endif

#ifdef SOLARIS
#define DIR_NAMELEN (strlen(dirInfo->d_name))
#define OLD_READDIR(D)  readdir_r( D, (struct direct *)buff )
#define READDIR(D)      readdir( D )
#endif

#ifndef DIR_NAMELEN
#define DIR_NAMELEN (dirInfo->d_namlen)
#endif

#ifndef READDIR
#define READDIR(D)  readdir(D)
#endif

/***/
/* Defines for I/O by system calls:
/***/

#ifdef FB_SYS_IO

#define USLEEP( T ) if ((T)!=0) (void) uSleep( (unsigned int)(T) )
                    
#define FILE_UNIT( UNIT ) int UNIT
#define UNIT_CAST( UNIT ) ((int) (UNIT))
#define UNIT_NULL 0

#define FB_STDIN  0
#define FB_STDOUT 1
#define FB_STDERR 2
typedef int UNIT;

/***/
/* Defines for I/O through C library:
/***/

#else

#define USLEEP( T ) (if ((T)!=0) (void) uSleep( (unsigned int)(T) ))

#define FILE_UNIT( UNIT ) FILE *UNIT
#define UNIT_CAST( UNIT ) ((FILE *) (UNIT))
#define FB_STDIN  stdin
#define UNIT_NULL (FILE *)NULL 

#define FB_STDOUT stdout
#define FB_STDERR stderr
typedef FILE *UNIT;
#endif

#define FB_READ        1
#define FB_WRITE       2
#define FB_READWRITE   3
#define FB_CREATE      4
#define FB_CREATE_KEEP 5
#define FB_UNKNOWN     6

#define dataRead(   FD, BUFF, N )          DataRead(   FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N) )
#define dataWrite(  FD, BUFF, N )          DataWrite(  FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N) )
#define dataReadB(  FD, BUFF, N, MT, TO )  DataReadB(  FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )
#define dataWriteB( FD, BUFF, N, MT, TO )  DataWriteB( FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )
#define dataSendB(  FD, BUFF, N, MT, TO )  DataSendB(  FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )
#define dataRecvB(  FD, BUFF, N, MT, TO )  DataRecvB(  FD, ((float *)(BUFF)), (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )
#define dataSendCB( FD, BUFF, N, MT, TO )  DataSendCB( FD, ((char *)(BUFF)),  (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )
#define dataRecvCB( FD, BUFF, N, MT, TO )  DataRecvCB( FD, ((char *)(BUFF)),  (NMR_INT)((NMR_INT)N), ((int)(MT)), ((int)(TO)) )

#define dataPos(  FD, POS )           DataPos(  FD, (NMR_INT)((NMR_INT)POS) )
#define dataPos2( FD, POS, BUFF, N )  DataPos2( FD, (NMR_INT)((NMR_INT)POS), ((VOID *)(BUFF)), (NMR_INT)((NMR_INT)N) )

#define dataTrunc( NAME, N, OPTR, V ) DataTrunc( NAME, (NMR_INT)((NMR_INT)N), OPTR, V )
#define dataTruncU( UNIT, N )         DataTruncU( UNIT, (NMR_INT)((NMR_INT)N) )

#define ByteSwapV( BUFF, N ) ByteSwapV( BUFF, (NMR_INT)((NMR_INT)N) )

int initDataIO();
int DataRead(), DataWrite(), DataReadB(), DataWriteB();
int DataSendB(), DataRecvB(), DataSendCB(), DataRecvCB();
int DataPos(), DataPos2();
int DataTrunc();

int dataOpen(), dataOpen2(), dataPOpen(), dataReadC(), dataWriteC(), dataReadCB();
int dataClose(), dataCreate(), dirCreate();

int getAutoSwapFlag(), setAutoSwapFlag(), getByteSwapFlag(), setByteSwapFlag();
int getDefTimeout(), getFSCheck();
int delFile(), delDirContents();

int uSleep();

int byteSwapV();


//int safe_stat( char* path, struct stat *buff );
//int dataWriteCB( FILE_UNIT(outUnit), char *array, NMR_INT bytesRequested, int maxTries, int timeOut );



  
