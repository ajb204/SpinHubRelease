#define DTABLESIZE getdtablesize()

#ifdef IBM
#include <sys/select.h>
#endif

#ifdef SOLARIS
#include <limits.h>
#undef  DTABLESIZE
#define DTABLESIZE OPEN_MAX
#endif

#ifdef WIN95 
#include <limits.h>
#undef  DTABLESIZE
#define DTABLESIZE 256
#endif

#ifdef WINNT
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/select.h>
#undef  DTABLESIZE
#define DTABLESIZE FD_SETSIZE
#endif

#ifdef HP 
#include <sys/stat.h>
#undef  DTABLESIZE
#define DTABLESIZE (size_t) 256
#endif

#ifdef LINUX
#include <sys/stat.h>
#endif

#ifdef WINXP 
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef CONVEX
#include <sys/stat.h>
#endif

#ifdef ALPHA
#include <sys/stat.h>
#endif

#ifdef SGI
#include <sys/stat.h>
#endif

#ifdef SOLARIS
#include <sys/stat.h>
#endif

#ifdef CRAY
#include <sys/stat.h>
#endif

static  fd_set readfds, writefds;

#ifndef MAC_OSX_NEW
//void *(void *, int, unsigned long, unsigned long) //copied from sdk
//I don't think Frank ever uses this? It won't compile on High Sierra as is 14thSept2022
//void    bZero(void *, int, unsigned long, unsigned long), bzero(void *, int, unsigned long, unsigned long);
#endif

#define FDNULL (fd_set *) NULL
#define TONULL (struct timeval *) NULL

#if defined (SOLARIS) || defined (WIN95)
#define FD_ZERO2(p) bZero((char *)(p), (int)sizeof(*(p)))
#else
#define FD_ZERO2(p) bzero((VOID *)(p), (int)sizeof(*(p)))
#endif
