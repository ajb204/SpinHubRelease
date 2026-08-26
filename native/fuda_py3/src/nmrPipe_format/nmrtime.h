/***/
/* nmrTime.h definitions for conditional performance timing of nmrPipe.
/***/

#ifdef NMRTIME

#include <sys/time.h>
#include <sys/timeb.h>
#include <values.h>

struct Timer
   {
    long  sec;
    unsigned short msec;
    float total; 
   };

static struct timeb tp;

#define TIMER_ON( T )   (void)ftime(&tp); \
                        T.sec=tp.time; T.msec=tp.millitm

#define TIMER_OFF( T )  (void)ftime(&tp);          \
                        T.total+=(tp.time-T.sec) + \
                        0.001*((float)tp.millitm-(float)T.msec)

#define TIMER_INCR( I ) if (I < MAXINT) I++

struct Timer procTime, readTime, writeTime, readWTime, writeWTime,
             clientReadTime, clientWriteTime, clientReadWTime,
             clientWriteWTime, connectTime, systemTime, msgTime,
             totalTime;

struct Timer serverReadTime[MAXSOCK], serverWriteTime[MAXSOCK],
             serverReadWTime[MAXSOCK], serverWriteWTime[MAXSOCK];

int          readRetries, writeRetries, sReadRetries, sWriteRetries;

#else

#define TIMER_ON( T )   /* No timer ON.        */
#define TIMER_OFF( T )  /* No timer OFF.       */
#define TIMER_INCR( I ) /* No timer increment. */

#endif
