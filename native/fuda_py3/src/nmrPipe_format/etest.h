/***/
/* etest.h: macro for abrieviated error checking.
/***/

#ifdef ETEST_PROG_NAME
#define ETEST( ID, EXPRESSION )                            \
   if ((error = (EXPRESSION)))                             \
      {                                                    \
       (void) fprintf( stderr, "%s ETEST error %d.\n", ETEST_PROG_NAME, ID );  \
       goto shutdown;                                      \
      }
#else
#define ETEST( ID, EXPRESSION )                            \
   if ((error = (EXPRESSION)))                             \
      {                                                    \
       (void) fprintf( stderr, "ETEST error %d.\n", ID );  \
       goto shutdown;                                      \
      }
#endif
