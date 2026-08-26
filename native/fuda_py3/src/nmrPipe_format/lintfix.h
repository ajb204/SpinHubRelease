
#ifndef _nmr_lintfix_h
#define _nmr_lintfix_h

#ifdef S_SPLINT_S
   int strcasecmp( const char *s1, const char *s2 );

   FILE *popen(const char *command, const char *type);
   FILE *fdopen(int fd, const char *mode);

   int fileno(FILE *stream);

   ssize_t recv(int sockfd, void *buf, size_t len, int flags);
   ssize_t send(int sockfd, const void *buf, size_t len, int flags);

   int wait(int *status);

   int re_exec(char *string);

#ifndef ENOENT
#define	ENOENT 2 /* No such file or directory */
#endif

#endif

#endif

