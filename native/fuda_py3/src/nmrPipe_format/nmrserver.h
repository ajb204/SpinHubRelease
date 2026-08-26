/***/
/* nmrserver.h: definitions for client/server processing.
/***/

#define MAXSOCK 64

typedef struct ServerInfo
   {
     char  *name;            /* Host name of remote server.                   */
     char  *cmnd;            /* Special command to start server, if any.      */
     int   cpuID;            /* CPU number (for multi-cpu machines).          */
     int   pid;              /* Process ID of server on remote machine.       */
     int   sock;             /* Client's socket to server.                    */
     int   partition;        /* Slice count for server to process.            */
     int   capacity;         /* Buffer capacity for server.                   */
     float speed;            /* Relative speed of server (higher is faster).  */
   } Server;

int freeServers();
int getServerLoc();
int getServerRange();
int nullServer();
int rdServers();
int setCapacity();
int setPartition();
int showServers();

