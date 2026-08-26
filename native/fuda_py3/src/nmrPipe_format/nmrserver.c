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
/* nmrserver: utilities for parallel processing services.
/*
/* rdServers:      read name and speed info from server list.
/* setCapacity:    set capacities for a list of servers. 
/* setPartition:   set partition sizes for a list of servers. 
/* getServerLoc:   get a server index given the server name.
/* getServerRange: get a partition range for a server given a server list.
/* freeServers:    deallocate server list contents. 
/***/

#include <stdio.h>
#include <string.h>

#include "nmrserver.h"
#include "memory.h"
#include "token.h"
#include "prec.h"
#include "rdtext.h"

#define FPR (void) fprintf

void gethostname();

int rdServers( serverListName, serverList, maxServers, serverCount )

    int    *serverCount, maxServers;
    char   *serverListName;
    Server *serverList;
{
    int   i, j, lineCount, argc, error;
    char  **nameList, serverName[NAMELEN+1];
    float serverSpeed;

/***/
/* Initialize.
/* Stop if no list name is given.
/* Read list of server names.
/***/

    error        = 0;
    nameList     = (char **) NULL;
    *serverCount = 0;

    for( i = 0; i < maxServers; i++ )
       {
        (void) nullServer( &serverList[i] );
       }

    if (!*serverListName) return( 0 );

    if (!(nameList = rdText( serverListName, &lineCount, &error )))
       {
        FPR( stderr,
             "Server Error reading server list %s\n",
             serverListName );

        goto shutdown;
       }

/***/
/* Scan list for server names and speeds:
/***/

    for( i = 0; i < lineCount; i++ )
       {
        argc = cntToken( nameList[i] );

        if (nameList[i][0] != '#' && argc >= 2)
           {
            if (2 != sscanf( nameList[i], "%s %f", serverName, &serverSpeed ))
               {
                FPR( stderr, "Server Error in Server List.\n" );

                error = 1;
                goto shutdown;
               }

            if (serverSpeed > 0.0)
               {
                serverList[*serverCount].name  = strDup( serverName );
                serverList[*serverCount].speed = serverSpeed;

                if (!serverList[*serverCount].name)
                   {
                    FPR( stderr, "Server Name Memory Error.\n" );

                    error = 1;
                    goto shutdown;
                   }

                if (argc > 2)
                   {
                    serverList[*serverCount].cmnd = 
                       strDup( rightToken( 3, nameList[i] ));

                    if (!serverList[*serverCount].cmnd)
                       {
                        FPR( stderr, "Server Command Memory Error.\n" );
    
                        error = 1;
                        goto shutdown;
                       }
                   }

                (*serverCount)++;

                if (*serverCount >= maxServers)
                   {
                    FPR( stderr, "Server Error in List Length.\n" );

                    error = 1;
                    goto shutdown;
                   }
               }
           }
       }

/***/
/* Test for valid server info.
/* Install CPU ID for server names which appear more than once. 
/***/

    if (!*serverCount)
       {
        FPR( stderr, "Server Error: Null Server List.\n" );
        error = 1;
        goto shutdown;
       }

    for( i = 0; i < *serverCount - 1; i++ )
       {
        for( j = i + 1; j < *serverCount; j++ )
           {
            if (!strcmp( serverList[i].name, serverList[j].name ))
               {
                serverList[j].cpuID++;
               }
           }
       }

/***/
/* Exit point:
/***/

shutdown:

    (void) freeText( nameList );

    return( error );
}

/***/
/* nullServer: NULL-initialize elements in a server entry.
/***/

int nullServer( serverInfo )

   Server *serverInfo;
{
#ifndef S_SPLINT_S
    serverInfo->name      = (char *) NULL;
    serverInfo->cmnd      = (char *) NULL;
#endif
    serverInfo->pid       = 0;
    serverInfo->cpuID     = 0;
    serverInfo->capacity  = 0;
    serverInfo->partition = 0;
    serverInfo->speed     = 0.0;

    return( 0 );
}

/***/
/* setCapacity: set buffer capacity for a list of servers.
/***/

int setCapacity( serverList, serverCount, capacity )

   int serverCount, capacity;
   Server *serverList;
{
    int   i;
    float minSpeed;

/***/
/* Find minimum speed.
/* Abort if minimum speed is no good.
/* Set capacities proportionally to speed.
/***/

    minSpeed = serverList[0].speed;

    for( i = 0; i < serverCount; i++ )
       {
        if (minSpeed > serverList[i].speed) minSpeed = serverList[i].speed;
       }

    if (minSpeed <= 0.0)
       {
        FPR( stderr, "Server Error in Server Speeds.\n" );
        return( 1 );
       }
   
    for( i = 0; i < serverCount; i++ )
       {
        serverList[i].capacity = 
           (float)capacity*serverList[i].speed/minSpeed; 

        if (serverList[i].capacity < 1) serverList[i].capacity = 1;
       }

    return( 0 );
}

/***/
/* setPartion: set partitions for a list of servers.
/***/

int setPartition( serverList, serverCount, itemCount )

   int serverCount, itemCount;
   Server *serverList;
{
    int   i, itemsUsed;
    float totalSpeed;

/***/
/* Find total speed.
/* Set partition size proportionally to relative speed. 
/***/

    itemsUsed  = 0;
    totalSpeed = 0.0;

    for( i = 0; i < serverCount; i++ )
       {
        totalSpeed += serverList[i].speed;
       }

    if (totalSpeed <= 0.0)
       {
        FPR( stderr, "Server Error in Server Speeds.\n" );
        return( 1 );
       }

    for( i = 0; i < serverCount; i++ )
       {
        serverList[i].partition =
               (float)itemCount*serverList[i].speed/totalSpeed;

        itemsUsed += serverList[i].partition;
       }

/***/
/* Adjust so that total of partition sizes equals number of items:
/***/
   
    if (itemsUsed > itemCount)
       {
        serverList[0].partition -= itemsUsed - itemCount;

        if (serverList[0].partition < 0)
           {
            FPR( stderr, "Server Error adjusting partitions.\n" );
            return( 1 );
           }
       }
    else if (itemsUsed < itemCount)
       {
        serverList[0].partition += itemCount - itemsUsed;
       }

    return( 0 );
}

/***/
/* getServerLoc: get a server index given the server name.
/*               If serverName is a null pointer, current host name is used.
/***/

int getServerLoc( serverList, serverCount, 
                  serverName, cpuID, 
                  serverLoc, verbose )

    int    *serverLoc, serverCount, cpuID, verbose;
    Server *serverList;
    char   *serverName;
{
    int i;
    char hostName[NAMELEN+1];

#ifdef WIN95
    if (serverName)
       (void) strcpy( hostName, serverName ); 
    else
       (void) strcpy( hostName, "windows95" );
#else
    if (serverName)
       (void) strcpy( hostName, serverName ); 
    else
       (void) gethostname( hostName, NAMELEN );
#endif

    if (verbose)
       {
        FPR( stderr, "Searching for host %s CPU %d:\n", hostName, cpuID );
       }

    for( i = 0; i < serverCount; i++ )
       {
        if (verbose)
           {
            FPR( stderr, 
                 " Testing host %s CPU %d:\n", 
                 serverList[i].name, serverList[i].cpuID );
           }

        if (!strcmp( serverList[i].name, hostName ) && 
              serverList[i].cpuID == cpuID)
           {
            if (verbose)
               {
                FPR( stderr, 
                     " Identified host %s CPU %d at Location %d:\n", 
                     serverList[i].name, serverList[i].cpuID, i );
               }

            *serverLoc = i;
            return( 0 );
           }
       }

    return( 1 );
}

/***/
/* getServerRange: get an origin=1 partition range for a server.
/***/

int getServerRange( serverList, serverLoc, first, last )

    int    serverLoc, *first, *last;
    Server *serverList;
{
    int i;

    *last = 0;

    for( i = 0; i <= serverLoc; i++ )
       {
        *first  =  *last + 1;
        *last   =  *first + serverList[i].partition - 1;
       }

    return( 0 );
}

/***/
/* showServers: display information about the listed servers.
/***/

int showServers( serverList, serverCount )

   Server *serverList;
   int    serverCount;
{
    int i;

    for( i = 0; i < serverCount; i++ )
       {
        FPR( stderr,
             "Host %2d: %10s CPU=%d Partition=%d Speed=%g\n",
             i+1, serverList[i].name,
             serverList[i].cpuID,
             serverList[i].partition,
             serverList[i].speed );
       }

    return( 0 );
}

/***/
/* freeServers: deallocate data associated with servers.
/***/

int freeServers( serverList, serverCount )

   Server *serverList;
   int    serverCount;
{
    int i;

    for( i = 0; i < serverCount; i++ )
       {
        if (serverList[i].name)
           {
            (void) deAlloc( "server",
                            serverList[i].name,
                            1+strlen(serverList[i].name) );
           }

        if (serverList[i].cmnd)
           {
            (void) deAlloc( "server",
                            serverList[i].cmnd,
                            1+strlen(serverList[i].cmnd) );
           }

        (void) nullServer( &serverList[i] );
       }

    return( 0 );
}

/***/
/* Bottom.
/***/
