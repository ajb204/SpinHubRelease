/***/
/* byteSwap2: perform 2-byte swap on words in buff.
/***/

#include "prec.h"
#include "bswap.h"

int byteSwap2_64( buff, length )

   char    *buff;
   NMR_INT length;
{
    char c0, c1;

    if (!buff || length < 1) return( 0 );

    while( length > 1 )
       {
        c0 = buff[0];
        c1 = buff[1];

        buff[0] = c1;
        buff[1] = c0;

        length -= 2;
        buff   += 2;
       }

    return( 0 );
}

/***/
/* byteSwap3: perform 3-byte swap on words in buff.
/***/

int byteSwap3_64( buff, length )

   char    *buff;
   NMR_INT length;
{
    char c0, c1, c2;

    if (!buff || length < 1) return( 0 );

    while( length > 2 )
       {
        c0 = buff[0];
        c1 = buff[1];
        c2 = buff[2];

        buff[0] = c2;
        buff[1] = c1;
        buff[2] = c0;

        length -= 3;
        buff   += 3;
       }

    return( 0 );
}

/***/
/* byteSwap4: perform 4-byte swap on words in buff.
/***/

int byteSwap4_64( buff, length )

   char    *buff;
   NMR_INT length;
{
    char c0, c1, c2, c3;
 
    if (!buff || length < 1) return( 0 );

    while( length > 3 )
       {
        c0 = buff[0];
        c1 = buff[1];
        c2 = buff[2];
        c3 = buff[3];

        buff[0] = c3;
        buff[1] = c2;
        buff[2] = c1;
        buff[3] = c0;

        length -= 4; 
        buff   += 4;
       }

    return( 0 );
}

/***/
/* byteSwap8: perform 8-byte swap on words in buff.
/***/

int byteSwap8_64( buff, length )

   char    *buff;
   NMR_INT length;
{
    char c0, c1, c2, c3, c4, c5, c6, c7;
 
    if (!buff || length < 1) return( 0 );

    while( length > 7 )
       {
        c0 = buff[0];
        c1 = buff[1];
        c2 = buff[2];
        c3 = buff[3];
        c4 = buff[4];
        c5 = buff[5];
        c6 = buff[6];
        c7 = buff[7];

        buff[0] = c7;
        buff[1] = c6;
        buff[2] = c5;
        buff[3] = c4;
        buff[4] = c3;
        buff[5] = c2;
        buff[6] = c1;
        buff[7] = c0;

        length -= 8; 
        buff   += 8;
       }

    return( 0 );
}

/***/
/* byteSwapN: perform N-byte swap of data in buff.
/***/

int byteSwapN_64( buff, length, n )

   char    *buff;
   NMR_INT length;
   int     n;
{
    char ctemp, *head, *tail;
    int  i, count;

    if (!buff || length < 1 || n < 2) return( 0 );

    if (n == 2)
       {
        (void) byteSwap2( buff, length );
	return( 0 );
       }
    else if (n == 3)
       {
        (void) byteSwap3( buff, length );
	return( 0 );
       }
    else if (n == 4)
       {
        (void) byteSwap4( buff, length );
	return( 0 );
       }
    else if (n == 8)
       {
        (void) byteSwap8( buff, length );
	return( 0 );
       }

    count = n/2;

    while( length >= n )
       {
        head = buff;
        tail = buff + n - 1;

        for( i = 0; i < count; i++ )
           {
            ctemp = *head; *head = *tail; *tail = ctemp;
            head++;
            tail--;
           }

        length -= n;
        buff   += n;
       }
 
    return( 0 );
}

int vDbl2Flt_64( s, n )

   char    *s;
   NMR_INT n;
{
    union{ char c[8]; double d; } ud;
    union{ char c[4]; float  f; } uf;

    char *src, *dest;
    int  i;

    if (!s || n < 1) return( 0 );

    src  = s;
    dest = s;

    for( i = 0; i < n; i++ )
       {
        ud.c[0] = *src++;
        ud.c[1] = *src++;
        ud.c[2] = *src++;
        ud.c[3] = *src++;
        ud.c[4] = *src++;
        ud.c[5] = *src++;
        ud.c[6] = *src++;
        ud.c[7] = *src++;

	uf.f = (float)ud.d;

	*dest++ = uf.c[0];
	*dest++ = uf.c[1];
	*dest++ = uf.c[2];
	*dest++ = uf.c[3];
       }

    return( 0 );
}
