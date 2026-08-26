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
/* raiseStr: change case of lower-case characters in string.
/*
/* Adjusted 6/24/92 for possible bug in INTEL toupper library function.
/***/

int raiseStr( string )

   char *string;
{
   if (!string) return( 0 );

   while( *string )
      {
       if (*string >= 'a' && *string <= 'z') *string -= ('a' - 'A');
       string++;
      }

   return( 0 );
}
