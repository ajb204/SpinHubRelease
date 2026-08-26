/***/
/* Include file for string tokenizing functions:
/***/

#include "prec.h"

    int   badToken;

    char    *getToken();
    char    *getTokenS();
    char    *getTokenC();

    float   getTokenR();
    int     getTokenI();
    NMR_INT getTokenN();
    
    char  *leftToken();
    char  *rightToken();

    int   getTokenLoc();
    int   cntToken();
    int   cntTokenS();
    int   cntTokenC();
    int   freeToken();
    int   getTokenMaxLength();

    int   findToken();

    int   isInteger();
    int   isFloat();

    int   istrlen();

    int   str2hash();
