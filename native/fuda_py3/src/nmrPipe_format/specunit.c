/***/
/* NMRPipe System (C) Frank Delaglio 1195-2016.
/*
/* specunit.c: utilities for parsing spectral units.
/*
/* spec2rPnt: converts a location in spectral units to a real point index.
/* spec2iPnt: converts a location in spectral units to an integer point index.
/* iPnt2spec: converts an integer point index to the desired spectral unit.
/* rPnt2spec: converts a real point index to the desired spectral unit.
/*
/* specWidth2rPnt: converts a spectral width to a width in points.
/* rPnt2specWidth: converts a point width to a width in spectral units. 
/*
/* specList2Pnt: allocates a list of integer point locations corresponding to
/*               a list of spectral locations.
/*
/* updateOrigin: recalculates spectral origin on basis of carrier frequency.
/* hasSpecLabel: returns 1 if the given text has a valid unit label.
/*
/* Nota Bene: 1. Integer point index values vary from 1 to NDSIZE.
/*            2. Global variable "badSpecUnits" is set to indicate errors.
/*
/* Conversion formulas by S. Grzesiek and D. Garrett.
/***/

#include <string.h>
#include <stdio.h>

#include "cmndargs.h"
#include "specunit.h"
#include "memory.h"
#include "fdatap.h"
#include "ctype.h"
#include "prec.h"
#include "atof.h"

/***/
/* spec2rPnt: convert spectral units to real point index.
/***/

int raiseStr();
int hasUnitLabel();
int str2SpecVal();
int isSpatialDim();
int str2SpecVal();
int isSpatialDim();
int str2LabVal();

float spec2rPnt( fdata, dimCode, specLabel )

   float fdata[FDATASIZE];
   char  *specLabel;
   int   dimCode;
{
    float value;

    value = spec2rPnt2( fdata, dimCode, specLabel, CLIP_FOLD );

    return( value );
}

/***/
/* spec2rPnt2: convert spectral units to real point index.
/***/

float spec2rPnt2( fdata, dimCode, specLabel, clipMode )

   float fdata[FDATASIZE];
   int   dimCode, clipMode; 
   char  *specLabel;
{
    float acqTime, specVal, pntVal, hzVal, df, sw, uw, obs, orig, car, s1, sN, delta, first;
    int   ilFlag, trueSize, size;
    char  label[NAMELEN+1];

/***/
/* Initialize.
/* Extract header parameters for this axis.
/* Extract elements of spectral label.
/***/

    badSpecUnits = 0;
    specVal      = 0.0;
    pntVal       = 0.0;
    *label       = '\0';

    if (!strcasecmp( specLabel, "df" ))
       {
        df = getParm( fdata, FDDMXVAL, NULL_DIM ); 
        if (df < 0.0) df = -df;
        return( df );
       }

    sw     = getParm( fdata, NDSW,     dimCode );
    obs    = getParm( fdata, NDOBS,    dimCode );
    car    = getParm( fdata, NDCAR,    dimCode );
    orig   = getParm( fdata, NDORIG,   dimCode );
       
    if (sw  == 0.0) sw  = 1.0;
    if (obs == 0.0) obs = 1.0;

    if (orig > car)
       {
        s1 = orig - sw;
        sN = orig;
       }
    else
       {
        s1 = orig + sw;
        sN = orig;
       }

    uw = sN - s1;

    (void) str2SpecVal( specLabel, &specVal, label );

/***/
/* Adjust size if this dimension has interleaved complex data.
/* Compute first hertz and hertz/point ratios accordingly.
/***/

    ilFlag   = isInterleaved( fdata, dimCode );
    trueSize = getParm( fdata, NDSIZE, dimCode );
    size     = trueSize;

    if (ilFlag) size /= 2;

    delta   = -sw/(float)size;
    first   = orig - delta*(size - 1);
    acqTime = size/sw;

/***/
/* Handle special case for NM, wavelength or spatial distance.
/* Identify spectral units.
/* Perform conversion.
/***/

    if (!strcasecmp( label, LAB_NM ))
       {
        if (!isSpatialDim( fdata, dimCode )) (void) strcpy( label, LAB_WNM );
       }

    if (!strcasecmp( label, LAB_PPM ))
       {
        pntVal   = 1.0 + (specVal*obs - first)/delta;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > size ) pntVal -= size;
            while( pntVal < 1 )    pntVal += size;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)    pntVal = 1;
            if (pntVal > size) pntVal = size;
           }

        if (ilFlag) pntVal = 2.0*pntVal - 1.0;
       }
    else if (!strcasecmp( label, LAB_HZ ))
       {
        pntVal   = 1.0 + (specVal - first)/delta;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > size ) pntVal -= size;
            while( pntVal < 1 )    pntVal += size;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)    pntVal = 1;
            if (pntVal > size) pntVal = size;
           }
 
        if (ilFlag) pntVal = 2.0*pntVal - 1.0;
       }
    else if (!strcasecmp( label, LAB_WN ))
       {
        pntVal   = 1.0 + (specVal*C_CM_SEC - first)/delta;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > size ) pntVal -= size;
            while( pntVal < 1 )    pntVal += size;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)    pntVal = 1;
            if (pntVal > size) pntVal = size;
           }

        if (ilFlag) pntVal = 2.0*pntVal - 1.0;
       }
    else if (!strcasecmp( label, LAB_WNM ))
       {
        sw     = C_NM_SEC/sw;
        orig   = C_NM_SEC/orig - sw;

        pntVal = 1.0 + (size - 1)*(specVal - orig)/sw;
 
        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > size ) pntVal -= size;
            while( pntVal < 1 )    pntVal += size;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)    pntVal = 1;
            if (pntVal > size) pntVal = size;
           }

        if (ilFlag) pntVal = 2.0*pntVal - 1.0;
       }
    else if (!strcasecmp( label, LAB_PCT ))
       {
        pntVal = 1.0 + specVal*(trueSize - 1)/100.0;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_SEC ))
       {
        pntVal = acqTime == 0.0 ? 0.0 : 1.0 + specVal*(size - 1)/acqTime;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_CM ))
       {
        pntVal = 1.0 + specVal*(trueSize - 1)/(0.1*sw);

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_MM ))
       {
        pntVal = 1.0 + specVal*(trueSize - 1)/sw;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_MIC ))
       {
        pntVal = 1.0 + specVal*(trueSize - 1)/(1000.0*sw);

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_NM ))
       {
        pntVal = 1.0 + specVal*(trueSize - 1)/(1000000.0*sw);

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!strcasecmp( label, LAB_USER ))
       {
        pntVal = 1.0 + (specVal - s1)*(trueSize - 1)/uw;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal < 1)        pntVal = 1;
            if (pntVal > trueSize) pntVal = trueSize;
           }
       }
    else if (!*label || !strcasecmp( label, LAB_PTS ))
       {
        pntVal = specVal;

        if (clipMode == CLIP_FOLD)
           {
            while( pntVal > trueSize ) pntVal -= trueSize;
            while( pntVal < 1 )        pntVal += trueSize;
           }
        else if (clipMode == CLIP_TRUNC)
           {
            if (pntVal == 0)        pntVal = 1;
            if (pntVal <  0)        pntVal = trueSize + pntVal;
            if (pntVal <  0)        pntVal = 1;
            if (pntVal >  trueSize) pntVal = trueSize;
           }
       }
    else
       {
        badSpecUnits = 1;
       }

    return( pntVal );
}

/***/
/* spec2iPnt: convert spectral units to integer point index.
/***/

int spec2iPnt( fdata, dimCode, specLabel )

   float fdata[FDATASIZE];
   char  *specLabel;
   int   dimCode; 
{
    float rPoint;
    int   pntVal;

    rPoint = spec2rPnt( fdata, dimCode, specLabel );
    pntVal = rPoint > 0.0 ? rPoint + 0.5 : rPoint - 0.5;

    return( pntVal );
}

/***/
/* iPnt2spec: convert an integer point index to a location in spectral units.
/***/

float iPnt2spec( fdata, dimCode, pntVal, specLabel )

   float fdata[FDATASIZE];
   char  *specLabel;
   int   pntVal, dimCode; 
{
    float rPntVal, specVal;

    rPntVal = pntVal;
    specVal = rPnt2spec( fdata, dimCode, rPntVal, specLabel );

    return( specVal );
}

/***/
/* rPnt2spec: convert a real point index to a location in spectral units.
/***/

float rPnt2spec( fdata, dimCode, pntVal, label )

   float pntVal, fdata[FDATASIZE];
   char  *label;
   int   dimCode; 
{
    float acqTime, specVal, hzVal, sw, uw, orig, obs, car, delta, first, s1, sN, oldPntVal;
    int   trueSize, size;

    char  specLabel[NAMELEN+1];

/***/
/* Initialize.
/* Extract header parameters for this axis.
/***/

    (void) strcpy( specLabel, label );

    badSpecUnits = 0;
    specVal      = 0.0;

    sw     = getParm( fdata, NDSW,     dimCode );
    obs    = getParm( fdata, NDOBS,    dimCode );
    car    = getParm( fdata, NDCAR,    dimCode );
    orig   = getParm( fdata, NDORIG,   dimCode );

    if (sw == 0.0)  sw  = 1.0;
    if (obs == 0.0) obs = 1.0;

    if (orig > car)
       {
        s1 = orig - sw;
        sN = orig;
       }
    else
       {
        s1 = orig + sw;
        sN = orig;
       }

    uw = sN - s1;

/***/
/* Handle special case for NM, wavelength or spatial distance.
/* Adjust size and point value if this dimension has interleaved complex data.
/* Compute first hertz and hertz/point ratios accordingly.
/***/

    if (!strcasecmp( specLabel, LAB_NM ))
       {
        if (!isSpatialDim( fdata, dimCode )) (void)strcpy( specLabel, LAB_WNM );
       }

    trueSize  = getParm( fdata, NDSIZE, dimCode );
    size      = trueSize;
    oldPntVal = pntVal;

    if (isInterleaved( fdata, dimCode ))
       {
        if ((int)pntVal % 2)
           pntVal  = (pntVal + 1.0)/2.0; 
        else
           pntVal  = pntVal/2.0;

        size  /= 2;
       }

    delta   = -sw/(float)size;
    first   = orig - delta*(size - 1);
    acqTime = size/sw;

/***/
/* Identify spectral units.
/* Perform conversion.
/***/

    if (!strcasecmp( specLabel, LAB_PPM ))
       {
        specVal  = (first + (pntVal - 1.0)*delta)/obs;
       }
    else if (!strcasecmp( specLabel, LAB_HZ ))
       {
        specVal  = first + (pntVal - 1.0)*delta;
       }
    else if (!strcasecmp( specLabel, LAB_WN ))
       {
        specVal  = (first + (pntVal - 1.0)*delta)/C_CM_SEC;
       }
    else if (!strcasecmp( specLabel, LAB_WNM ))
       {
        sw      = C_NM_SEC/sw;
        orig    = C_NM_SEC/orig - sw;

        if (size > 1)
           specVal = orig + sw*(pntVal - 1.0)/(float)(size - 1);
        else
           specVal = orig;
       }
    else if (!strcasecmp( specLabel, LAB_PCT ))
       {
        if (trueSize != 1)
           specVal = 100.0*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = 100.0;
       }
    else if (!strcasecmp( specLabel, LAB_SEC ))
       {
        if (size != 1)
           specVal = acqTime*(pntVal - 1)/(float)(size - 1);
        else
           specVal = acqTime;
       }
    else if (!strcasecmp( specLabel, LAB_MM ))
       {
        if (trueSize > 1)
           specVal = sw*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = sw;
       }
    else if (!strcasecmp( specLabel, LAB_CM ))
       {
        if (trueSize > 1)
           specVal = 0.1*sw*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = 0.1*sw;
       }
    else if (!strcasecmp( specLabel, LAB_MIC ))
       {
        if (trueSize > 1)
           specVal = 1000.0*sw*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = 1000.0*sw;
       }
    else if (!strcasecmp( specLabel, LAB_NM ))
       {
        if (trueSize != 1)
           specVal = 1000000.0*sw*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = 1000000.0*sw;
       }
    else if (!strcasecmp( specLabel, LAB_USER ))
       {
        if (trueSize > 1)
           specVal = s1 + uw*(oldPntVal - 1)/(float)(trueSize - 1);
        else
           specVal = uw;
       }
    else if (!*specLabel || !strcasecmp( specLabel, LAB_PTS ))
       {
        specVal = oldPntVal;
       }
    else
       {
        badSpecUnits = 1;
       }

    return( specVal );
}

/***/
/* specWidth2rPnt: convert a width in spectral units to real points
/***/

float specWidth2rPnt( fdata, dimCode, specLabel )

   float fdata[FDATASIZE];
   char  *specLabel;
   int   dimCode; 
{
    float acqTime, specVal, pntVal, sw, obs;
    char  label[NAMELEN+1];
    int   trueSize, size;

/***/
/* Initialize.
/* Extract header parameters for this axis.
/* Extract elements of spectral label.
/***/

    badSpecUnits = 0;
    specVal      = 0.0;
    pntVal       = 0.0;
    *label       = '\0';

    sw   = getParm( fdata, NDSW,   dimCode );
    obs  = getParm( fdata, NDOBS,  dimCode );

    if (sw == 0.0)  sw  = 1.0;
    if (obs == 0.0) obs = 1.0;

    (void) str2SpecVal( specLabel, &specVal, label );

/***/
/* Adjust size if this dimension has interleaved complex data.
/***/

    trueSize  = getParm( fdata, NDSIZE, dimCode );
    size      = trueSize;

    if (isInterleaved( fdata, dimCode )) size /= 2;

    acqTime = size/sw;

/***/
/* Handle special case for NM, wavelength or spatial distance.
/* Identify spectral units.
/* Perform conversion.
/***/

    if (!strcasecmp( label, LAB_NM ))
       {
        if (!isSpatialDim( fdata, dimCode )) (void) strcpy( label, LAB_WNM );
       }

    if (!strcasecmp( label, LAB_PPM ))
       {
        pntVal = (float)size*specVal*obs/sw; 
       }
    else if (!strcasecmp( label, LAB_HZ ))
       {
        pntVal = (float)size*specVal/sw; 
       }
    else if (!strcasecmp( label, LAB_USER ))
       {
        pntVal = (float)size*specVal/sw;
       }
    else if (!strcasecmp( label, LAB_WN ))
       {
        pntVal = (float)size*specVal*C_CM_SEC/sw;
       }
    else if (!strcasecmp( label, LAB_WNM ))
       {
        sw      = C_NM_SEC/sw;
        pntVal = (float)size*specVal/sw; 
       }
    else if (!strcasecmp( label, LAB_PCT ))
       {
        pntVal = (float)trueSize*specVal/100.0; 
       }
    else if (!strcasecmp( label, LAB_SEC ))
       {
        pntVal = acqTime == 0.0 ? 0.0 : (float)size*specVal/acqTime;
       }
    else if (!strcasecmp( label, LAB_MM ))
       {
        pntVal = (float)trueSize*specVal/sw;
       }
    else if (!strcasecmp( label, LAB_CM ))
       {
        pntVal = (float)trueSize*specVal/(0.1*sw);
       }
    else if (!strcasecmp( label, LAB_MIC ))
       {
        pntVal = (float)trueSize*specVal/(1000.0*sw);
       }
    else if (!strcasecmp( label, LAB_NM ))
       {
        pntVal = (float)trueSize*specVal/(1000000.0*sw);
       }
    else if (!*label || !strcasecmp( label, LAB_PTS ))
       {
        pntVal = specVal;
       }
    else
       {
        badSpecUnits = 1;
       }

    return( pntVal );
}

/***/
/* rPnt2specWidth: convert a real spectral width in pts to spectral units.
/***/

float rPnt2specWidth( fdata, dimCode, pntVal, label )

   float pntVal, fdata[FDATASIZE];
   char  *label;
   int   dimCode; 
{
    float acqTime, specVal, sw, obs, oldPntVal;
    int   trueSize, size;

    char  specLabel[NAMELEN+1];

/***/
/* Initialize.
/* Extract header parameters for this axis.
/***/

    (void) strcpy( specLabel, label );

    badSpecUnits = 0;
    specVal      = 0.0;

    sw   = getParm( fdata, NDSW,   dimCode );
    obs  = getParm( fdata, NDOBS,  dimCode );

    if (sw  == 0.0) sw  = 1.0;
    if (obs == 0.0) obs = 1.0;

/***/
/* Adjust size and point value if this dimension has interleaved complex data.
/***/

    trueSize  = getParm( fdata, NDSIZE, dimCode );
    size      = trueSize;
    oldPntVal = pntVal;

    if (isInterleaved( fdata, dimCode ))
       {
        if ((int)pntVal % 2)
           pntVal  = (pntVal + 1.0)/2.0; 
        else
           pntVal  = pntVal/2.0;

        size /= 2;
       }

    acqTime = size/sw;

/***/
/* Handle special case for NM, wavelength or spatial distance.
/* Identify spectral units.
/* Perform conversion.
/***/

    if (!strcasecmp( specLabel, LAB_NM ))
       {
        if (!isSpatialDim( fdata, dimCode )) (void)strcpy( specLabel, LAB_WNM );
       }

    if (!strcasecmp( specLabel, LAB_PPM ))
       {
        specVal = pntVal*sw/(obs*size);
       }
    else if (!strcasecmp( specLabel, LAB_HZ ))
       {
        specVal = pntVal*sw/(float)size;
       }
    else if (!strcasecmp( specLabel, LAB_USER ))
       {
        specVal = pntVal*sw/(float)size;
       }
    else if (!strcasecmp( specLabel, LAB_WN ))
       {
        specVal = (pntVal*sw/(float)size)/C_CM_SEC;
       }
    else if (!strcasecmp( specLabel, LAB_WNM ))
       {
        sw      = C_NM_SEC/sw;
        specVal = pntVal*sw/(float)size;
       }
    else if (!strcasecmp( specLabel, LAB_PCT ))
       {
        specVal = 100.0*oldPntVal/(float)trueSize; 
       }
    else if (!strcasecmp( specLabel, LAB_SEC ))
       {
        specVal = size == 0 ? 0 : pntVal*acqTime/(float)size;
       }
    else if (!strcasecmp( specLabel, LAB_MM ))
       {
        specVal = sw*oldPntVal/(float)trueSize;
       }
    else if (!strcasecmp( specLabel, LAB_CM ))
       {
        specVal = 0.1*sw*oldPntVal/(float)trueSize;
       }
    else if (!strcasecmp( specLabel, LAB_MIC ))
       {
        specVal = 1000.0*sw*oldPntVal/(float)trueSize;
       }
    else if (!strcasecmp( specLabel, LAB_NM ))
       {
        specVal = 1000000.0*sw*oldPntVal/(float)trueSize;
       }
    else if (!*specLabel || !strcasecmp( specLabel, LAB_PTS ))
       {
        specVal = oldPntVal;
       }
    else
       {
        badSpecUnits = 1;
       }

    return( specVal );
}

/***/
/* str2SpecVal: extract number and label from spectral value text.
/***/

int str2SpecVal( str, val, lab )

   char  *str, *lab;
   float *val;
{
    int status;

    status = str2LabVal( str, val, lab, validSpecUnits );

    return( status );
}

/***/
/* str2LabVal: extract number and label from text.
/***/

int str2LabVal( str, val, lab, labList )

   char  *str, *lab, **labList;
   float *val;
{

   char ctemp[NAMELEN+1];
   int  loc;

   *val = 0.0;
   *lab = '\0';

   loc  = hasUnitLabel( str, labList );

   if (loc)
      {
       (void) strcpy( ctemp, str );
       (void) strcpy( lab, str+loc );

       ctemp[loc] = '\0';

       (void) sscanf( ctemp, "%f", val );
      }
   else
      {
       (void) sscanf( str, "%f", val );
      }

   return( 1 );
}

/***/
/* hasSpecLabel: return label's loc if specStr contains a valid label.
/***/

int hasSpecLabel( specStr )

   char *specStr;
{
   int status;

   status = hasUnitLabel( specStr, validSpecUnits );

   return( status );
}

/***/
/* hasUnitLabel: return label's loc if specStr contains a valid label.
/***/

int hasUnitLabel( specStr, labList )

   char *specStr, **labList;
{
    int  i, labLen, strLen, idealLoc, actualLoc;
    char *sPtr, str[NAMELEN+1], lab[NAMELEN+1]; 

    (void) strcpy( str, specStr );
    (void) raiseStr( str );

    strLen = strlen( str );

    for( i = 0; labList[i]; i++ )
       {
        (void) strcpy( lab, labList[i] );
        (void) raiseStr( lab );

        labLen    = strlen( lab );
        idealLoc  = strLen - labLen;
        sPtr      = strstr( str, lab ); 
        actualLoc = sPtr ? sPtr - str : -1;
    
        if (idealLoc  < 1) continue;
        if (actualLoc < 1) continue;
 
        if (idealLoc == actualLoc) return( actualLoc );
       }
 
    return( 0 );
}

/***/
/* specList2Pnt: allocate integer point list corresponding to a list
/*               of spectral poisitions.
/***/

int *specList2Pnt( specList, specCount, fdata, dimCode )

   char  **specList;
   int   specCount, dimCode;
   float fdata[FDATASIZE];
{
    int i, *pntList;

    if (!(pntList = intAlloc( "spec", specCount ))) return( (int *)NULL );

    for( i = 0; i < specCount; i++ )
       {
        pntList[i] = spec2iPnt( fdata, dimCode, specList[i] );

        if (badSpecUnits)
           {
            (void) deAlloc( "spec", pntList, specCount );
            return( (int *)NULL );
           }
       }

    return( pntList );
}

/***/
/* specArgD: spectral-units-position equivalent of intArgD command line parsing.
/***/

int specArgD( argc, argv, flag, value, fdata, dimCode )

   int   argc, *value, dimCode;
   char  **argv, *flag;
   float fdata[FDATASIZE];
{
    int  loc;
    char *cPtr;

    if (loc = flagLoc( argc, argv, flag ))
       {
        if (cPtr = getNthArg( argc, argv, loc+1 ))
           {
            *value = 0.5 + spec2rPnt2( fdata, dimCode, cPtr, CLIP_TRUNC );

            if (badSpecUnits)
               {
                (void) fprintf( stderr, "SpecUnit Error: bad units.\n" );
                return( 0 );
               }

            return( 1 );
           }
       }

    return( 0 );
}

/***/
/* specFltArgD: spectral-units-position equivalent of fltArgD command parsing.
/***/

int specFltArgD( argc, argv, flag, value, fdata, dimCode )

   int   argc, dimCode;
   char  **argv, *flag;
   float *value, fdata[FDATASIZE];
{
    int  loc;
    char *cPtr;

    if (loc = flagLoc( argc, argv, flag ))
       {
        if (cPtr = getNthArg( argc, argv, loc+1 ))
           {
            *value = spec2rPnt2( fdata, dimCode, cPtr, CLIP_TRUNC );

            if (badSpecUnits)
               {
                (void) fprintf( stderr, "SpecUnit Error: bad units.\n" );
                return( 0 );
               }

            return( 1 );
           }
       }

    return( 0 );
}

/***/
/* wideArgD: spectral-units-width equivalent of intArgD command line parsing.
/***/

int wideArgD( argc, argv, flag, value, fdata, dimCode )

   int   argc, *value, dimCode;
   char  **argv, *flag;
   float fdata[FDATASIZE];
{
    int   loc;
    float wide;
    char  *cPtr;

    if (loc = flagLoc( argc, argv, flag ))
       {
        if (cPtr = getNthArg( argc, argv, loc+1 ))
           {
            wide   = specWidth2rPnt( fdata, dimCode, cPtr );
            *value = wide < 0.0 ? wide - 0.5 : wide + 0.5;

            if (badSpecUnits)
               {
                (void) fprintf( stderr, "SpecUnit Error: bad units.\n" );
                return( 0 );
               }

            return( 1 );
           }
       }

    return( 0 );
}

/***/
/* wideFltArgD: spectral-units-width equivalent of fltArgD command line parsing.
/***/

int wideFltArgD( argc, argv, flag, value, fdata, dimCode )

   int   argc, dimCode;
   char  **argv, *flag;
   float *value, fdata[FDATASIZE];
{
    int  loc;
    char *cPtr;

    if (loc = flagLoc( argc, argv, flag ))
       {
        if (cPtr = getNthArg( argc, argv, loc+1 ))
           {
            *value = specWidth2rPnt( fdata, dimCode, cPtr );

            if (badSpecUnits)
               {
                (void) fprintf( stderr, "SpecUnit Error: bad units.\n" );
                return( 0 );
               }

            return( 1 );
           }
       }

    return( 0 );
}

/***/
/* updateExtOrig: update the spectral origin of an extracted region.
/***/

int updateExtOrig( fdata, dimCode, ix1, ixn, fSize )

   float fdata[FDATASIZE];
   int   dimCode, ix1, ixn, fSize;
{
    int   midOld, midNew, mid, size;
    float orig, obs, car, sw;

/***/
/* Extract spectral parameters:
/***/

    sw     = getParm( fdata, NDSW,   dimCode );
    obs    = getParm( fdata, NDOBS,  dimCode );
    car    = getParm( fdata, NDCAR,  dimCode );
    size   = 1 + ixn - ix1;

    if (sw   == 0.0) sw   = 1.0;
    if (obs  == 0.0) obs  = 1.0;

    if (fSize == 0)  fSize = 1;
    if (size  == 0)  size  = 1;

    mid    = 1 + size/2;
    midOld = 1 + fSize/2;
    midNew = ix1 + size/2;

    orig = obs*car - sw*(fSize - mid)/fSize;

    car   += sw*(midOld - midNew)/(fSize*obs);
    orig  += sw*(midOld - midNew)/fSize;

    sw   *= (float)size/fSize;
    
    (void) setParm( fdata, NDCENTER, (float)mid, dimCode );
    (void) setParm( fdata, NDSW,     sw,         dimCode );
    (void) setParm( fdata, NDCAR,    car,        dimCode );
    (void) setParm( fdata, NDORIG,   orig,       dimCode );

    return( 0 );
}

/***/
/* updateOrig: update the spectral origin according to the carrier frequency.
/***/

int updateOrig( fdata, dimCode )

   float fdata[FDATASIZE];
   int   dimCode;
{
    int   ix1, ixn, mid, oldMid, size, oldSize;
    float orig, obs, car, sw;

/***/
/* Adjustment for previously extracted data:
/***/

    size = getParm( fdata, NDSIZE, dimCode );
    ix1  = getParm( fdata, NDX1, dimCode );
    ixn  = getParm( fdata, NDXN, dimCode );

    if (ix1 || ixn)
       {
        oldSize = 1 + ixn - ix1;
        oldMid  = 1 + oldSize/2;

        (void) setParm( fdata, NDSIZE, (float)oldSize, dimCode );

        car = iPnt2spec( fdata, dimCode, oldMid, "ppm" );

        (void) setParm( fdata, NDCAR, car, dimCode );
        (void) setParm( fdata, NDSIZE, (float)size, dimCode );
       }

/***/
/* Extract spectral parameters:
/***/

    mid = 1 + size/2;

    if (isInterleaved( fdata, dimCode )) size /= 2;

    sw  = getParm( fdata, NDSW,   dimCode );
    obs = getParm( fdata, NDOBS,  dimCode );
    car = getParm( fdata, NDCAR,  dimCode );

    if (sw   == 0.0) sw   = 1.0;
    if (obs  == 0.0) obs  = 1.0;
    if (size == 0)   size = 1;

/***/
/* Compute the origin:
/***/

    orig = obs*car - sw*(size - mid)/size;

    (void) setParm( fdata, NDCENTER, (float)mid, dimCode );
    (void) setParm( fdata, NDORIG,   orig,       dimCode );

    return( 0 );
}

/***/
/* getSpecUnits: return the units code associated with a spectral label.
/***/

int getSpecUnits( specLabel )

   char *specLabel;
{
    if (!strcasecmp( specLabel, LAB_PPM ))
       return( LAB_PPM_ID );
    else if (!strcasecmp( specLabel, LAB_HZ ))
       return( LAB_HZ_ID );
    else if (!strcasecmp( specLabel, LAB_PCT ))
       return( LAB_PCT_ID );
    else if (!strcasecmp( specLabel, LAB_WN ))
       return( LAB_WN_ID );
    else if (!strcasecmp( specLabel, LAB_MM ))
       return( LAB_MM_ID );
    else if (!strcasecmp( specLabel, LAB_CM ))
       return( LAB_CM_ID );
    else if (!strcasecmp( specLabel, LAB_MIC ))
       return( LAB_MIC_ID );
    else if (!strcasecmp( specLabel, LAB_NM ))
       return( LAB_NM_ID );
    else if (!strcasecmp( specLabel, LAB_WNM ))
       return( LAB_WNM_ID );
    else if (!strcasecmp( specLabel, LAB_SEC ))
       return( LAB_SEC_ID );
    else if (!strcasecmp( specLabel, LAB_USER ))
       return( LAB_USER_ID );

    return( LAB_PTS_ID );
}

/***/
/* getSpecLabel: return the units label associated with the units code.
/***/

int getSpecLabel( specID, specLabel )

   int  specID;
   char *specLabel;
{
    switch( specID )
       {
        case LAB_PPM_ID:
           (void) strcpy( specLabel, LAB_PPM );
           break;
        case LAB_PCT_ID:
           (void) strcpy( specLabel, LAB_PCT );
           break;
        case LAB_HZ_ID:
           (void) strcpy( specLabel, LAB_HZ );
           break;
        case LAB_WN_ID:
           (void) strcpy( specLabel, LAB_WN );
           break;
        case LAB_MM_ID:
           (void) strcpy( specLabel, LAB_MM );
           break;
        case LAB_CM_ID:
           (void) strcpy( specLabel, LAB_CM );
           break;
        case LAB_MIC_ID:
           (void) strcpy( specLabel, LAB_MIC );
           break;
        case LAB_NM_ID:
           (void) strcpy( specLabel, LAB_NM );
           break;
        case LAB_WNM_ID:
           (void) strcpy( specLabel, LAB_WNM );
           break;
        case LAB_SEC_ID:
           (void) strcpy( specLabel, LAB_SEC );
           break;
        case LAB_USER_ID:
           (void) strcpy( specLabel, LAB_USER );
           break;
        default:
           (void) strcpy( specLabel, LAB_PTS );
           break;
       }

    return( 0 );
}

int appendDefUnits( valStr, unitStr )

   char *valStr, *unitStr;
{
   int n;

   if (!valStr || *unitStr) return( 0 );

   n = strlen( valStr );
   if (!n) return( 0 );

   if (isalpha( valStr[n-1] )) return( 0 );

   (void) strcat( valStr, unitStr );

   return( 1 );
}

/***/
/* Bottom.
/***/
