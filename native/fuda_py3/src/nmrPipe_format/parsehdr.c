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
/* parseHdr: interpret command line arguments regarding file header.
/***/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define APLIST
#define ALTLIST
#define AQLIST
#define FTLIST
#define FTDOMLIST
#define ACQMETHLIST
#define MODELIST
#define UNITS_LIST
#define FDATA_LOCLIST

#include "cmndargs.h"
#include "fdatap.h"
#include "apod.h"
#include "prec.h"
#include "atof.h"
#include "specunit.h"
#include "gettime.h"
#include "namelist.h"

#define FPR (void)fprintf

int isInteger(); //AJB 14/09/22

int parseHdr( argc, argv, fdata )

   float fdata[FDATASIZE];
   char  **argv;
   int   argc;
{
    float rtemp;
    int   i, itemp, locX, locY, locZ, locA, ix, iy, iz, ia, status;

#define ANAMES "SP EM GM TM TRI or NONE"

/***/
/* Help Mode:
/***/

if (argc == 1 || flagLoc( argc, argv, "-help" ))
{
 FPR( stderr, "Adjust Header by Integer Location 1-512 or Keyword:\n" );
 FPR( stderr, " -fdata loc val ...  (See Below for Keywords).\n" );
 FPR( stderr, "General Header Arguments:\n" );
 FPR( stderr, " -u1 -u2...-u5   User Parameters.\n" );
 FPR( stderr, " -ndim           Number of Dimensions.\n" );
 FPR( stderr, " -tau            Tau Value.\n" );
 FPR( stderr, " -day            Day.\n" );
 FPR( stderr, " -month          Month [1-12].\n" );
 FPR( stderr, " -year           Year.\n" );
 FPR( stderr, " -hours          Hours [0-23].\n" );
 FPR( stderr, " -mins           Minutes [0-59].\n" );
 FPR( stderr, " -secs           Seconds [0-59].\n" );
 FPR( stderr, " -temperature    Temperature (Also: -temp).\n" );
 FPR( stderr, " -pressure       Pressure (Also: -pres).\n" );
 FPR( stderr, " -nusDim         NUS Dimensions (For Unexpanded Data).\n" );
 FPR( stderr, " -aq2D           2D Acquisition Code or Keyword:\n" );
 FPR( stderr, "                    0 Magnitude.\n" );
 FPR( stderr, "                    1 Real.\n" );
 FPR( stderr, "                    1 TPPI.\n" );
 FPR( stderr, "                    2 Complex.\n" );
 FPR( stderr, "                    2 States.\n" );
 FPR( stderr, "                    3 Image.\n" );
 FPR( stderr, " -fileCount      Number of Files in Data Series.\n" );
 FPR( stderr, " -pipeFlag       Mult-File or Stream Format:\n" );
 FPR( stderr, "                    0 = Single File or 2D Series Data.\n" );
 FPR( stderr, "                    1 = Pipeline Stream or Single File (Also 2 3 4).\n" );
 FPR( stderr, " -cubeFlag       4D Cube Format:\n" );
 FPR( stderr, "                    1 = 3D Series (4D Only).\n" );
 FPR( stderr, "                    0 = 2D Series or 4D Stream.\n" );
 FPR( stderr, " -tpFlag         2D Transpose State of Planes:\n" );
 FPR( stderr, "                    0 = Not Transposed.\n" );
 FPR( stderr, "                    1 = Transposed.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Header Labels:\n" );
 FPR( stderr, " -srcName          Source File Name.\n" );
 FPR( stderr, " -operName         Operator Name.\n" );
 FPR( stderr, " -userName         User Name.\n" );
 FPR( stderr, " -comment          Comment.\n" );
 FPR( stderr, " -title            Title.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Axis-Specific Header Flags -x.. -y.. etc;\n" );
 FPR( stderr, "Default Values are zero:\n" );
 FPR( stderr, " -xN      Actual Size in File           Pts Real + Imag.\n" );
 FPR( stderr, " -xT      Current Time Domain Size      Pts Real or Imag.\n" );
 FPR( stderr, " -xSW     Full Spectral Width           Hz.\n" );
 FPR( stderr, " -xOBS    Observe Frequency at 0.0ppm   MHz (Used for ppm).\n" );
 FPR( stderr, " -xOBSMID Observe Frequency at Center   MHz (Optional).\n" );
 FPR( stderr, " -xCAR    Carrier Position              ppm.\n" );
 FPR( stderr, " -xLAB    Axis Label                    8 Chars, No Spaces.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Note: OBSMID is used for conversion schemes which\n" );
 FPR( stderr, "      adjust the observe frequency to correspond to the\n" );
 FPR( stderr, "      value at 0.0ppm rather than at the center of the\n" );
 FPR( stderr, "      spectrum, for rigorous ppm calculation. In most\n" );
 FPR( stderr, "      cases, the original and adjusted values are more\n" );
 FPR( stderr, "      or less the same. The OBS value is used in ppm\n" );
 FPR( stderr, "      calculations.  The OBSMID value is used only as a\n" );
 FPR( stderr, "      record of the unadjusted OBS value.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, " -xMODE  Acquisition Mode Code or Keyword:\n" );
 FPR( stderr, "          0=Complex,States,Complex-N,States-N.\n" );
 FPR( stderr, "          0=States-TPPI,States-TPPI-N.\n" );
 FPR( stderr, "          1=Real,TPPI.\n" );
 FPR( stderr, "          2=Sequential,Bruker.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, " -xALT   Sign Adjustment Code or Keyword:\n" );
 FPR( stderr, "          0=Complex,States    No Alternation Needed.\n" );
 FPR( stderr, "          0=Real,TPPI         No Alternation Needed.\n" );
 FPR( stderr, "          1=Sequential,Bruker Real Alternation.\n" );
 FPR( stderr, "          2=States-TPPI       Complex Alternation.\n" );
 FPR( stderr, "         With Negation of Imaginaries:\n" );
 FPR( stderr, "          16=Complex-N        No Alternation Needed.\n" );
 FPR( stderr, "          16=States-N         No Alternation Needed.\n" );
 FPR( stderr, "          18=States-TPPI-N    Complex Alternation.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Note: Default Sign Adjustment Mode is Automatically\n" );
 FPR( stderr, "      Set by the Corresponding Acquisition Mode Keyword.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Axis-Specific Header Flags -xQ -yQ etc;\n" );
 FPR( stderr, "Not required for conversion; used to adjust header:\n" );
 FPR( stderr, " -xTD    Original Time Domain Size.  Pts.\n" );
 FPR( stderr, " -xFD    Final Max Freq Domain Size. Pts.\n" );
 FPR( stderr, " -xZFP   Zero Fill Parameter. Pts (Negative).\n" );
 FPR( stderr, " -xMID   Carrier Virtual Loc.        Pts.\n" );
 FPR( stderr, " -xORIG  Last Point Coord.           Hz.\n" );
 FPR( stderr, " -xDIM   Dimension Order.\n" );
 FPR( stderr, "          2=Original X-Axis\n" ); 
 FPR( stderr, "          1=Original Y-Axis\n" ); 
 FPR( stderr, "          3=Original Z-Axis\n" ); 
 FPR( stderr, "          4=Original A-Axis\n" ); 
 FPR( stderr, " -xFT    Fourier Flag         Code or Keyword:\n" );
 FPR( stderr, "          0=Time              Time-Domain Data.\n" );
 FPR( stderr, "          1=Freq              Frequency-Domain Data.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Special Acquisition-Axis Flags (no longer required):\n" );
 FPR( stderr, " -aqN     Actual Size.         Pts Real or Imag.\n" );
 FPR( stderr, " -aqMODE  Acquisition Mode.    0=Complex.\n" );
 FPR( stderr, "                               1=Real.\n" );
 FPR( stderr, "                               2=Sequential.\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Axis-Specific Processing Values -xQ -yQ etc:\n" );
 FPR( stderr, " -xAPOD  Apodize Function      %s.\n", ANAMES );
 FPR( stderr, " -xQ1    Apodize Parameter 1.  Depends on Function.\n" );
 FPR( stderr, " -xQ2    Apodize Parameter 2.  Depends on Function.\n" );
 FPR( stderr, " -xQ3    Apodize Parameter 3.  Depends on Function.\n" );
 FPR( stderr, " -xC1    First Point Scale.    Default is 1.0.\n" );
 FPR( stderr, " -xP0    Zero Order Phase.     Degrees.\n" );
 FPR( stderr, " -xP1    First Order Phase.    Degrees.\n" );
 FPR( stderr, "Values for Interpretation of Axis Units:\n" );
 FPR( stderr, " -xUNITS Axis Units            Pts Hz ppm %%\n" );
 FPR( stderr, "                               WN cm mm mic nm WL-nm\n" );
 FPR( stderr, "                               in c1 c2 pix Hi\n" );
 FPR( stderr, " -xMETH  Acquistion Method     0=Fourier (FT)\n" );
 FPR( stderr, "                               1=Direct.\n" );
 FPR( stderr, " -xFTDOM Domain of FT Data     0=Freq (Spectral)\n" );
 FPR( stderr, "                               1=Dist (Spatial)\n" );
 FPR( stderr, "\n" );
 FPR( stderr, "Keywords for -fdata Locations (1 - 512):\n" );

 for( i = 0; fdataLocList[i].name; i++ )
    {
     itemp = 1 + (int)fdataLocList[i].val;
     if (itemp < 0 || itemp > FDATASIZE) continue;

     FPR( stderr, " %-14s %3d\n", fdataLocList[i].name, itemp );
    }

 FPR( stderr, "\n" );

 return( 1 );
}

/***/
/* Special aquisition dimension parameters:
/***/

    if (0 > fltStrArgD( argc, argv, "-aqMODE", fdata + FDQUADFLAG, modeList ))
       {
        return( 1 );
       }

    if (0 > fltStrArgD( argc, argv, "-aqMODE", fdata + FDF2AQSIGN, altList ))
       {
        return( 1 );
       }

    if (0 > fltStrArgD( argc, argv, "-aqALT",  fdata + FDF2AQSIGN, altList ))
       {
        return( 1 );
       }

    (void) fltArgD( argc, argv, "-aqN", fdata + FDREALSIZE );

/***/
/* Number of points per dimension in file:
/***/

    (void) fltArgD( argc, argv, "-xN", fdata + FDSIZE );
    (void) fltArgD( argc, argv, "-yN", fdata + FDSPECNUM );
    (void) fltArgD( argc, argv, "-zN", fdata + FDF3SIZE );
    (void) fltArgD( argc, argv, "-aN", fdata + FDF4SIZE );

/***/
/* Current/Original time-domain length, final/max freq-domain length, zf:
/***/

    (void) fltArgD( argc, argv, "-xT",  fdata + FDF2APOD );
    (void) fltArgD( argc, argv, "-yT",  fdata + FDF1APOD );
    (void) fltArgD( argc, argv, "-zT",  fdata + FDF3APOD );
    (void) fltArgD( argc, argv, "-aT",  fdata + FDF4APOD );

    (void) fltArgD( argc, argv, "-xTD", fdata + FDF2TDSIZE );
    (void) fltArgD( argc, argv, "-yTD", fdata + FDF1TDSIZE );
    (void) fltArgD( argc, argv, "-zTD", fdata + FDF3TDSIZE );
    (void) fltArgD( argc, argv, "-aTD", fdata + FDF4TDSIZE );

    (void) fltArgD( argc, argv, "-xFD", fdata + FDF2FTSIZE );
    (void) fltArgD( argc, argv, "-yFD", fdata + FDF1FTSIZE );
    (void) fltArgD( argc, argv, "-zFD", fdata + FDF3FTSIZE );
    (void) fltArgD( argc, argv, "-aFD", fdata + FDF4FTSIZE );

    (void) fltArgD( argc, argv, "-xZFP", fdata + FDF2ZF );
    (void) fltArgD( argc, argv, "-yZFP", fdata + FDF1ZF );
    (void) fltArgD( argc, argv, "-zZFP", fdata + FDF3ZF );
    (void) fltArgD( argc, argv, "-aZFP", fdata + FDF4ZF );

/***/
/* Spectral Width:
/***/

    (void) fltArgD( argc, argv, "-xSW", fdata + FDF2SW );
    (void) fltArgD( argc, argv, "-ySW", fdata + FDF1SW );
    (void) fltArgD( argc, argv, "-zSW", fdata + FDF3SW );
    (void) fltArgD( argc, argv, "-aSW", fdata + FDF4SW );

/***/
/* Observe frequency:
/***/

    (void) fltArgD( argc, argv, "-xOBS", fdata + FDF2OBS );
    (void) fltArgD( argc, argv, "-yOBS", fdata + FDF1OBS );
    (void) fltArgD( argc, argv, "-zOBS", fdata + FDF3OBS );
    (void) fltArgD( argc, argv, "-aOBS", fdata + FDF4OBS );

    (void) fltArgD( argc, argv, "-xOBSMID", fdata + FDF2OBSMID );
    (void) fltArgD( argc, argv, "-yOBSMID", fdata + FDF1OBSMID );
    (void) fltArgD( argc, argv, "-zOBSMID", fdata + FDF3OBSMID );
    (void) fltArgD( argc, argv, "-aOBSMID", fdata + FDF4OBSMID );

/***/
/* Point location of carrier:
/***/

    (void) fltArgD( argc, argv, "-xMID", fdata + FDF2CENTER );
    (void) fltArgD( argc, argv, "-yMID", fdata + FDF1CENTER );
    (void) fltArgD( argc, argv, "-zMID", fdata + FDF3CENTER );
    (void) fltArgD( argc, argv, "-aMID", fdata + FDF4CENTER );

/***/
/* Desired phase corrections:
/***/

    (void) fltArgD( argc, argv, "-xP0", fdata + FDF2P0 );
    (void) fltArgD( argc, argv, "-yP0", fdata + FDF1P0 );
    (void) fltArgD( argc, argv, "-zP0", fdata + FDF3P0 );
    (void) fltArgD( argc, argv, "-aP0", fdata + FDF4P0 );

    (void) fltArgD( argc, argv, "-xP1", fdata + FDF2P1 );
    (void) fltArgD( argc, argv, "-yP1", fdata + FDF1P1 );
    (void) fltArgD( argc, argv, "-zP1", fdata + FDF3P1 );
    (void) fltArgD( argc, argv, "-aP1", fdata + FDF4P1 );

/***/
/* Apodization Parameters:
/* Note special case of first point scale, which is stored as C1 - 1.0.
/***/

    if (0 > fltStrArgD(argc, argv, "-xAPOD", fdata+FDF2APODCODE, apList))
       {
        return( 1 );
       }

    if (0 > fltStrArgD(argc, argv, "-yAPOD", fdata+FDF1APODCODE, apList))
       {
        return( 1 );
       }

    if (0 > fltStrArgD(argc, argv, "-zAPOD", fdata+FDF3APODCODE, apList))
       {
        return( 1 );
       }

    if (0 > fltStrArgD(argc, argv, "-aAPOD", fdata+FDF4APODCODE, apList))
       {
        return( 1 );
       }

    (void) fltArgD( argc, argv, "-xQ1", fdata + FDF2APODQ1 );
    (void) fltArgD( argc, argv, "-yQ1", fdata + FDF1APODQ1 );
    (void) fltArgD( argc, argv, "-zQ1", fdata + FDF3APODQ1 );
    (void) fltArgD( argc, argv, "-aQ1", fdata + FDF4APODQ1 );

    (void) fltArgD( argc, argv, "-xQ2", fdata + FDF2APODQ2 );
    (void) fltArgD( argc, argv, "-yQ2", fdata + FDF1APODQ2 );
    (void) fltArgD( argc, argv, "-zQ2", fdata + FDF3APODQ2 );
    (void) fltArgD( argc, argv, "-aQ2", fdata + FDF4APODQ2 );

    (void) fltArgD( argc, argv, "-xQ3", fdata + FDF2APODQ3 );
    (void) fltArgD( argc, argv, "-yQ3", fdata + FDF1APODQ3 );
    (void) fltArgD( argc, argv, "-zQ3", fdata + FDF3APODQ3 );
    (void) fltArgD( argc, argv, "-aQ3", fdata + FDF4APODQ3 );

    if (flagLoc( argc, argv, "-xC1" ))
       {
        fdata[FDF2C1] = fltArg( argc, argv, "-xC1" ) - 1.0;
       }

    if (flagLoc( argc, argv, "-yC1" ))
       {
        fdata[FDF1C1] = fltArg( argc, argv, "-yC1" ) - 1.0;
       }

    if (flagLoc( argc, argv, "-zC1" ))
       {
        fdata[FDF3C1] = fltArg( argc, argv, "-zC1" ) - 1.0;
       }

    if (flagLoc( argc, argv, "-aC1" ))
       {
        fdata[FDF4C1] = fltArg( argc, argv, "-aC1" ) - 1.0;
       }

/***/
/* FT Domain Flags Spatial/Distance or Freq/Wavelength:
/***/

    ix = ((int)fdata[FDDOMINFO]) & 1;
    iy = ((int)fdata[FDDOMINFO]) & 2;
    iz = ((int)fdata[FDDOMINFO]) & 4;
    ia = ((int)fdata[FDDOMINFO]) & 8;

    locX = intStrArgD( argc, argv, "-xFTDOM", &ix, ftDomList );
    if (locX < 0) return( 1 );

    locY = intStrArgD( argc, argv, "-yFTDOM", &iy, ftDomList );
    if (locY < 0) return( 1 );

    locZ = intStrArgD( argc, argv, "-zFTDOM", &iz, ftDomList );
    if (locZ < 0) return( 1 );

    locA = intStrArgD( argc, argv, "-aFTDOM", &ia, ftDomList );
    if (locA < 0) return( 1 );

    if (locX || locY || locZ || locA)
       {
        fdata[FDDOMINFO] = ix + 2*iy + 4*iz + 8*ia;
       }

/***/
/* Acq Method Flags, Fourier or Direct:
/***/

    ix = ((int)fdata[FDMETHINFO]) & 1;
    iy = ((int)fdata[FDMETHINFO]) & 2;
    iz = ((int)fdata[FDMETHINFO]) & 4;
    ia = ((int)fdata[FDMETHINFO]) & 8;

    locX = intStrArgD( argc, argv, "-xMETH", &ix, acqMethList );
    if (locX < 0) return( 1 );

    locY = intStrArgD( argc, argv, "-yMETH", &iy, acqMethList );
    if (locY < 0) return( 1 );

    locZ = intStrArgD( argc, argv, "-zMETH", &iz, acqMethList );
    if (locZ < 0) return( 1 );

    locA = intStrArgD( argc, argv, "-aMETH", &ia, acqMethList );
    if (locA < 0) return( 1 );

    if (locX || locY || locZ || locA)
       { 
        fdata[FDMETHINFO] = ix + 2*iy + 4*iz + 8*ia;
       }

/***/
/* Fourier Transform Flags (Time/Frequency Domain):
/***/

    if (0 > fltStrArgD( argc, argv, "-xFT", fdata + FDF2FTFLAG, ftList ))
       {
        return( 1 );
       }

    if (0 > fltStrArgD( argc, argv, "-yFT", fdata + FDF1FTFLAG, ftList ))
       {
        return( 1 );
       }

    if (0 > fltStrArgD( argc, argv, "-zFT", fdata + FDF3FTFLAG, ftList ))
       {
        return( 1 );
       }

    if (0 > fltStrArgD( argc, argv, "-aFT", fdata + FDF4FTFLAG, ftList ))
       {
        return( 1 );
       }

/***/
/* Aprox. origin in terms of Spectral Width and Observe Frequency.
/* Then, explicit spectral axis origin:
/***/

    (void) fltArgD( argc, argv, "-xCAR", fdata + FDF2CAR );
    (void) fltArgD( argc, argv, "-yCAR", fdata + FDF1CAR );
    (void) fltArgD( argc, argv, "-zCAR", fdata + FDF3CAR );
    (void) fltArgD( argc, argv, "-aCAR", fdata + FDF4CAR );

    if (flagLoc( argc, argv, "-xCAR" ))
       {
        rtemp           = fltArg( argc, argv, "-xCAR" );
        fdata[FDF2ORIG] = -0.5*fdata[FDF2SW] + fdata[FDF2OBS]*rtemp;
       }

    if (flagLoc( argc, argv, "-yCAR" ))
       {
        rtemp           = fltArg( argc, argv, "-yCAR" );
        fdata[FDF1ORIG] = -0.5*fdata[FDF1SW] + fdata[FDF1OBS]*rtemp;
       }

    if (flagLoc( argc, argv, "-zCAR" ))
       {
        rtemp           = fltArg( argc, argv, "-zCAR" );
        fdata[FDF3ORIG] = -0.5*fdata[FDF3SW] + fdata[FDF3OBS]*rtemp;
       }

    if (flagLoc( argc, argv, "-aCAR" ))
       {
        rtemp           = fltArg( argc, argv, "-aCAR" );
        fdata[FDF4ORIG] = -0.5*fdata[FDF4SW] + fdata[FDF4OBS]*rtemp;
       }

    (void) fltArgD( argc, argv, "-xORIG", fdata + FDF2ORIG );
    (void) fltArgD( argc, argv, "-yORIG", fdata + FDF1ORIG );
    (void) fltArgD( argc, argv, "-zORIG", fdata + FDF3ORIG );
    (void) fltArgD( argc, argv, "-aORIG", fdata + FDF4ORIG );

/***/
/* Quadrature aquisition mode and corresponding sign-alternation:
/***/

    status = fltStrArgD( argc, argv, "-xMODE", fdata + FDF2QUADFLAG, modeList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -xMODE should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }
 
    status = fltStrArgD( argc, argv, "-xMODE", fdata + FDF2AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }

    status = fltStrArgD( argc, argv, "-xALT",  fdata + FDF2AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -xALT  should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

    status = fltStrArgD( argc, argv, "-yMODE", fdata + FDF1QUADFLAG, modeList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -yMODE should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

    status = fltStrArgD( argc, argv, "-yMODE", fdata + FDF1AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }

    status = fltStrArgD( argc, argv, "-yALT",  fdata + FDF1AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -yALT  should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

    status = fltStrArgD( argc, argv, "-zMODE", fdata + FDF3QUADFLAG, modeList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -zMODE should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

    status = fltStrArgD( argc, argv, "-zMODE", fdata + FDF3AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }

    status = fltStrArgD( argc, argv, "-zALT",  fdata + FDF3AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -zALT  should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }


    status = fltStrArgD( argc, argv, "-aMODE", fdata + FDF4QUADFLAG, modeList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -aMODE should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

    status = fltStrArgD( argc, argv, "-aMODE", fdata + FDF4AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }

    status = fltStrArgD( argc, argv, "-aALT",  fdata + FDF4AQSIGN,   altList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -aALT  should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'Complex', not as a number.\n" );
       }

/***/
/* Special Axis Units:
/***/

    status = fltStrArgD( argc, argv, "-xUNITS",  fdata + FDF2UNITS, unitLabelList );
    if (status == KEYWORD_BAD) return( 1 );

    status = fltStrArgD( argc, argv, "-yUNITS",  fdata + FDF1UNITS, unitLabelList );
    if (status == KEYWORD_BAD) return( 1 );

    status = fltStrArgD( argc, argv, "-zUNITS",  fdata + FDF3UNITS, unitLabelList );
    if (status == KEYWORD_BAD) return( 1 );

    status = fltStrArgD( argc, argv, "-aUNITS",  fdata + FDF4UNITS, unitLabelList );
    if (status == KEYWORD_BAD) return( 1 );

/***/
/* User parameters:
/***/

    (void) fltArgD( argc, argv, "-u1", fdata + FDUSER1 );
    (void) fltArgD( argc, argv, "-u2", fdata + FDUSER2 );
    (void) fltArgD( argc, argv, "-u3", fdata + FDUSER3 );
    (void) fltArgD( argc, argv, "-u4", fdata + FDUSER4 );
    (void) fltArgD( argc, argv, "-u5", fdata + FDUSER5 );

/***/
/* Miscellaneous:
/***/

    status = fltStrArgD( argc, argv, "-aq2D",   fdata + FD2DPHASE, aq2DList );

    if (status == KEYWORD_BAD)
       {
        return( 1 );
       }
    else if (status == KEYWORD_NUM)
       {
        FPR( stderr, "Warning from %s:\n", getProgName() );
        FPR( stderr, "  Parameter -aq2D should now be specified as a\n" );
        FPR( stderr, "  keyword such as 'States', not as a number.\n" );
       }

    (void) fltArgD( argc, argv, "-fileCount", fdata + FDFILECOUNT );
    (void) fltArgD( argc, argv, "-pipeFlag",  fdata + FDPIPEFLAG );
    (void) fltArgD( argc, argv, "-cubeFlag",  fdata + FDCUBEFLAG );

    (void) fltArgD( argc, argv, "-tpFlag",      fdata + FDTRANSPOSED );
    (void) fltArgD( argc, argv, "-tau",         fdata + FDTAU );
    (void) fltArgD( argc, argv, "-day",         fdata + FDDAY );
    (void) fltArgD( argc, argv, "-month",       fdata + FDMONTH );
    (void) fltArgD( argc, argv, "-year",        fdata + FDYEAR );
    (void) fltArgD( argc, argv, "-temp",        fdata + FDTEMPERATURE );
    (void) fltArgD( argc, argv, "-temperature", fdata + FDTEMPERATURE );
    (void) fltArgD( argc, argv, "-pres",        fdata + FDPRESSURE );
    (void) fltArgD( argc, argv, "-pressure",    fdata + FDPRESSURE );
    (void) fltArgD( argc, argv, "-hours",       fdata + FDHOURS );
    (void) fltArgD( argc, argv, "-mins",        fdata + FDMINS );
    (void) fltArgD( argc, argv, "-secs",        fdata + FDSECS );
    (void) fltArgD( argc, argv, "-scans",       fdata + FDSCANS );
    (void) fltArgD( argc, argv, "-dmxval",      fdata + FDDMXVAL );
    (void) fltArgD( argc, argv, "-dmxflag",     fdata + FDDMXFLAG );

    (void) fltArgD( argc, argv, "-ndim",     fdata + FDDIMCOUNT );
    (void) fltArgD( argc, argv, "-nusDim",   fdata + FDNUSDIM );

    (void) fltArgD( argc, argv, "-xDIM",   fdata + FDDIMORDER1 );
    (void) fltArgD( argc, argv, "-yDIM",   fdata + FDDIMORDER2 );
    (void) fltArgD( argc, argv, "-zDIM",   fdata + FDDIMORDER3 );
    (void) fltArgD( argc, argv, "-aDIM",   fdata + FDDIMORDER4 );

/***/
/* Text parameters:
/***/

    (void) fdTxtD( argc, argv, "-xLAB", SIZE_F2LABEL, fdata + FDF2LABEL );
    (void) fdTxtD( argc, argv, "-yLAB", SIZE_F1LABEL, fdata + FDF1LABEL );
    (void) fdTxtD( argc, argv, "-zLAB", SIZE_F3LABEL, fdata + FDF3LABEL );
    (void) fdTxtD( argc, argv, "-aLAB", SIZE_F4LABEL, fdata + FDF4LABEL );

    (void) fdTxtD( argc, argv, "-srcName",  SIZE_SRCNAME,  fdata + FDSRCNAME );
    (void) fdTxtD( argc, argv, "-userName", SIZE_USERNAME, fdata + FDUSERNAME );
    (void) fdTxtD( argc, argv, "-operName", SIZE_OPERNAME, fdata + FDOPERNAME );
    (void) fdTxtD( argc, argv, "-comment",  SIZE_COMMENT,  fdata + FDCOMMENT );
    (void) fdTxtD( argc, argv, "-title",    SIZE_TITLE,    fdata + FDTITLE );

    return( 0 );
}

/***/
/* parseHdr2: interprets command line parameters for direct header values.
/***/

int parseHdr2( argc, argv, fdata )

   int   argc;
   char  **argv;
   float fdata[FDATASIZE];
{
    int   itemp, loc1, loc2, status;
    float rtemp;
    char  *sPtr;

    rtemp = 0.0;

    if ((loc1 = flagLoc( argc, argv, "-fdata" )))
       {
        loc2 = nextFlag( argc, argv, loc1 );

        while( loc1 < loc2 - 1 )
           {
            sPtr  = getNthArg( argc, argv, loc1+1 );
            itemp = -1;

            if (isInteger( sPtr ))
               {
                itemp = (ATOI( sPtr )) - 1;
               }
            else
               {
                status = getValByName( "-fdata", getNthArg( argc, argv, loc1+1 ), fdataLocList, &rtemp );

                if (status == KEYWORD_BAD) 
                   FPR( stderr, "Error in -fdata Arguments: Bad FDATA Location %s\n", sPtr ); 
                else
                   itemp = (int)rtemp;
               }

            if (itemp < 0 || itemp > FDATASIZE - 1) 
               {
                FPR( stderr, "Error in -fdata Arguments: Bad FDATA Location or Keyword %s\n", sPtr ); 
               }
            else
               {
                rtemp = ATOF( getNthArg( argc, argv, loc1+2 ));
                (void) setParm( fdata, itemp, rtemp, 0 );
               }

            loc1 += 2;
           }
       }

    return( 0 );
}

/***/
/* setHdrTime: sets the header values with the current date and time.
/***/

int setHdrTime( fdata )

   float fdata[FDATASIZE];
{
    int month, day, year, hour, min, sec;

    (void) getTime( &month, &day, &year, &hour, &min, &sec );

    (void) setParm( fdata, FDMONTH, (float)month, NULL_DIM );
    (void) setParm( fdata, FDDAY,   (float)day,   NULL_DIM );
    (void) setParm( fdata, FDYEAR,  (float)year,  NULL_DIM );
    (void) setParm( fdata, FDHOURS, (float)hour,  NULL_DIM );
    (void) setParm( fdata, FDMINS,  (float)min,   NULL_DIM );
    (void) setParm( fdata, FDSECS,  (float)sec,   NULL_DIM );

    return( 0 );
}

/***/
/* Bottom.
/***/
