/**************************************************/
/* nmrPipe class         */
/* A.Baldwin     */
/* 19th March 2018*/
/**************************************************/

#ifndef PIPECLASS_HPP
#define PIPECLASS_HPP


//#include <iostream>
//#include <sstream>
//#include <bitset>
//#include <cstdio>
//#include <string>
//#include <vector>

#include "slice.hpp"



//class slice4D;
//class slice3D;
//class slice2D;
//class slice1D;

class pipe
{
  public:

  int dim;
  int *size,*dimord;
  float *sw,*car,*obs,*orig;
  int *xn,*x1,*center;
  float **xvals; //pointer to frequency ranges
  std::string iName;
  std::string oName= "raw\\tmp";
  vector<peakEntry> peakList;
  int peaks;
  FILE*A;
  FILE*B;
  double noise=1E6;
  vector<string> labels;

  float prehead[512]; //space for the nmrpipe preheader

  int LetterTest(char test);
  string GetLabel(int posy);
  void readHeader();
  int GetIndex(float val,int col);
  void SetPos(int &i,int &j,int &k,int &l);
  void SetPosInt(int *loc);
  void SetPos3D(int &i,int &j,int &k);
  void SetPos2D(int &i,int &j);
  void SetPos1D(int &i);
  double GetData(int *loc);
  void LocalMax(int p);
  void WritePeak();
  void WritePipe4D(float *arr,int si,int sj,int sk, int sl);

#ifdef DOUBLE3D
  void WritePipe3D(double *arr,int si,int sj,int sk);
  #else
  void WritePipe3D(float *arr,int si,int sj,int sk);
   #endif
  void WritePipe2D(double *arr,int si,int sj,const std::string &tag = "decon");
  void WritePipe2D(float *arr,int si,int sj,const std::string &tag = "decon");

#ifdef DOUBLE1D
  void WritePipe1D(double *arr,int si);
#else
  void WritePipe1D(float *arr,int si);
  #endif
  void optPeak4D();
  float GetVal3D(int &i,int &j, int &k);
  float GetVal2D(int &i,int &j);
  float GetVal1D(int &i);
  void extract4D();
  void diag4D();
  void project4D();

  ~pipe()
  {

    delete [] size;
    delete [] dimord;
    delete [] sw;
    delete [] car;
    delete [] obs;
    delete [] orig;
    delete [] xn;
    delete [] x1;
    delete [] center;

    for(int i=0;i<4;++i)
      delete [] xvals[i];
    delete [] xvals;

    
    //int *size,*dimord;
  
    //float *sw,*car,*obs,*orig;
    //int *xn,*x1,*center;

    //float **xvals; //pointer to frequency ranges

    
  }
  
  /*
    //from NMRglue, pipe.py, dictionary indicating
    //where in the header the various parameters live
fdata_dic = {
'FDF4CENTER': '82', 'FDF2P0': '109', 'FDF2P1': '110', 'FDF1P1': '246',
'FDF2X1': '257', 'FDF1P0': '245', 'FDF3AQSIGN': '476', 'FDDISPMAX': '251',
'FDF4FTFLAG': '31', 'FDF3X1': '261', 'FDRANK': '180', 'FDF2C1': '418',
'FDF2QUADFLAG': '56', 'FDSLICECOUNT': '443', 'FDFILECOUNT': '442',
'FDMIN': '248', 'FDF3OBS': '10', 'FDF4APODQ2': '407', 'FDF4APODQ1': '406',
'FDF3FTSIZE': '200', 'FDF1LB': '243', 'FDF4C1': '409', 'FDF4QUADFLAG': '54',
'FDF1SW': '229', 'FDTRANSPOSED': '221', 'FDSECS': '285', 'FDF1APOD': '428',
'FDF2APODCODE': '413', 'FDPIPECOUNT': '75', 'FDOPERNAME': '464',
'FDF3LABEL': '20', 'FDPEAKBLOCK': '362', 'FDREALSIZE': '97', 'FDF4SIZE': '32',
'FDF4SW': '29', 'FDF4ORIG': '30', 'FDF3XN': '262', 'FDF1OBS': '218',
'FDDISPMIN': '252', 'FDF2XN': '258', 'FDF3P1': '61', 'FDF3P0': '60',
'FDF1ORIG': '249', 'FDF2FTFLAG': '220', 'FDF1TDSIZE': '387',
'FDLASTPLANE': '78',
'FDF1ZF': '437', 'FDF4FTSIZE': '201', 'FDF3C1': '404', 'FDFLTFORMAT': '1',
'FDF4CAR': '69', 'FDF1FTFLAG': '222', 'FDF2OFFPPM': '480', 'FDF1LABEL': '18',
'FDSIZE': '99', 'FDYEAR': '296', 'FDF1C1': '423', 'FDUSER3': '72',
'FDF1FTSIZE': '98', 'FDMINS': '284', 'FDSCALEFLAG': '250', 'FDF3TDSIZE': '388',
'FDTITLE': '297', 'FDPARTITION': '65', 'FDF3FTFLAG': '13', 'FDF2APODQ1': '415',
'FD2DVIRGIN': '399', 'FDF2APODQ3': '417', 'FDF2APODQ2': '416',
'FD2DPHASE': '256', 'FDMAX': '247', 'FDF3SW': '11', 'FDF4TDSIZE': '389',
'FDPIPEFLAG': '57', 'FDDAY': '295', 'FDF2UNITS': '152', 'FDF4APODQ3': '408',
'FDFIRSTPLANE': '77', 'FDF3SIZE': '15', 'FDF3ZF': '438', 'FDDIMORDER': '24',
'FDF3ORIG': '12', 'FD1DBLOCK': '365', 'FDF1AQSIGN': '475', 'FDF2OBS': '119',
'FDF1XN': '260', 'FDF4UNITS': '59', 'FDDIMCOUNT': '9', 'FDF4XN': '264',
'FDUSER2': '71', 'FDF4APODCODE': '405', 'FDUSER1': '70', 'FDMCFLAG': '135',
'FDFLTORDER': '2', 'FDUSER5': '74', 'FDCOMMENT': '312', 'FDF3QUADFLAG': '51',
'FDUSER4': '73', 'FDTEMPERATURE': '157', 'FDF2APOD': '95', 'FDMONTH': '294',
'FDF4OFFPPM': '483', 'FDF3OFFPPM': '482', 'FDF3CAR': '68', 'FDF4P0': '62',
'FDF4P1': '63', 'FDF1OFFPPM': '481', 'FDF4APOD': '53', 'FDF4X1': '263',
'FDLASTBLOCK': '359', 'FDPLANELOC': '14', 'FDF2FTSIZE': '96',
'FDUSERNAME': '290',
'FDF1X1': '259', 'FDF3CENTER': '81', 'FDF1CAR': '67', 'FDMAGIC': '0',
'FDF2ORIG': '101', 'FDSPECNUM': '219', 'FDF2LABEL': '16', 'FDF2AQSIGN': '64',
'FDF1UNITS': '234', 'FDF2LB': '111', 'FDF4AQSIGN': '477', 'FDF4ZF': '439',
'FDTAU': '199', 'FDF4LABEL': '22', 'FDNOISE': '153', 'FDF3APOD': '50',
'FDF1APODCODE': '414', 'FDF2SW': '100', 'FDF4OBS': '28', 'FDQUADFLAG': '106',
'FDF2TDSIZE': '386', 'FDHISTBLOCK': '364', 'FDSRCNAME': '286',
'FDBASEBLOCK': '361', 'FDF1APODQ2': '421', 'FDF1APODQ3': '422',
'FDF1APODQ1': '420', 'FDF1QUADFLAG': '55', 'FDF3UNITS': '58', 'FDF2ZF': '108',
'FDCONTBLOCK': '360', 'FDDIMORDER4': '27', 'FDDIMORDER3': '26',
'FDDIMORDER2': '25', 'FDDIMORDER1': '24', 'FDF2CAR': '66',
'FDF3APODCODE': '400',
'FDHOURS': '283', 'FDF1CENTER': '80', 'FDF3APODQ1': '401', 'FDF3APODQ2': '402',
'FDF3APODQ3': '403', 'FDBMAPBLOCK': '363', 'FDF2CENTER': '79'}
*/

};

#endif
