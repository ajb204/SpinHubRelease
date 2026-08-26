/* 
   re-shuffle.cpp
  
   Combine two fids into one and take linear combinations of both.
   Setup for analysing methyl groups.
   
   February 2006
   (c) D.F. Hansen; flemming@pound.med.utoronto.ca

   September 2006   Modified by D.F. Hansen
   Include define flg, so that also 'varian-int' format are accepted.
   (Later this might be read from the varian header in the fid; along with NP,NI,etc.)
   Thus, if dps='r' define VARIANINT
         if dps='i' define VARIANFLT

   July 2007 Modified by A.Baldwin
   Feb 2013 Generalised by A.Baldwin

   Taking the linear combination of CH3

*/

#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdio>
#include <iostream>

#define VARIANINT

//////////////////////////////////////////////////////////////////////////


void swap4(int n, float f[]) {
  /* Swap a 4 byte float */
  int i;
  union swap4 { float f; char s[4]; } in, out;
  
  for(i=0;i<n;i++) {
    in.f = f[i];
    out.s[0] = in.s[3];
    out.s[1] = in.s[2];
    out.s[2] = in.s[1];
    out.s[3] = in.s[0];
    f[i] = out.f;
  }
  return;
  };

void swap4(int n, int f[]) {
  /* Swap a 4 byte int */
  int i;
  union swap4 { int f; char s[4]; } in, out;
  
  for(i=0;i<n;i++) {
    in.f = f[i];
    out.s[0] = in.s[3];
    out.s[1] = in.s[2];
    out.s[2] = in.s[1];
    out.s[3] = in.s[0];
    f[i] = out.f;
  }
  return;
};

int main(int argc, char* argv[] )  {

  if(argc!=7)
    {
      printf("USAGE: %s NP NI NZ mode(0-ncyc,phase 1-phase,ncyc) [FixedFID] [VarianFID] \n\n",argv[0]);
      exit(2);
    }
 
  int NP=atoi(argv[1]);
  int NI=atoi(argv[2]);
  int NZ=atoi(argv[3]);
  int MODE=atoi(argv[4]);

  std::string oName=std::string(argv[5]);
  std::string iName=std::string(argv[6]);

  std::cout << "NP     : " << NP << std::endl;
  std::cout << "NI     : " << NI << std::endl;
  std::cout << "NZ     : " << NZ << std::endl;
  std::cout << "Mode   : " << MODE << std::endl;
  std::cout << "Infile : " << iName << std::endl;
  std::cout << "Outfile: " << oName << std::endl;

#define maxNP 2000    // Number of points in the direct dimension
#define maxNI 200      // Number of indirect increments to be used 
#define maxNZ 50     // Number of z array 

  if(NP>maxNP){
    std::cout << "NP is too large! Increase limit" << std::endl;
    exit(2);
  }
  if(NI>maxNI){
    std::cout << "NI is too large! Increase limit" << std::endl;
    exit(2);
  }
  if(NZ>maxNZ){
    std::cout << "NZ is too large! Increase limit" << std::endl;
    exit(2);
  }

  static int Data[maxNZ][4][2][maxNP]; //Need to use static to get extra me
  static int Data2[maxNZ][4][2][maxNP]; //Need to use static to get extra mem
  static int NewData[8][maxNP];

  /*  #define MAXDATA NP
  #ifdef VARIANINT
    static int Data[NZ][4][2][MAXDATA];
    static int Data2[NZ][4][2][MAXDATA];
    static int NewData[8][MAXDATA];
  #endif
#ifdef VARIANFLT
    static float Data[NZ][4][2][MAXDATA];
    static float Data2[NZ][4][2][MAXDATA];
    static float NewData[8][MAXDATA];

#endif
  */


  char prehead[50];
  char Head[NZ][8][2][50];
  char prehead2[50];
  char Head2[NZ][8][2][50];

  std::string oName0=std::string("fid.0");
  std::string oName1=std::string("fid.1");
  std::string oName2=std::string("fid.2");
  std::string oName3=std::string("fid.3");

  FILE* A   = fopen("fid","r");

  FILE* fo0   = fopen(oName0.c_str(),"w");
  FILE* fo1   = fopen(oName1.c_str(),"w");
  FILE* fo2   = fopen(oName2.c_str(),"w");
  FILE* fo3   = fopen(oName3.c_str(),"w");

  if (A==NULL ) {
    std::cerr << "Cannot open file" << std::endl;
    exit(1);
  };
  int np=    NP;
  int ni=    NI;
  int nz=    NZ;

  /* read and write overall header */
  fread(prehead,1,32,A);
  fwrite(prehead,1,32,fo0);
  rewind(A);
  fread(prehead,1,32,A);
  fwrite(prehead,1,32,fo1);
  rewind(A);
  fread(prehead,1,32,A);
  fwrite(prehead,1,32,fo2);
  rewind(A);
  fread(prehead,1,32,A);
  fwrite(prehead,1,32,fo3);

    //DEAL WITH FIRST FID
  for (int d2=0;d2<ni;d2++){ 
    //
    //Read the data
    for ( int gzlvl=0;gzlvl<nz;gzlvl++){
      for ( int sel=0;sel<4;sel++) {      
	for ( int phase=0;phase<2;phase++) {      
	  //
	  //Header
	  fread(Head[gzlvl][sel][phase],1,28,A);
	  //Data
	  fread(Data[gzlvl][sel][phase],4,np,A);
	  //
	}; 
      };
    };
    //
    // Write the data
    for ( int gzlvl=0;gzlvl<nz;gzlvl++){
      for ( int phase=0;phase<2;phase++){
	//
	// Store data in temporary array
	for ( int j=0;j<4;j++ ){
	  for ( int i=0;i<np;i++){
	    NewData[j][i]=  Data[gzlvl][j][phase][i];
	  };
	  swap4(np,NewData[j]);
	};
	
	for ( int i=0;i<np;i++){
	  NewData[4][i]=NewData[0][i]+2*NewData[1][i]+2*NewData[2][i] + NewData[3][i];
	  NewData[5][i]=NewData[0][i]-2*NewData[1][i]+2*NewData[2][i] - NewData[3][i];
	  NewData[6][i]=NewData[0][i] + NewData[1][i]  - NewData[2][i] - NewData[3][i];
	  NewData[7][i]=NewData[0][i] - NewData[1][i]  - NewData[2][i] + NewData[3][i];
	};
	swap4(np,NewData[4]);	
	swap4(np,NewData[5]);	
	swap4(np,NewData[6]);	
	swap4(np,NewData[7]);	
	
	fwrite(Head[gzlvl][0][phase],1,28,fo0);
	fwrite(Head[gzlvl][1][phase],1,28,fo1);
	fwrite(Head[gzlvl][2][phase],1,28,fo2);
	fwrite(Head[gzlvl][3][phase],1,28,fo3);
	
	fwrite(NewData[4],4,np,fo0);
	fwrite(NewData[5],4,np,fo1);
	fwrite(NewData[6],4,np,fo2);
	fwrite(NewData[7],4,np,fo3);	
      }
    }
    
    
  };
  
  
  fclose(A);
  fclose(fo0);
  fclose(fo1);
  fclose(fo2);
  fclose(fo3);
  return 0;
};
