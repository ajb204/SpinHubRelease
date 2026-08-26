/* 
   AddFIDvarian.cpp
  
   Take two FIDs and combine.
   
   (c) A . Baldwin  March 2024

   
   Modified by DF Hansen, August 2006
   Modified by AJ Baldwin, Feb 2013 to take command line args

    - Include define flg, so that also 'varian-int' format are accepted.
      (Later this might be read from the varian header in the fid; along with NP,NI,and NZ )
      Thus, if dps='r' define VARIANINT
            if dps='i' define VARIANFLT
    - Include a mode,
      Thus, if array='ncyc,phase' define MODE 0
            if array='phase,ncyc' define MODE 1 

*/


#define VARIANINT
//////////////////////////////////////////////////////////////////////////////
//
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdio>
#include <iostream>
//
// Here we go ..

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
 
 
  if(argc!=6)
    {
      printf("USAGE: %s NP NFids [VarianFID1] [VarianFID2] [NewFID]\n\n",argv[0]);
      exit(2);
    }
  
  int NP=atoi(argv[1]);
  int NI=atoi(argv[2]);
  std::string iName1=std::string(argv[3]);
  std::string iName2=std::string(argv[4]);
  std::string oName=std::string(argv[5]);


  std::cout << "NP     : " << NP << std::endl;
  std::cout << "NI     : " << NI << std::endl;
  std::cout << "Infile1 : " << iName1 << std::endl;
  std::cout << "Infile2 : " << iName2 << std::endl;
  std::cout << "Outfile: " << oName << std::endl;


#ifdef VARIANINT 
  //static int Data[NI][NZ][2][MAXDATA];
  int Data1[NP];
  int Data2[NP];
  int NewData[NP];
#endif
#ifdef VARIANFLT 
  float Data1[NP];
  float Data2[NP];
  float NewData[NP];
  //static float Data[NI][NZ][2][MAXDATA];
#endif
  
  //static int Data[maxNI][maxNZ][2][maxNP]; //Need to use static to get extra mem

  char prehead1[50];
  char prehead2[50];
  char Head1[50];
  char Head2[50];
  
  FILE* OneA  = fopen(iName1.c_str(),"r");
  FILE* OneB  = fopen(iName2.c_str(),"r");
  FILE* fo   = fopen(oName.c_str(),"w");

  if (OneA==NULL) {
    std::cerr << " Cannot open input file: " << iName1 << std::endl;
    std::cerr << "                    PROGRAM ABORTED" << std::endl;
    exit(1);
  };
  if (OneB==NULL) {
    std::cerr << " Cannot open input file: " << iName2 << std::endl;
    std::cerr << "                    PROGRAM ABORTED" << std::endl;
    exit(1);
  };
  if ( fo == NULL ){
    std::cerr << " Cannot open output file: " << oName << std::endl;
    std::cerr << "                    PROGRAM ABORTED" << std::endl;
    exit(1);
  };
  int np=      NP;
  int ni=      NI;
  /* read and write overall header */
  fread(prehead1,1,32,OneA);
  fread(prehead2,1,32,OneB);
  fwrite(prehead1,1,32,fo);  
  for (int d2=0;d2<ni;d2++)
    {
	  //
	  //Header
	  fread(Head1,1,28,OneA);
	  fread(Head2,1,28,OneB);
	  //
	  //Data
	  fread(Data1,4,np,OneA);
	  fread(Data2,4,np,OneB);
	  
	  swap4(np,Data1);
	  swap4(np,Data2);
	  
	  for ( int i=0;i<np;i++)
	    {
	      NewData[i]=Data1[i]+Data2[i];
	    };
	  swap4(np,NewData);	
	  //
	  // Write the data (in the right order!)
	  fwrite(Head1,1,28,fo);
	  fwrite(NewData,4,np,fo);
    };
  fclose(OneA);
  fclose(OneB);
  fclose(fo);
  return 0;

};
