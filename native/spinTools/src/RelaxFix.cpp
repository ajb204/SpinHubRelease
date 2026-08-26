/* 
   RelaxFix.cpp
  
   Re-shuffle data from a 2D relaxation experiment
   
   D Flemming Hansen;  September 2005
   (C) flemming@pound.med.utoronto.ca
   
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


  //#ifdef VARIANINT 
  //static int Data[NI][NZ][2][MAXDATA];
  //#endif
  //#ifdef VARIANFLT 
  //static float Data[NI][NZ][2][MAXDATA];
  //#endif

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

  static int Data[maxNI][maxNZ][2][maxNP]; //Need to use static to get extra mem

  char prehead[50];
  char Head[NI][NZ][2][50];
  
  FILE* One  = fopen(iName.c_str(),"r");
  FILE* fo   = fopen(oName.c_str(),"w");

  if (One==NULL) {
    std::cerr << " Cannot open input file: " << iName << std::endl;
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
  int nz=      NZ;
  int mode = MODE;
  /* read and write overall header */
  fread(prehead,1,32,One);
  fwrite(prehead,1,32,fo);  
  for (int d2=0;d2<ni;d2++){ 
    if ( mode == 0 ) {
      for ( int ncyc=0;ncyc<nz;ncyc++){
	for ( int phase=0;phase<2;phase++) {
	  //
	  //Header
	  fread(Head[d2][ncyc][phase],1,28,One);
	  //
	  //Data
	  fread(Data[d2][ncyc][phase],4,np,One);
	}; 
      };
    } else if ( mode == 1 ) {
      for ( int phase=0;phase<2;phase++) {
	for ( int ncyc=0;ncyc<nz;ncyc++){
	  //
	  //Header
	  fread(Head[d2][ncyc][phase],1,28,One);
	  //
	  //Data
	  fread(Data[d2][ncyc][phase],4,np,One);
	}; 
      };
    } else {
      std::cerr << " MODE=" << mode << " is not allowed!\n";
      std::cerr << " MODE must be 0 or 1" << std::endl;
    };
  };
  //
  // Write the data (in the right order!)
  for ( int ncyc=0;ncyc<nz;ncyc++){    
    for ( int d2=0;d2<ni;d2++){ 
      for ( int phase=0;phase<2;phase++) {
	fwrite(Head[d2][ncyc][phase],1,28,fo);
	fwrite(Data[d2][ncyc][phase],4,np,fo);
      };
    };
  };
  fclose(One);
  fclose(fo);
  return 0;
};
