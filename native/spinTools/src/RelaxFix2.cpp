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
      Thus, if array='ncyc,flg,phase' define MODE 0
            if array='ncyc,phase,flg' define MODE 1 
            if array='phase,ncyc,flg' define MODE 2

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
 
 
  if(argc!=8)
    {
      printf("USAGE: %s NP NI NZ FG mode(0,1,2) [FixedFID] [VarianFID] \n\n",argv[0]);
    exit(2);
    }
 
  int NP=atoi(argv[1]);
  int NI=atoi(argv[2]);
  int NZ=atoi(argv[3]);
  int FG=atoi(argv[4]);
  int MODE=atoi(argv[5]);

  std::string oName=std::string(argv[6]);
  std::string iName=std::string(argv[7]);

  std::cout << "NP     : " << NP << std::endl;
  std::cout << "NI     : " << NI << std::endl;
  std::cout << "NZ     : " << NZ << std::endl;
  std::cout << "FG     : " << FG << std::endl;
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
#define maxFG 2      // Number of flags

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
  if(FG>maxFG){
    std::cout << "FG is too large! Increase limit" << std::endl;
    exit(2);
  }

  static int Data[maxNI][maxNZ][maxFG][2][maxNP]; //Need to use static to get extra mem

  char prehead[50];
  char Head[NI][NZ][FG][2][50];
  
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
  int fg=      FG;
  int mode = MODE;
  /* read and write overall header */
  fread(prehead,1,32,One);
  fwrite(prehead,1,32,fo);  
  for (int d2=0;d2<ni;d2++){ 
    if ( mode == 0 ) {
      for ( int ncyc=0;ncyc<nz;ncyc++){
	for ( int flg=0;flg<fg;flg++){
	  for ( int phase=0;phase<2;phase++) {
	    //
	    //Header
	    fread(Head[d2][ncyc][flg][phase],1,28,One);
	    //
	    //Data
	    fread(Data[d2][ncyc][flg][phase],4,np,One);
	  }; 
	};
      };
    } else if ( mode == 1 ) {
      for ( int ncyc=0;ncyc<nz;ncyc++){
	for ( int phase=0;phase<2;phase++) {
	  for ( int flg=0;flg<fg;flg++){
	    //
	    //Header
	    fread(Head[d2][ncyc][flg][phase],1,28,One);
	    //
	    //Data
	    fread(Data[d2][ncyc][flg][phase],4,np,One);
	  }; 
	};
      };
    } else if ( mode == 2 ) {
      for ( int phase=0;phase<2;phase++) {
	for ( int ncyc=0;ncyc<nz;ncyc++){
	  for ( int flg=0;flg<fg;flg++){
	    //
	    //Header
	    fread(Head[d2][ncyc][flg][phase],1,28,One);
	    //
	    //Data
	    fread(Data[d2][ncyc][flg][phase],4,np,One);
	  }; 
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
    for ( int flg=0;flg<fg;flg++){    
      for ( int d2=0;d2<ni;d2++){ 
	for ( int phase=0;phase<2;phase++) {
	  fwrite(Head[d2][ncyc][flg][phase],1,28,fo);
	  fwrite(Data[d2][ncyc][flg][phase],4,np,fo);
	};
      };
    };
  };
  fclose(One);
  fclose(fo);
  return 0;
};
