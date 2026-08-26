/* 
   re-shuffle.cpp
  
   Re-shuffle data from a 2D diffusion experiment

   Want to get both the up and downfield peaks for a 2D diffusion expt.
   Plus we want to correct for a 90o phase error in the pulse sequence.
   
   D Flemming Hansen; February 2006
   (C) flemming@pound.med.utoronto.ca 

   - Allow for both integer and float format in the varian fid.
     Flemming, March 20th, 2006

N.B. using int format has bit size of 4. 
Do not use long int if both 32 and 64 bit processors are required
AJB. 9/6/08

Apr 2013 Generalised by A.Baldwin


*/

#define PTYPE 0    // 0 for the upfield peak; 1 for the downfield; 2 for both

#define VARIANINT  //Can take long integer 'VARIANINT' (data from old sun, s.a. pence)
                   //         short float  'VARIANFLT' (data from newer machines, s.a. curie )
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdio>
#include <iostream>


void swap4(int n, float f[]) {
  /* Swap a 4 byte float .. shitty SUN machines! */
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
  /* Swap a 4 byte float .. shitty SUN machines! */
//CHANGE INT TO LONG INT FOR 32 BIT MACHINES

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
//
// Here we go .. 
int main(int argc, char* argv[] )  
{

  if(argc!=7)
    {
      printf("USAGE: %s NP NI NZ mode(0-ncyc,ip_flg,phase 1-ip_flg,phase,ncyc) [FixedFID] [VarianFID] \n\n",argv[0]);
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

  static int Data[maxNZ][2][2][maxNP]; //Need to use static to get extra me
  static int Data2[maxNZ][2][2][maxNP]; //Need to use static to get extra mem
  static int NewData[5][maxNP];



  char prehead[50];
  char Head[NZ][4][2][50];
  char prehead2[50];
  char Head2[NZ][4][2][50];

  //  std::string oName0=std::string("fid.final");
  FILE* One   = fopen("fid","r");
  FILE* fo   = fopen(oName.c_str(),"w");

  if (One==NULL) {
    std::cerr << " Cannot open input file: " << "fid" << std::endl;
    std::cerr << "                    PROGRAM ABORTED" << std::endl;
    exit(1);
  };
  if ( fo == NULL ){
    std::cerr << " Cannot open output file: " << oName << std::endl;
    std::cerr << "                    PROGRAM ABORTED" << std::endl;
    exit(1);
  };
  int np=     NP;
  int ni=     NI;
  int NoGrads=NZ;
  /* read and write overall header */
  fread(prehead,1,32,One);
  fwrite(prehead,1,32,fo);  
  for (int d2=0;d2<ni;d2++){ 
    //
    //Read the data (In this example the array was set to: "gzlvl,phase,IP_flg"
    for ( int gzlvl=0;gzlvl<NoGrads;gzlvl++){
      for ( int IP_flg=0;IP_flg<2;IP_flg++){
	for ( int phase=0;phase<2;phase++) {
	  //
	  //Header
	  fread(Head[gzlvl][IP_flg][phase],1,28,One);
	  //
	  //Data
	  fread(Data[gzlvl][IP_flg][phase],4,np,One);
	}; 
      };
    };
    //
    // Write the data (in the right order!)
    
      
    for ( int gzlvl=0;gzlvl<NoGrads;gzlvl++){
      for ( int IP_flg=0;IP_flg<2;IP_flg++){	
        for ( int phase=0;phase<2;phase++){		

	  //
	  // Peak selection
	  if ( PTYPE != 2 ) {
	    if ( PTYPE != IP_flg ) {
	      continue;
	    };
	  };
	  //
	  // Store data in temporary array
	  for ( int i=0;i<np;i++){
	    NewData[0][i]=  Data[gzlvl][0][phase][i];
	    NewData[1][i]=  Data[gzlvl][1][phase][i];
	  };
	  //
	  swap4(np,NewData[0]);
	  swap4(np,NewData[1]);
	  //
	  // 90 degrees phase correction for IP data .. complex number are beautiful :-p
	  for ( int i=0;i<np/2;i++){	    	    
	    NewData[2][2*i]=   -1*NewData[0][2*i+1];
	    NewData[2][2*i+1]=  1*NewData[0][2*i];	    
	  };
	  // No phase correction for AP data
	  for ( int i=0;i<np;i++){
	    NewData[3][i]=NewData[1][i];
	  };	  
          #ifdef VARIANINT
            //CHANGE INT TO LONG INT FOR 32 BIT MACHINES
	    int pf=1;
          #endif
	  #ifdef VARIANFLT 
	    double pf=1.;
          #endif
	  if ( IP_flg == 1 ) {
	    pf=-1*pf;
	  };
	  for ( int i=0;i<np;i++ ) {
	    //NewData[4][i]=NewData[2][i]+pow(-1,IP_flg)*NewData[3][i];
	    NewData[4][i]=NewData[2][i]+pf*NewData[3][i];
	  };
	  swap4(np,NewData[4]);	
	  //
	  fwrite(Head[gzlvl][IP_flg][phase],1,28,fo);
	  fwrite(NewData[4],4,np,fo);
	};
      };
    };
  };
  fclose(One);
  fclose(fo);
  return 0;
};
