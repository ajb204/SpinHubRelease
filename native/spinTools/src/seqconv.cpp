/*
  SEQUENCE CONVERTER

  Do useful things with amino acid sequence data

  A.Baldwin June2006
  ajb204@nmrmed.med.utoronto.ca

Reads in file of the form: residue number, amino acid (3 letter)
Outputs residue number, amino acid (3 letter), amino acid (1 letter)

*/
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdio>
#include <iostream>
#include <string>


using namespace std;
using std::string;


int main(int argc, char* argv[])
{
  if ( argc != 2 ) 
    {
      std::cerr << " USAGE: \n";
      std::cerr << " " << argv[0] << " inputfile" << std::endl;
      std::cerr << std::endl;
      exit(2); 
    }
  std::string iName=std::string(argv[1]);
  
  FILE* One  = fopen(iName.c_str(),"r");
  if (One==NULL) 
    {
      std::cerr << " Cannot open input file: " << iName << std::endl;
      std::cerr << "                    PROGRAM ABORTED" << std::endl;
      exit(1);
    };

  fclose(One);



  /*
  ofstream myfile ("example.txt");
  if (myfile.is_open())
  {
    myfile << "This is a line.\n";
    myfile << "This is another line.\n";
    myfile.close();
  }
  else cout << "Unable to open file";
  
  
  cout << "READ IN\n";

   string line;
   ifstream myfile2 (iName.c_str());
   if (myfile2.is_open())
     {
       while (! myfile2.eof() )
	 {
	   getline (myfile2,line);
	   cout << line << endl;
	 }
       myfile2.close();
     }
   
   else cout << "Unable to open file"; 
   
   return 0;
  */








  int g=0;  
  string amino[500];
  int seq[500];
  ifstream ANDY20 (iName.c_str());
  while (!ANDY20.eof() )
    {
      g++;
      ANDY20 >> seq[g] >> amino[g];
      //      cout << seq[g] << amino[g] << endl;
    }
  ANDY20.close();
  
  string amino2[500];
  cout << endl;  
  int i=0;
  for(i=0;i<g;i++)
    {
      if(amino[i]=="ALA")
	amino2[i]="A";
      if(amino[i]=="ARG")
	amino2[i]="R";
      if(amino[i]=="ASN")
	amino2[i]="N";
      if(amino[i]=="ASP")
	amino2[i]="D";
      if(amino[i]=="CYS")
	amino2[i]="C";
      if(amino[i]=="GLN")
	amino2[i]="Q";
      if(amino[i]=="GLU")
	amino2[i]="E";
      if(amino[i]=="GLY")
	amino2[i]="G";
      if(amino[i]=="HIS")
	amino2[i]="H";
      if(amino[i]=="ILE")
	amino2[i]="I";
      if(amino[i]=="LEU")
	amino2[i]="L";
      if(amino[i]=="LYS")
	amino2[i]="K";
      if(amino[i]=="MET")
	amino2[i]="M";
      if(amino[i]=="PHE")
	amino2[i]="F";
      if(amino[i]=="PRO")
	amino2[i]="P";
      if(amino[i]=="SER")
	amino2[i]="S";
      if(amino[i]=="THR")
	amino2[i]="T";
      if(amino[i]=="TRP")
	amino2[i]="W";
      if(amino[i]=="TYR")
	amino2[i]="Y";
      if(amino[i]=="VAL")
	amino2[i]="V";
      cout << seq[i] << "\t" << amino[i] << "\t" << amino2[i] << "\t" << endl;
      if(amino2[i]=="")
	cout << "problem with residue " << i << endl;
    }
  return 0;
  
}
