/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE1DUSTA_HPP
#define SLICE1DUSTA_HPP

//#include <iostream>
#include <sstream>
#include <fstream>
#include "slice.hpp"
#include "fftw3.h"
#include "pipeClass.hpp"
#include "general.hpp"
#include "slice1D.hpp"


class slice1Dusta
{
 public:

  //make non-copyable.
  slice1Dusta() = default;
  slice1Dusta(const slice1Dusta&) = delete;
  slice1Dusta& operator=(const slice1Dusta&) = delete;

  
  vector<raw1D>spec;
  double refx,refy; //location of slice in spectrum
  string refn;
  int refp,refi,refj; //get locations of peak (indicies)
  int si; //size x/y
  int size;
  double imin,imax; //limits
  double sig1; //sig1 will be the 'core' dimension
  int symmode;
  int indirect=0;  //0 take carbon, 1 take proton
  vector<peakEntry> peakList;
  int peaks;
  int squash_window_i;
  double noiseVal=0;

 
  bool BASE=false;

#ifdef DOUBLE1D
 //fft peakshape
  fftw_complex *P,*C;
  fftw_plan p1,pinv;
  double *ivals; //1D array of x and y vals

 double voigt1; //0 Gaus 1 Lorentz
  double lor1;
  double *DI,*DB,*DS,*DBA,*DBR, *DIb2; //will hold x,y,I,Band S
 
  double *DIb,*DBb,*DSb,baseCentre;  //hold the baseline
  fftw_complex *Pb,*Cb;
  fftw_plan p1b,pinvb;
  int SIZEMEM=sizeof(double);
  #else
 //fft peakshape
  fftwf_complex *P,*C;
  fftwf_plan p1,pinv;
  float *ivals; //1D array of x and y vals


 float voigt1; //0 Gaus 1 Lorentz
  float lor1;
  float *DI,*DB,*DS,*DBA,*DBR, *DIb2; //will hold x,y,I,Band S
  float *DIb,*DBb,*DSb,baseCentre;  //hold the baseline
  fftwf_complex *Pb,*Cb;
  fftwf_plan p1b,pinvb;
   int SIZEMEM=sizeof(float);
  #endif

  int windowB,bval;
  vector <int> windowBvals;

  double maxInt,maxIntBase;

  double ApplyIter();

  double ApplyIterBase();

  double GetChi2();


  void SetBinit(string initFile);
  void Read();


  void ReadPipe(pipe &pipefile);


  void ReadPipe();
  void ReadBase(pipe &pipefile);

  void MakeSquare();

  void MakeSquareBase();

#ifdef DOUBLE1D
  int DoIndex(double ref,double *vals,int ii);
#else
  int DoIndex(float ref,float *vals,int ii);
#endif


 //get integer index of closest point in ppm
  void SetIndex();

  void BlurFrom3D(float *ref);
  void MapBlur();
  void UnMapBlur();


  void InitBlur();

  void InitBlurBase();
  void SetBlur(double *ref);

  //function to return a gaussian at a specific point


  void SetPeak(const slice1D &inst);

  void GetPeak();

  void GetPeakBase();

  void CalcSpec();


  void PrintSpec();
  void PrintSpecPure();
  void Squash();

  void Cull(float frac = 1.0);

  void CullCentral(double Lim);

  
  int CountElements();

  int CountElementsSym();

  int correlate(FILE *out_pt);

  int correlateBase(FILE *out_pt);

    ~slice1Dusta()
  {
    delete [] ivals;
    delete [] DI;
    delete [] DB;
    delete [] DS;
    if(symmode)
      {
	delete [] DBA;
	delete [] DBR;
      }

#ifdef DOUBLE1D
    if(P)
      delete [] P;
    if(C)
      delete [] C;
    if(p1)
      fftw_destroy_plan(p1);
    if(pinv)
      fftw_destroy_plan(pinv);
#else
    if(P)
      delete [] P;
    if(C)
      delete [] C;
    if(p1)
      fftwf_destroy_plan(p1);
    if(pinv)
      fftwf_destroy_plan(pinv);
#endif

    if(BASE)
      {
	delete [] DIb;
	delete [] DBb;
	delete [] DSb;
	delete [] DIb2;

#ifdef DOUBLE1D
	if(Pb)
	  delete [] Pb;
	if(Cb)
	  delete [] Cb;
	if(p1b)
	  fftw_destroy_plan(p1b);
	if(pinvb)
	  fftw_destroy_plan(pinvb);
#else
	if(Pb)
	  delete [] Pb;
	if(Cb)
	  delete [] Cb;
	if(p1b)
	  fftwf_destroy_plan(p1b);
	if(pinvb)
	  fftwf_destroy_plan(pinvb);
#endif
      }
  }
  
  
};


#endif
