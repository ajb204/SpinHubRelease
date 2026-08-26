/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE4D_HPP
#define SLICE4D_HPP

//#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "slice.hpp"
#include "fftw3.h"
#include "pipeClass.hpp"
#include "general.hpp"



/*4D FTs running at floating point precision to preserve memory*/
//1st April 2019: DS removed to save memory: DB now does all the work
//works because we don't need to save DBs in this mode.

class slice4D
{
 public:

  //make non-copyable.
  slice4D() = default;
  slice4D(const slice4D&) = delete;
  slice4D& operator=(const slice4D&) = delete;

  
  vector<raw4D>spec;
  double refx,refy,refz; //location of slice in spectrum
  string refn; //name of slice
  int refi,refj,refk,refl; //get locations of peak (indicies)
  int si,sj,sk,sl; //size x/y
  int size,size2;
  double imin,jmin,kmin,lmin,imax,jmax,kmax,lmax; //limits
  double noiseVal;
  bool SPARSE=true;
  int *n;
  vector<peakEntry> peakList;
  int peaks;
 
#ifdef DOUBLE4D
  double *ivals,*jvals,*kvals,*lvals; //1D array of x and y vals
  double sig1,sig2,sig3,sig4;
  double lor1,lor2,lor3,lor4;
  double voigt1, voigt2, voigt3, voigt4; //0 gaus 1 lorentz
  //double *DI,*DB,*DS; //will hold x,y,I,Band S
  //float *DB,*DS;
  double *DB;
  //fft peakshape
  fftw_complex *P,*C;
  fftw_plan p2,p1,pinv;
  int SIZEMEM = sizeof(double);
#else
  float *ivals,*jvals,*kvals,*lvals; //1D array of x and y vals
  float sig1,sig2,sig3,sig4;
  float lor1,lor2,lor3,lor4;
  float voigt1, voigt2, voigt3, voigt4; //0 gaus 1 lorentz
  //double *DI,*DB,*DS; //will hold x,y,I,Band S
  //float *DB,*DS;
  float *DB;
  //fft peakshape
  fftwf_complex *P,*C;
  fftwf_plan p2,p1,pinv;
    int SIZEMEM = sizeof(float);
#endif

  void Read();


  void MakeSquare();

 #ifdef DOUBLE4D
  int DoIndex(double ref,double *vals,int ii);
 #else
  int DoIndex(double ref,float *vals,int ii);
 #endif
 //get integer index of closest point in ppm
  void SetIndex();


  void BlankDBFull(); //completely blank DB

  void BlankDB(); //blank just the 2D slice


  int CountElements();


#ifdef DOUBLE2D
  void SetBlur(int p,double *DB_2D); //copy in 2D blurs to 4D
  #else
  void SetBlur(int p,float *DB_2D); //copy in 2D blurs to 4D
 #endif

#ifdef DOUBLE2D
  void ReadBlur(int p,double *DS_2D); //copy blurs back to 2D DS
  #else
  void ReadBlur(int p,float *DS_2D); //copy blurs back to 2D DS
  #endif


  void SetPeak(slice4D inst);



  void GetPeak();


  void RunFFT_1();


  void RunFFT_inv();


  void CalcSpec();
  void BuildSparseDB(double cutoff);

  struct SparsePt4D
  {
    int i, j, k, l, ii;
    double val;
  };
  std::vector<double> pki, pkj, pkk, pkl;
  std::vector<SparsePt4D> sparseDB;


  void PrintSpec();

  ~slice4D()
  {

    delete [] ivals;
    delete [] jvals;
    delete [] kvals;
    delete [] lvals;

    delete [] DB;


    if(P)
      delete [] P;
    if(C)
      delete [] C;
    
#ifdef DOUBLE4D
    if(p2)
      fftw_destroy_plan(p2);
    if(p1)
      fftw_destroy_plan(p1);
    if(pinv)
      fftw_destroy_plan(pinv);
#else
    if(p2)
      fftwf_destroy_plan(p2);
    if(p1)
      fftwf_destroy_plan(p1);
    if(pinv)
      fftwf_destroy_plan(pinv);
#endif

    
    delete [] n;

  }

};

#endif
