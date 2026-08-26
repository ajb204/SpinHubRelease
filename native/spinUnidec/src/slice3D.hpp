/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE3D_HPP
#define SLICE3D_HPP

//#include <iostream>
#include <sstream>
#include <fstream>
#include "slice.hpp"
#include "fftw3.h"
#include "pipeClass.hpp"
#include "general.hpp"

using namespace std;



class slice3D
{
 public:

  //make non-copyable.
  slice3D() = default;
  slice3D(const slice3D&) = delete;
  slice3D& operator=(const slice3D&) = delete;

  
  vector<raw3D>spec;
  double refx,refy,refz; //location of slice in spectrum
  string refn;
  int refi,refj,refk; //get locations of peak (indicies)
  int si,sj,sk; //size x/y
  int size,size2;
  double imin,jmin,kmin,imax,jmax,kmax; //limits
  double sig1,sig2,sig3;
  double noiseVal;
  int indirect=0;//0 for carbon, 1 for proton
  double maxInt;

  int squash_window_i,squash_window_j,squash_window_k;

  vector<peakEntry> peakList;
  int peaks;
  int pipey=0; //read in from nmrPipe? ie is DI alive?
  bool SPARSE=true;
  
#ifdef DOUBLE3D
  double *DI,*DB,*DS; //will hold x,y,I,Band S
  //fft peakshape
  fftw_complex *P,*C;
  fftw_plan p1,pinv;
  double *ivals,*jvals,*kvals; //1D array of x and y vals
  double *pki,*pkj,*pkk;
  double voigt1, voigt2, voigt3; //0 gaus 1 lorentz
  double lor1, lor2, lor3;
  int SIZEMEM=sizeof(double);
#else
  float *DI,*DB,*DS; //will hold x,y,I,Band S
  //fft peakshape
  fftwf_complex *P,*C;
  fftwf_plan p1,pinv;
  float *ivals,*jvals,*kvals; //1D array of x and y vals
  float *pki,*pkj,*pkk;
  float voigt1, voigt2, voigt3; //0 gaus 1 lorentz
  float lor1, lor2, lor3;
  int SIZEMEM=sizeof(float);
#endif



  void SetMem();
  void ReadPipe(pipe &pipefile);
  void InitBlur();
  void InitBlurPeakRestricted();

  void MakeSquare();


#ifdef DOUBLE3D
  int DoIndex(double ref,double *vals,int ii,int p);
#else
  int DoIndex(double ref,float *vals,int ii,int p);
#endif

  //get integer index of closest point in ppm
  void SetIndex();
  void SetIndexRestricted();

#ifdef DOUBLE1D
  void SetBlur(int p,double *ref);
#else
  void SetBlur(int p,float *ref);
#endif


  double ApplyIter();

  void GetPeak();

  void BlankDB();
  void BlankDBfull();
  void CalcSpec();

  void BuildSparseDB(double cutoff);
  void CalcSpecSparse();
  

  void PrintProj();
  void PrintSpec();
  double GetChi2();

  void Squash();
  void Cull(float frac);

  int CountElements();

  void AddIfNew(vector<vector <int> > &peakSum,int j,int k);

  int correlate(string projListFile,FILE *out_pt);
  int correlateRestricted(FILE *out_pt);

  struct SparsePt3D
  {
    int i, j, k, ii;
    double val;
  };
  std::vector<SparsePt3D> sparseDB;

  
  ~slice3D()
  {

    delete [] DI;
    delete [] DS;
    delete [] DB;
    delete [] ivals;
    delete [] jvals;
    delete [] kvals;

    delete [] pki;
    delete [] pkj;
    delete [] pkk;
    
    //fftw variables
    if(P)
      delete [] P;
    if(C)
      delete [] C;

#ifdef DOUBLE3D
    if(p1)
      fftw_destroy_plan(p1);
    if(pinv)
      fftw_destroy_plan(pinv);
#else
    if(p1)
      fftwf_destroy_plan(p1);
    if(pinv)
      fftwf_destroy_plan(pinv);
#endif

    
  }


  
};


#endif
