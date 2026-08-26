/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef DECONMAIN_HPP
#define DECONMAIN_HPP


//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <iostream>
//#include <fstream>
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_fft_complex.h>
//#include <time.h>
//#include <iterator>
//#include <cstring>
//#include <algorithm>
//#include <string>

//#include <complex.h>
//#include <vector>

//#include <sstream>
//#include <queue>
//#include <map>
//#include <random>

//Essential for Linux
//#include <numeric>
//#include <iterator>
//#include <cstring>



#include <iostream>
#include <sstream>
#include <fstream>
#include "slice.hpp"
#include "fftw3.h"
#include "pipeClass.hpp"
#include "general.hpp"

#include<string.h>

#include "slice.hpp"     //basic data objects
#include "pipeClass.hpp" //nmrPipe files
#include "slice2D.hpp"   //for 2D sets (versus peak)
#include "slice3D.hpp"   //for 3D sets
#include "slice1D.hpp"   //for 1D sets (versus peak)
#include "slice4D.hpp"   //for 4D sets



class decon
{
 public:

  string infile;
  string peakfile;
  string baseFile;

  int parsed;


  double rand;
  double dec3d;
  int ncpus;
  double symmy;
  double noiseVal; //noisevalue
  double sig1;     //peak width factor
  double sig2;     //peak width factor
  double sig3;     //peak width factor
  double sig4;     //peak width factor
  double fac;     //factor to multiply width
  double squash;  //region to remove intensitiy from diagonals
  //int indirect=1; //0=C, 1=H. for symmetric spectra, which dimension in peaklist to slice
  int bore=0;  //use 2D peak-list bore mode in 3D
  bool recon=false; //restricted nD reconstruction at supplied full peak positions
  bool pseudo3D=false; //3P: two spectral dimensions plus one real/pseudo axis
  bool pseudo2DFit=false; //2P FIT: one spectral dimension plus one real/pseudo axis
  bool pseudo2DOutput=false; //pure 1D projection job: write Number, PPM, Intensity peak list
  bool enhance=false; //single-pass enhanced-resolution spectrum output (1D-3D only)

  bool FIT=false;
  // Pseudo2D FIT only: after the conventional zero-phase amplitude fit,
  // allow one shared spectral phase parameter per resonance to relax.
  bool FitPhase=false;
  double FitRad=5;
  // Optional explicit 2D FIT extraction radii (ppm).  Presence of either
  // FitF1 or FitF2 in decon.init switches from FitRad-derived radii to these.
  bool FitRadFix=false;
  double FitF1=2.0;
  double FitF2=0.4;

  // Constrain Gaussian/Lorentzian widths during the LM refinement to remain
  // local to the robust guess-and-check solution.  This prevents an almost
  // flat, enormously broad component from acting as an artificial baseline.
  // FUDA-like linewidth mode for 2D fitting.  When true, Gaussian and
  // Lorentzian FWHM are exactly the same within each dimension in both
  // guess/check and LM; g controls only the Gaussian/Lorentzian mixture.
  // Set false in code (or FitWidthRestrict=0 in decon.init) for the original
  // independent Gaussian/Lorentzian linewidth model.
  bool FitWidthRestrict=true;

  
  int STD;
  double squashH;
  double squashC;
  double centrefac=1;

  int iterShow=100; //number of iterations to show results

  // Convergence monitoring. systemIter is deliberately separate from the
  // per-DoRun iteration counter so maxIter retains its existing semantics.
  long long systemIter=0;
  std::ofstream convergenceFile;
  bool convergenceFileOpened=false;
  int convergenceFlushInterval=100;
  void LogConvergenceEvent(const std::string &label);
  
  int squash_window_i,squash_window_j,squash_window_k,squash_window_l;

  int symmode; //symmode
  int mode=1; //1 2 or 3D

  int peaks; //number of peaks to study
  float voigt1=0.0; //0 gaus 1 lorentz
  float voigt2=0.0; //0 gaus 1 lorentz
  float voigt3=0.0; //0 gaus 1 lorentz
  float voigt4=0.0;
  float lor1=0.0;
  float lor2=0.0;
  float lor3=0.0;
  float lor4=0.0;

  int dim; //number of dimensions in spectrum
  int si,sj,sk,sl; //sizes of the 3 dimensions
  double uimin,uimax,ujmin,ujmax,ukmin,ukmax,ulmin,ulmax; //spec limits (ppm)

  //value should depend on precision of machine and FFT type.
  //if FFT values are saved as floats then the conv score is the sum of floats,
  //but saved as a double so shouldnt be too precise. too low, calc will be
  //unstable. too high, might not be enough converged.
  double convVal=1E-7;   //fractional convergence criteria. 
  //double peakThresh=1E-4; //threshold for including peak shape in sparse calc
  int maxIter=10000;      //maximum number of iterations
  int maxIter3D=100;
  int peakN; //number of points above threshold for peak shape function

  double dmax; //make data height
  double thresh;//dmax/noiseval threshold
  string initFile;
  double noiseValSTD; //noisevalue

  //vector<slice1D> sliceLib1D; //library of 2D slices
  //vector<slice2D> sliceLib2D; //library of 2D slices
  //vector<slice3D> sliceLib3D; //library of 2D slices
  //vector<slice3D> sliceLib3Dbore; //library of 3D slices
  //vector<slice4D> sliceLib4D; //library of 2D slices

  deque<slice1D> sliceLib1D; //library of 2D slices
  deque<slice2D> sliceLib2D; //library of 2D slices
  deque<slice3D> sliceLib3D; //library of 2D slices
  deque<slice3D> sliceLib3Dbore; //library of 3D slices
  deque<slice4D> sliceLib4D; //library of 2D slices

  
  vector<peakEntry> peakList;


  //peak shape functions
  //function to return a gaussian at a specific point
//  double Gaus(double x,double x0,double sig){
//    return exp(-pow(x-x0,2.)/(2.*sig*sig));}
//  //function to return a gaussian at a specific point
//  double Lorentz(double x,double x0,double sig){
//    return pow(sig/2,2)/(pow(x-x0,2)+pow(sig/2,2));}
//  //function to return a 2D peak value at a given point
//  double Peak(double x1,double x2,double sig){
//    switch(peaky){
//    case 0:
//      return Gaus(x1,x2,sig);
//      break;
//    case 1:
//      return Lorentz(x1,x2,sig);
//      break;}
//      }


  //read in the file and load data into arrays

  void splash(int argc,char *argv[],bool parallel);
  void parse(string file_name);
  
  void readPeak();
  void readPeakND();
  int CountElements();
  int CountElements(const int p);

  //COREFUNCTIONS
  void calcspec();
  //void calcspec(const int p);

  //no longer userful - using nmrPipe to readin/out
  void PrintSpec();
  //FOR SYMMODE ONLY
  void RunMapBlur();
  //FOR SYMMODE ONLY
  void RunUnMapBlur();
  //signbit gives 1 if negative, 0 if positive
  int signy(double value);

  void DoRun();
  //void DoRun(const int p);
  //Run B_(k+1)=I.B_k / (B conv P)
  double ApplyIter();
  //double ApplyIter(const int p);
  void GetChi2();
  void CullSym();
  //culling in DB space
  void Cull(float frac = 1.0);
  void Squash();
  //void Squash(const int p);
  void correlate(string tag);

  //stick a file in a vector of vectors
  void SetPeaks();


  /**********************************/ //start of 4D
  void Setup4D();
  void calcspec4D();
  /**********************************/ //end of 4D

  /**********************************/ //start of 3D
  void Setup3D();
  void calcspec3D();
  /**********************************/ //end of 3D


  /**********************************/ //start of 2D
  void Setup2Dfrom4D();
  void Setup1Dpure(pipe &pipefile);
  void SetupSTD(pipe &pipefile);
  void Setup2Dpure(pipe &pipefile);
  void calcspec2D();
  /**********************************/ //start of 2D



  /**********************************/ //start of 1D
  void Setup1Dfrom3D(pipe &pipefile);
  void calcspec1D();
  /**********************************/ //end of 1D



  //perform 3D deconvolution
  //start with taking 1D splices
  //then perform option 3D decon
  //reconstruct 3D spectrum, slice and project.


  //a Binit file can contain a list of ppms: this means don't initialise in all places.
  //1. if does not equal False, read file and use this to set Binits for deltas.

  //also: add baseline mode.
  //if we are baselining, setup window.
  //if doing a protein baseline, need to read in protein baseline file.
  //need to fft and save for us.

  //when calculating spectrum, add contribution from convolution of baseline
  //need to figure out which deltas to use for baseline.


  //in the UDC iterator, need to have option for averaging I when dealing with baseline./
  //udc needs to know if we have a baseline in play or not.

  //strategy: lets add the baseline mode as an addition.
  //1. read file. make PSF. (note addbackground function, with protein centre) //DONE
  //2. read window to use, and corresponding mask.                             //not done. currently single val
  //3. in simspec, add baseline component if needed.                           //DONE
  //4. In UDC, if mode is on, increment baseline deltas, using window.         //

  void ProtocolSTD(int argc, char *argv[]);

  
  // 31st March 2025: 
  //synced up with protocol1D in uSTA project
  //doesn't give exactly the same values as decon_darwin_i386_uSTA
  //differences: we now used floats, uSTA seems to use doubles
  //some differences in the slice1D functions.
  //disconcerting.
     
  void Protocol1D(int argc,char *argv[]);

  //perform 3D deconvolution
  //start with taking 1D splices
  //then perform option 3D decon
  //reconstruct 3D spectrum, slice and project.
  void Protocol2D(int argc,char *argv[]);
  void Protocol3D(int argc,char *argv[]);
  void Protocol3P(int argc,char *argv[]);
  void Protocol2PFit(int argc,char *argv[]);
  void Protocol4D(int argc,char *argv[]);

}; //end of class



// parse definitions are in deconParse.hpp

#endif
