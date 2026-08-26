/*************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE2D_HPP
#define SLICE2D_HPP

//#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "slice.hpp"
#include "fftw3.h"
#include "pipeClass.hpp"
#include "general.hpp"


#include <numeric>

using namespace std;


struct FitPeak2DLocal
{
    int i, j;          // integer indices
    int ii;            // linear index

    int group;         // fitting group number

    double x;          // ppm position
    double y;

    double raw;        // original DB intensity
    double intensity;  // fitted intensity
    double fitted;     // fitted value
    std::string name;   // peak name for per-peak FUDA-style output

    // dimension 1
    double sig1;
    double lor1;
    double voigt1;

    // dimension 2
    double sig2;
    double lor2;
    double voigt2;
};

struct FitGroup2DLocal
{
    std::vector<int> members;

    int minI, maxI;
    int minJ, maxJ;
};



class slice2D
{
 public:

  //make non-copyable.
  slice2D() = default;
  slice2D(const slice2D&) = delete;
  slice2D& operator=(const slice2D&) = delete;


  
  vector<raw2D>spec;
  double refx,refy; //location of slice in spectrum
  string refn;
  int refp,refi,refj; //get locations of peak (indicies)
  int si,sj; //size x/y
  int size;
  double imin,jmin,imax,jmax; //limits
  int indirect=0; //0 for carbon, 1 for proton
  bool SPARSE=true;
  
#ifdef DOUBLE2D
 double *ivals,*jvals; //1D array of x and y vals
  double sig1,sig2,sig3;
  double lor1, lor2;
 
  double *DI=nullptr,*DB=nullptr,*DS=nullptr,*DBA=nullptr,*DBR=nullptr; //will hold x,y,I,Band S
  int SIZEMEM=sizeof(double);
  //fft peakshape
  fftw_complex *P=nullptr,*C=nullptr;
  fftw_plan p1,pinv;
#else
  float *ivals,*jvals; //1D array of x and y vals
  float sig1,sig2,sig3;
  float lor1, lor2;
 
  float *DI=nullptr,*DB=nullptr,*DS=nullptr,*DBA=nullptr,*DBR=nullptr; //will hold x,y,I,Band S
  int SIZEMEM=sizeof(float);
  //fft peakshape
  fftwf_complex *P=nullptr,*C=nullptr;
  fftwf_plan p1,pinv;

#endif

  vector<peakEntry> peakList;
  double maxInt;
  float voigt1; //0 gaus 1 lorentz
  float voigt2; //0 gaus 1 lorentz
  int peaks;
  int symmode;
  int pipey=0;
  double noiseVal=0;

  bool BORE=false; //if we have a restricted peaklist.
  
  double DB_sum;
  int squash_window_i,squash_window_j;


  

  void LocalMax(double &hm,double &cm,int j);

//reads in human readable data file. outdated.
  void Read();

  void MapBlur();
  void UnMapBlur();
  void ReadPipe(pipe &pipefile);
  void InitBlurPeak();
  void InitBlur();

  double ApplyIter();
  double ApplyIter2();

#ifdef DOUBLE2D
  int DoIndex(double ref,double *vals,int ii);
#else
  int DoIndex(double ref,float *vals,int ii);
  #endif

 //get integer index of closest point in ppm
  void SetIndex();


  void MakeSquare();
  void SetBlur(double *ref);


  void FitPeaks2D(double radius1,
		  double radius2,
		  const std::string &paramOut = "fitparams.out",
		  const std::string &gnuplotOut = "fitted.out",
		  int maxIter = 50,
		  double threshold = 0.0,
                  bool useReferencePeakList = false,
                  bool restrictLMWidths = true);

    void BuildPeakListFromDB(double threshold,
                             std::vector<FitPeak2DLocal> &peaks);

    void BuildPeakListFromReference(std::vector<FitPeak2DLocal> &peaks);

    bool WriteFudaFitOutputs(const std::string &fitDir,
                             double radius1, double radius2,
                             double obs1MHz, double obs2MHz);


    // Pseudo-3D (3P) support: the shape is fitted once from the summed 2D
    // projection, then independent amplitudes are solved for every real-axis
    // slice using that fixed shape.  Data are flattened slice-major.
    bool FitPseudo3DIntensities(const std::vector<double>& stack, int nslices,
                                double radius1, double radius2,
                                std::vector<std::vector<double> >& intensities,
                                std::vector<std::vector<double> >& intensityEsd,
                                const std::string& fitDir,
                                const std::vector<double>& zvals,
                                const std::string& zlabel,
                                double obs1MHz, double obs2MHz);

    std::vector<FitGroup2DLocal>
    BuildGroups(const std::vector<FitPeak2DLocal> &peaks,
                double radius1,
                double radius2);

    double FitOneGroup2D(const FitGroup2DLocal &group,
                         std::vector<FitPeak2DLocal> &peaks,
                         int maxIter,
                         std::vector<double> &model,
                         std::vector<double> &raw,double radIppm,double radJppm);


  std::vector<FitPeak2DLocal> fittedPeaks2D;

  void RebuildDSFromFit(double nsig = 6.0);
  
   
  //function to return a gaussian at a specific point
//  double Gaus(double x,double x0,double sig){
//    return exp(-pow(x-x0,2.)/(2.*sig*sig));}
//  //function to return a gaussian at a specific point
//  double Lorentz(double x,double x0,double sig){
//    return pow(sig/2,2)/(pow(x-x0,2)+pow(sig/2,2));}
//  //function to return a 2D peak value at a given point
//  double Peak(double x1,double x2,double sig1)
//  {
//    switch(peaky){
//    case 0:
//      return Gaus(x1,x2,sig1);
//      break;
//    case 1:
//      return Lorentz(x1,x2,sig1);
//      break;}
//  }

  void SetPeak(const slice2D &inst);
  void GetPeak();
  void CalcSpec();
  void BuildSparseDB(double cutoff);

  struct SparsePt2D
  {
    int i, j, ii;
    double val;
  };
  std::vector<double> pki, pkj;
  std::vector<SparsePt2D> sparseDB;

  void PrintSpec();

  void PrintSpecPure();

  double GetChi2();

  void Squash();

  void Cull(float frac = 1.0);

  void CullCentral(double Hlim,double Clim);

  int CountElements();
  int correlate(FILE* out_pt,string infile);

  ~slice2D()
  {

    delete [] ivals;
    delete [] jvals;
    delete [] DI;
    delete [] DS;
    delete [] DB;

    if(BORE)
      {
	delete [] DBA;
	delete [] DBR;
      }

    if(P)
      delete [] P;
    if(C)
      delete [] C;

#ifdef DOUBLE2D
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
