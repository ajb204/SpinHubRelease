/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE_HPP
#define SLICE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
//#include <string>
//#include <vector>

//#include <cmath>


using namespace std;



class leakEntry
{
 public:

  int k;
  int l;
  double val;
};

class peakEntry
{
 public:
  double x,y,z;          //store peak positions
  double xOld,yOld,zOld; //used for local optimisation
  string name;           //name of peak
  //int index=0;  //index with respect to spec

  int indexI=0; //index with respect to spec
  int indexJ=0; //index with respect to spec
  int indexK=0; //index with respect to spec
  int indexL=0; //index with respect to spec
};

inline float Gaus(float x,float x0,float sig){
  return exp(-pow(x-x0,2.)/(2.*sig*sig));}
//function to return a gaussian at a specific point
inline float Lorentz(float x, float x0,float sig){
   if (sig ==0.0) {
      return 0;
    }
  return pow(sig/2,2)/(pow(x-x0,2)+pow(sig/2,2));}
//function to return a 2D peak value at a given point
inline float Peak(float x1,float x2,float sig1, float r1, float peaky)
{
  float voigt = (1-peaky)*Gaus(x1,x2,sig1/2.355)+(peaky)*Lorentz(x1,x2, r1);
  return voigt;
}
inline float PeakPV(float x1,float x2,float sig1, float r1, float peaky)
{
  float voigt = (1-peaky)*Gaus(x1,x2,sig1/2.355)+(peaky)*Lorentz(x1,x2, r1);
  return voigt;
}

inline float GaussPV(float xx,float x0,float Gamma)
{
  float sigma=Gamma/(2*sqrt(2*log(2)));
  return exp(-pow(xx-x0,2)/(2*sigma*sigma))*1./(sqrt(2*3.14159)*sigma);
}
inline float LorentzPV(float xx,float x0,float Gamma)
{
  return (Gamma/2.)/(pow(xx-x0,2.)+pow(Gamma/2.,2.) )/(3.14159);
}

//https://docs.mantidproject.org/nightly/fitting/fitfunctions/PseudoVoigt.html
inline float PseudoVoigt(float xx,float x0,float Gamma,float nu)
{
  float yvals=nu*GaussPV(xx,x0,Gamma)+(1-nu)*(LorentzPV(xx,x0,Gamma));
  float ymax=nu*GaussPV(x0,x0,Gamma)+(1-nu)*(LorentzPV(x0,x0,Gamma));
  return yvals/ymax;
}

class raw1D
{
 public:
  double x,y;
};

class raw2D
{
 public:
  double x,y,z;
};

class raw3D
{
 public:
  double x,y,z,a;
};

class raw4D
{
 public:
  double x,y,z,a,b;
};



static inline double PeakFraction1D(double dx,
                                    double sig,
                                    double lor,
                                    double voigt)
{
  return static_cast<double>(Peak(static_cast<float>(dx),
                                  0.0f,
                                  static_cast<float>(sig),
                                  static_cast<float>(lor),
                                  static_cast<float>(voigt)));
}


// Dawson's integral F(x) = exp(-x^2) integral_0^x exp(t^2) dt.
// This compact Numerical-Recipes style approximation avoids an FFT Hilbert
// transform inside the pseudo2D fitter while retaining high accuracy over the
// full line shape.  The function is odd by construction.
static inline double Dawson1D(double x)
{
  const int NMAX=6;
  const double H=0.4;
  static bool init=false;
  static double c[NMAX];
  if(!init) {
    for(int i=0;i<NMAX;++i) {
      const double a=(2.0*i+1.0)*H;
      c[i]=std::exp(-a*a);
    }
    init=true;
  }
  const int n0=2*static_cast<int>(0.5*x/H + (x>=0.0 ? 0.5 : -0.5));
  const double xp=x-n0*H;
  double e1=std::exp(2.0*xp*H), e2=e1*e1, d1=n0+1.0, d2=d1-2.0;
  double sum=0.0;
  for(int i=0;i<NMAX;++i, d1+=2.0, d2-=2.0) {
    sum += c[i]*(e1/d1 + 1.0/(d2*e1));
    e1 *= e2;
  }
  return 0.56418958354775628695 * std::exp(-xp*xp) * sum;
}

// Hilbert/quadrature partner of the existing unit-height pseudo-Voigt.
// Sign convention is chosen so positive phase gives positive dispersion on
// the +dx side.  At phase==0 PeakFraction1DPhased is bit-for-bit the old
// absorptive model apart from the final double precision multiply/add.
static inline double PeakDispersion1D(double dx,
                                      double sig,
                                      double lor,
                                      double voigt)
{
  const double safeSig=std::max(std::fabs(sig),1e-15);
  const double safeLor=std::max(std::fabs(lor),1e-15);
  const double sigma=safeSig/2.355;
  const double z=dx/(std::sqrt(2.0)*sigma);
  const double dg=(2.0/std::sqrt(3.14159265358979323846))*Dawson1D(z);
  const double gamma=0.5*safeLor;
  const double dl=(gamma*dx)/(dx*dx+gamma*gamma);
  return (1.0-voigt)*dg + voigt*dl;
}

static inline double PeakFraction1DPhased(double dx,
                                          double sig,
                                          double lor,
                                          double voigt,
                                          double phase)
{
  const double a=PeakFraction1D(dx,sig,lor,voigt);
  if(phase==0.0) return a;
  return a*std::cos(phase) + PeakDispersion1D(dx,sig,lor,voigt)*std::sin(phase);
}

static double SolveRadiusForThreshold1D(double target,
                                        double sig,
                                        double lor,
                                        double voigt)
{
  if (!(target > 0.0) || target >= 1.0)
    return 0.0;

  const double safeSig = std::max(sig, 1e-12);
  const double safeLor = std::max(lor, 1e-12);

  const double sigma = safeSig / 2.355;
  const double gausRadius =
      sigma * std::sqrt(std::max(0.0, -2.0 * std::log(target)));

  const double halfLor = 0.5 * safeLor;
  const double lorRadius =
      halfLor * std::sqrt(std::max(0.0, 1.0 / target - 1.0));

  double hi = std::max(gausRadius, lorRadius);
  if (!(hi > 0.0) || !std::isfinite(hi))
    hi = std::max(safeSig, safeLor);
  if (!(hi > 0.0) || !std::isfinite(hi))
    hi = 1.0;

  double lo = 0.0;
  double fhi = PeakFraction1D(hi, safeSig, safeLor, voigt) - target;

  int expand = 0;
  while (fhi > 0.0 && expand < 48)
  {
    hi *= 2.0;
    fhi = PeakFraction1D(hi, safeSig, safeLor, voigt) - target;
    ++expand;
  }

  if (fhi > 0.0)
    return hi;

  for (int it = 0; it < 80; ++it)
  {
    const double mid = 0.5 * (lo + hi);
    const double fm = PeakFraction1D(mid, safeSig, safeLor, voigt) - target;
    if (fm > 0.0)
      lo = mid;
    else
      hi = mid;
  }

  return 0.5 * (lo + hi);
}


static bool SolveLinearSystem(std::vector<double> A,
                              std::vector<double> b,
                              std::vector<double>& x,
                              int n)
{
  const double eps = 1e-14;
  x.assign(n, 0.0);

  for (int k = 0; k < n; ++k)
  {
    int piv = k;
    double best = std::fabs(A[k * n + k]);
    for (int r = k + 1; r < n; ++r)
    {
      const double v = std::fabs(A[r * n + k]);
      if (v > best)
      {
        best = v;
        piv = r;
      }
    }

    if (best < eps)
      return false;

    if (piv != k)
    {
      for (int c = k; c < n; ++c)
        std::swap(A[k * n + c], A[piv * n + c]);
      std::swap(b[k], b[piv]);
    }

    const double diag = A[k * n + k];
    for (int c = k; c < n; ++c)
      A[k * n + c] /= diag;
    b[k] /= diag;

    for (int r = k + 1; r < n; ++r)
    {
      const double f = A[r * n + k];
      if (std::fabs(f) < eps)
        continue;
      for (int c = k; c < n; ++c)
        A[r * n + c] -= f * A[k * n + c];
      b[r] -= f * b[k];
    }
  }

  for (int k = n - 1; k >= 0; --k)
  {
    double v = b[k];
    for (int c = k + 1; c < n; ++c)
      v -= A[k * n + c] * x[c];
    x[k] = v;
  }

  return true;
}


struct Peak1DParts
{
  double value;
  double d_x0;
  double d_sig;
  double d_lor;
  double d_voigt;
};


static inline Peak1DParts EvaluatePeak1D(double x, double x0, double sig, double lor, double voigt)
{
  const double k = 2.355;
  const double safeSig = std::max(std::fabs(sig), 1e-12);
  const double s = safeSig / k;
  const double dx = x - x0;
  const double dx2 = dx * dx;
  const double s2 = s * s;

  const double g = std::exp(-dx2 / (2.0 * s2));

  const double a = 0.5 * std::max(std::fabs(lor), 0.0);
  const double a2 = a * a;
  const double denom = dx2 + a2;
  double l = 0.0;
  double dldx0 = 0.0;
  double dldr = 0.0;

  if (a > 0.0)
  {
    l = a2 / denom;
    const double denom2 = denom * denom;
    dldx0 = 2.0 * dx * a2 / denom2;
    dldr = a * dx2 / denom2;
  }

  Peak1DParts out;
  out.value = (1.0 - voigt) * g + voigt * l;
  out.d_x0 = (1.0 - voigt) * (g * dx / s2) + voigt * dldx0;
  out.d_sig = (1.0 - voigt) * (g * dx2 / (s2 * safeSig));
  out.d_lor = voigt * dldr;
  out.d_voigt = l - g;
  return out;
}


static inline double clampd(double v, double lo, double hi)
{
  return std::max(lo, std::min(v, hi));
}



#endif
