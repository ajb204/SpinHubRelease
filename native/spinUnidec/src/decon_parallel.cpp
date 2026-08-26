#include "deconMain.hpp"
#include "omp.h"

/*************************************************************/
/* Bayesian Deconvolver for NMR data                          */
/*************************************************************/
// Parallel entry point. The thread API must match the compile-
// time precision used by the dimension being executed.

namespace {

int threads_from_input(const decon &dec)
{
  return (dec.parsed == 1 && dec.ncpus > 0) ? dec.ncpus : 1;
}

void init_fftw_threads_1d(int nthreads)
{
#ifdef DOUBLE1D
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#else
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nthreads);
#endif
}

void init_fftw_threads_2d(int nthreads)
{
#ifdef DOUBLE2D
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#else
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nthreads);
#endif
}

void init_fftw_threads_3d(int nthreads)
{
#ifdef DOUBLE3D
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#else
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nthreads);
#endif
}

void init_fftw_threads_4d(int nthreads)
{
#ifdef DOUBLE4D
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#else
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nthreads);
#endif
}

void cleanup_fftw_threads_1d()
{
#ifdef DOUBLE1D
  fftw_cleanup_threads();
#else
  fftwf_cleanup_threads();
#endif
}

void cleanup_fftw_threads_2d()
{
#ifdef DOUBLE2D
  fftw_cleanup_threads();
#else
  fftwf_cleanup_threads();
#endif
}

void cleanup_fftw_threads_3d()
{
#ifdef DOUBLE3D
  fftw_cleanup_threads();
#else
  fftwf_cleanup_threads();
#endif
}

void cleanup_fftw_threads_4d()
{
#ifdef DOUBLE4D
  fftw_cleanup_threads();
#else
  fftwf_cleanup_threads();
#endif
}

} // namespace

int main(int argc,char *argv[])
{
  decon dec;
  dec.splash(argc,argv,true);

  const int nthreads = threads_from_input(dec);
  omp_set_num_threads(nthreads);

  switch(dec.dim){
  case 1: //1D unidec
    //init_fftw_threads_1d(nthreads);
    dec.Protocol1D(argc,argv);
    //cleanup_fftw_threads_1d();
    break;
  case 2: //2D unidec / physical pseudo2D restrained FIT
    //init_fftw_threads_2d(nthreads);
    if (dec.pseudo2DFit) dec.Protocol2PFit(argc,argv);
    else dec.Protocol2D(argc,argv);
    //cleanup_fftw_threads_2d();
    break;
  case 3: //3D using 1D slices
    //init_fftw_threads_3d(nthreads);
    if (dec.pseudo3D) dec.Protocol3P(argc,argv);
    else dec.Protocol3D(argc,argv);
    //cleanup_fftw_threads_3d();
    break;
  case 4: //4D via 2D slices
    //init_fftw_threads_4d(nthreads);
    dec.Protocol4D(argc,argv);
    //cleanup_fftw_threads_4d();
    break;
  default:
    cout << "Error: unsupported dimension " << dec.dim << endl;
    return 100;
  }

  cout << "exiting cleanly." << endl;
  return 0;
}
