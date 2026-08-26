#ifndef FFTW3_COMPAT_H
#define FFTW3_COMPAT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double fftw_complex[2];
typedef float fftwf_complex[2];
typedef struct fftw_plan_s *fftw_plan;
typedef struct fftwf_plan_s *fftwf_plan;

enum { FFTW_MEASURE = 0, FFTW_ESTIMATE = 64, FFTW_PATIENT = 32, FFTW_DESTROY_INPUT = 1 << 27 };

void fftw_execute(const fftw_plan p);
void fftwf_execute(const fftwf_plan p);
void fftw_forget_wisdom(void);
int fftw_export_to_filename(const char *filename);

fftw_complex *fftw_alloc_complex(size_t n);
fftwf_complex *fftwf_alloc_complex(size_t n);

fftw_plan fftw_plan_dft_r2c_1d(int n0, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r_1d(int n0, fftw_complex *in, double *out, unsigned flags);
fftw_plan fftw_plan_dft_r2c_2d(int n0, int n1, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r_2d(int n0, int n1, fftw_complex *in, double *out, unsigned flags);
fftw_plan fftw_plan_dft_r2c_3d(int n0, int n1, int n2, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r_3d(int n0, int n1, int n2, fftw_complex *in, double *out, unsigned flags);
fftw_plan fftw_plan_dft_r2c(int rank, const int *n, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r(int rank, const int *n, fftw_complex *in, double *out, unsigned flags);

int fftwf_init_threads(void);
void fftwf_plan_with_nthreads(int nthreads);
fftwf_plan fftwf_plan_dft_r2c_1d(int n0, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r_1d(int n0, fftwf_complex *in, float *out, unsigned flags);
fftwf_plan fftwf_plan_dft_r2c_2d(int n0, int n1, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r_2d(int n0, int n1, fftwf_complex *in, float *out, unsigned flags);
fftwf_plan fftwf_plan_dft_r2c_3d(int n0, int n1, int n2, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r_3d(int n0, int n1, int n2, fftwf_complex *in, float *out, unsigned flags);
fftwf_plan fftwf_plan_dft_r2c(int rank, const int *n, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r(int rank, const int *n, fftwf_complex *in, float *out, unsigned flags);

  //added 29th July 2026.
int fftw_init_threads(void);
void fftw_plan_with_nthreads(int nthreads);
void fftw_cleanup_threads(void);

int fftwf_init_threads(void);
void fftwf_plan_with_nthreads(int nthreads);
void fftwf_cleanup_threads(void);


/* Plan management */
void fftwf_destroy_plan(fftwf_plan plan);
void fftw_destroy_plan(fftw_plan plan);

#ifdef __cplusplus
}
#endif

#endif
