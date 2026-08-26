#ifndef IMSL_H
#define IMSL_H

#include "imslerr.h"

typedef struct {
    float 	re;
    float 	im;
} f_complex;

typedef struct {
    int   	domain_dim;
    int   	target_dim;
    int   	*order;
    int   	*num_coef;
    int   	*num_breakpoints;
    float 	**breakpoints;
    float 	**coef;
} Imsl_f_ppoly;

typedef struct {
    int   	domain_dim;
    int   	target_dim;
    int   	*order;
    int   	*num_coef;
    int   	*num_knots;
    float 	**knots;
    float 	**coef;
} Imsl_f_spline;

typedef struct {
    double	re;
    double 	im;
} d_complex;

typedef struct {
    int   	domain_dim;
    int   	target_dim;
    int   	*order;
    int   	*num_coef;
    int   	*num_breakpoints;
    double 	**breakpoints;
    double 	**coef;
} Imsl_d_ppoly;

typedef struct {
    int   	domain_dim;
    int   	target_dim;
    int   	*order;
    int   	*num_coef;
    int   	*num_knots;
    double 	**knots;
    double 	**coef;
} Imsl_d_spline;

typedef enum {
    IMSL_NOTE		    = 1,
    IMSL_ALERT		    = 2,
    IMSL_WARNING	    = 3,
    IMSL_FATAL		    = 4,
    IMSL_TERMINAL	    = 5,
    IMSL_WARNING_IMMEDIATE  = 6,
    IMSL_FATAL_IMMEDIATE    = 7
} Imsl_error;

typedef enum {
    IMSL_ODE_INITIALIZE	    = 1,
    IMSL_ODE_CHANGE	    = 2,
    IMSL_ODE_RESET	    = 3
} Imsl_ode;
    				/* enums for imsl_page */
typedef enum {
    IMSL_SET_PAGE_WIDTH     =  -1,
    IMSL_GET_PAGE_WIDTH     =   1,
    IMSL_SET_PAGE_LENGTH    =  -2,
    IMSL_GET_PAGE_LENGTH    =   2
} Imsl_page_options;
    				/* enums for imsl_write_options */
typedef enum {
    IMSL_SET_DEFAULTS       =   0,
    IMSL_SET_CENTERING      =  -1,
    IMSL_GET_CENTERING      =   1,
    IMSL_SET_ROW_WRAP       =  -2,
    IMSL_GET_ROW_WRAP       =   2,
    IMSL_SET_PAGING         =  -3,
    IMSL_GET_PAGING         =   3,
    IMSL_SET_NAN_CHAR       =  -4,
    IMSL_GET_NAN_CHAR       =   4,
    IMSL_SET_TITLE_PAGE     =  -5,
    IMSL_GET_TITLE_PAGE     =   5,
    IMSL_SET_FORMAT         =  -6,
    IMSL_GET_FORMAT         =   6
} Imsl_write_options;
				/* enums for quadrature routines */
typedef enum {
    IMSL_ALG		    = 1,
    IMSL_ALG_LEFT_LOG	    = 2,
    IMSL_ALG_RIGHT_LOG	    = 3,
    IMSL_ALG_LOG	    = 4,
    IMSL_INF_BOUND	    = 5,
    IMSL_BOUND_INF	    = 6,
    IMSL_INF_INF            = 7,
    IMSL_COS		    = 8,
    IMSL_SIN		    = 9
} Imsl_quad;

#ifdef USE_IMSL_C

#if defined(ANSI) || __STDC__
#define IMSL_PROTO(P,Q)    P Q
#else
#define IMSL_PROTO(P,Q)    P()
#endif

typedef void IMSL_PROTO((*Imsl_error_print_proc),
			(Imsl_error,long,char*,char*));

	/* Chapter 1 --- Linear System */

float	      * IMSL_PROTO(imsl_f_lin_sol_gen,
		      (int n, float* a, float *b, ...));
double	      * IMSL_PROTO(imsl_d_lin_sol_gen,
		      (int n, double* a, double *b, ...));
float         * IMSL_PROTO(imsl_f_lin_sol_posdef,
		      (int, float*, float*, ...));
double        * IMSL_PROTO(imsl_d_lin_sol_posdef,
		      (int, double*, double*, ...));
f_complex     * IMSL_PROTO(imsl_c_lin_sol_posdef,
		      (int, f_complex*, f_complex*, ...));
d_complex     * IMSL_PROTO(imsl_z_lin_sol_posdef,
		      (int, d_complex*, d_complex*, ...));
float         * IMSL_PROTO(imsl_f_lin_cgsol_posdef,
		      (int, void (*fcn)(float*,float*),float*, ...));
double        * IMSL_PROTO(imsl_d_lin_cgsol_posdef,
		      (int, void (*fcn)(double*,double*),double*, ...));
float	      * IMSL_PROTO(imsl_f_lin_sol_nonnegdef,
		      (int n, float *a, float *b, ...));
double        * IMSL_PROTO(imsl_d_lin_sol_nonnegdef,
		      (int n, double *a, double *b, ...));
f_complex     * IMSL_PROTO(imsl_c_lin_sol_gen,
		      (int, f_complex*, f_complex*, ...));
d_complex     * IMSL_PROTO(imsl_z_lin_sol_gen,
		      (int, d_complex*, d_complex*, ...));
float         * IMSL_PROTO(imsl_f_lin_least_squares_gen,
		      (int, int, float*, float*, ...));
double        * IMSL_PROTO(imsl_d_lin_least_squares_gen,
		      (int, int, double*, double*, ...));
float         * IMSL_PROTO(imsl_f_mat_mul_rect,
		      (char*, ...));
double        * IMSL_PROTO(imsl_d_mat_mul_rect,
		      (char*, ...));
f_complex     * IMSL_PROTO(imsl_c_mat_mul_rect,
		      (char*, ...));
d_complex     * IMSL_PROTO(imsl_z_mat_mul_rect,
		      (char*, ...));
float	      * IMSL_PROTO(imsl_f_lin_svd_gen,
		      (int, int, float*, ...));
double        * IMSL_PROTO(imsl_d_lin_svd_gen,
		      (int, int, double*, ...));
f_complex     * IMSL_PROTO(imsl_c_lin_svd_gen,
		      (int, int, f_complex*, ...));
d_complex     * IMSL_PROTO(imsl_z_lin_svd_gen,
		      (int, int, d_complex*, ...));

	/* Chapter 2 --- Eigenvalues */

float	      * IMSL_PROTO(imsl_f_eig_sym,
		      (int, float *, ...));
double	      * IMSL_PROTO(imsl_d_eig_sym,
		      (int, double *, ...));
f_complex     * IMSL_PROTO(imsl_f_eig_gen,
		      (int, float *, ...));
d_complex     * IMSL_PROTO(imsl_d_eig_gen,
		      (int, double *, ...));
f_complex     * IMSL_PROTO(imsl_c_eig_gen,
		      (int, f_complex *, ...));
d_complex     * IMSL_PROTO(imsl_z_eig_gen,
		      (int, d_complex *, ...));
float         * IMSL_PROTO(imsl_c_eig_herm,
		      (int, f_complex *, ...));
double        * IMSL_PROTO(imsl_z_eig_herm,
		      (int, d_complex *, ...));
float	      * IMSL_PROTO(imsl_f_geneig_sym_posdef,
		      (int, float *, float *, ...));
double	      * IMSL_PROTO(imsl_d_geneig_sym_posdef,
		      (int, double *, double *, ...));
float	      * IMSL_PROTO(imsl_f_eig_symgen,
		      (int, float *, float *, ...));
double	      * IMSL_PROTO(imsl_d_eig_symgen,
		      (int, double *, double *, ...));

	/* Chapter 3 --- Interpolation and Approximation */

Imsl_f_ppoly  * IMSL_PROTO(imsl_f_ppoly_create,
		      (int, int, int*, int*, ...));
Imsl_d_ppoly  * IMSL_PROTO(imsl_d_ppoly_create,
		      (int, int, int*, int*, ...));
Imsl_f_ppoly  * IMSL_PROTO(imsl_f_cub_spline_interp,
		      (int, float[], float[]));
Imsl_d_ppoly  * IMSL_PROTO(imsl_d_cub_spline_interp,
		      (int, double[], double[]));
Imsl_f_ppoly  * IMSL_PROTO(imsl_f_cub_spline_interp_e_cnd,
		      (int, float[], float[], ...));
Imsl_d_ppoly  * IMSL_PROTO(imsl_d_cub_spline_interp_e_cnd,
		      (int, double[], double[], ...));
Imsl_f_ppoly  * IMSL_PROTO(imsl_f_cub_spline_interp_shape,
		      (int, float[], float[],...));
Imsl_d_ppoly  * IMSL_PROTO(imsl_d_cub_spline_interp_shape,
		      (int, double[], double[],...));
float         	IMSL_PROTO(imsl_f_cub_spline_value,
		      (float, Imsl_f_ppoly*,...));
double        	IMSL_PROTO(imsl_d_cub_spline_value,
		      (double, Imsl_d_ppoly*,...));
float         	IMSL_PROTO(imsl_f_cub_spline_integral,
		      (float, float, Imsl_f_ppoly*));
double        	IMSL_PROTO(imsl_d_cub_spline_integral,
		      (double, double, Imsl_d_ppoly*));
Imsl_f_spline * IMSL_PROTO(imsl_f_spline_create,
		      (int, int, int*, int*, ...));
Imsl_d_spline * IMSL_PROTO(imsl_d_spline_create,
		      (int,int,int*, int*, ...));
Imsl_f_spline * IMSL_PROTO(imsl_f_spline_interp,
		      (int,float[],float[],...));
Imsl_d_spline * IMSL_PROTO(imsl_d_spline_interp,
		      (int, double[],double[],...));
float         * IMSL_PROTO(imsl_f_spline_knots,
		      (int, float[], ...));
double        * IMSL_PROTO(imsl_d_spline_knots,
		      (int, double[], ...));
float         	IMSL_PROTO(imsl_f_spline_value,
		      (float, Imsl_f_spline*,...));
double         	IMSL_PROTO(imsl_d_spline_value,
		      (double, Imsl_d_spline*,...));
float          	IMSL_PROTO(imsl_f_spline_integral,
		      (float, float, Imsl_f_spline*));
double         	IMSL_PROTO(imsl_d_spline_integral,
		      (double, double, Imsl_d_spline*));
Imsl_f_spline * IMSL_PROTO(imsl_f_spline_2d_interp,
		      (int,float[], int, float[], float[],...));
Imsl_d_spline * IMSL_PROTO(imsl_d_spline_2d_interp,
		      (int, double[],int, double[], double[],...));
float          	IMSL_PROTO(imsl_f_spline_2d_value,
		      (float,float, Imsl_f_spline*,...));
double         	IMSL_PROTO(imsl_d_spline_2d_value,
		      (double,double, Imsl_d_spline*,...));
float          	IMSL_PROTO(imsl_f_spline_2d_integral,
		      (float, float, float, float, Imsl_f_spline*));
double         	IMSL_PROTO(imsl_d_spline_2d_integral,
		      (double, double,double,double, Imsl_d_spline*));
Imsl_f_spline * IMSL_PROTO(imsl_f_spline_2d_least_squares,
		      (int, float[], int, float[], float[], int, int,...));
Imsl_d_spline * IMSL_PROTO(imsl_d_spline_2d_least_squares,
		      (int, double[], int, double[], double[], int, int,...));
Imsl_f_spline * IMSL_PROTO(imsl_f_spline_least_squares,
		      (int, float[], float[], int,...));
Imsl_d_spline * IMSL_PROTO(imsl_d_spline_least_squares,
		      (int, double[], double[], int,...));
float         * IMSL_PROTO(imsl_f_poly_regression,
		      (int, float[], float[], int, ...));
double        * IMSL_PROTO(imsl_d_poly_regression,
		      (int, double[], double[], int, ...));
void           	IMSL_PROTO(imsl_f_spline_print,
		      (Imsl_f_spline*));
void           	IMSL_PROTO(imsl_d_spline_print,
		      (Imsl_d_spline*));
void           	IMSL_PROTO(imsl_f_ppoly_print,
		      (Imsl_f_ppoly*));
void           	IMSL_PROTO(imsl_d_ppoly_print,
		      (Imsl_d_ppoly*));
float         * IMSL_PROTO(imsl_f_user_fcn_least_squares,
		      (float (*fcn) (int, float), int, int, float[], 
		       float[], ...));
double        * IMSL_PROTO(imsl_d_user_fcn_least_squares,
		      (double (*fcn) (int, double), int, int, double[], 
		       double[], ...));
Imsl_f_ppoly  * IMSL_PROTO(imsl_f_cub_spline_smooth,
		      (int, float[], float[], ...));
Imsl_d_ppoly  * IMSL_PROTO(imsl_d_cub_spline_smooth,
		      (int, double[], double[], ...));
float         * IMSL_PROTO(imsl_f_scattered_2d_interp,
		      (int, float[],  float[],  int, int, float[],
		       float[], ...));
double        * IMSL_PROTO(imsl_d_scattered_2d_interp,
		      (int, double*, double[], int, int, double[], 
		       double[], ...));

	/* Chapter 4 --- Intergration and Differentiation */

float        	IMSL_PROTO(imsl_f_int_fcn,
		      (float (*fcn) (float), float, float, ...));
double       	IMSL_PROTO(imsl_d_int_fcn,
		      (double (*fcn) (double), double, double, ...));
float        	IMSL_PROTO(imsl_f_int_fcn_sing_pts,
		      (float (*fcn) (float), float, float, int, float*, ...));
double       	IMSL_PROTO(imsl_d_int_fcn_sing_pts,
		      (double (*fcn) (double), double, double, int, 
		       double*, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_sing,
		      (float (*fcn) (float), float, float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_sing,
		      (double (*fcn) (double), double, double, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_inf,
		      (float (*fcn) (float), float, Imsl_quad, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_inf,
		      (double (*fcn) (double), double, Imsl_quad, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_trig,
		      (float (*fcn) (float), float, float, Imsl_quad, 
		       float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_trig,
		      (double (*fcn) (double), double, double,
		       Imsl_quad, double, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_alg_log,
		      (float (*fcn) (float), float, float,
		       Imsl_quad, float, float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_alg_log,
		      (double (*fcn) (double), double,
		       double, Imsl_quad, double, double, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_fourier,
		      (float (*fcn) (float), float, Imsl_quad, float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_fourier,
		      (double (*fcn) (double), double,
		       Imsl_quad, double, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_smooth,
		      (float (*fcn) (float), float, float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_smooth,
		      (double (*fcn) (double), double, double, ...));
float	     	IMSL_PROTO(imsl_f_int_fcn_cauchy,
		      (float (*fcn) (float), float, float, float, ...));
double	     	IMSL_PROTO(imsl_d_int_fcn_cauchy,
		      (double (*fcn) (double), double, double, double, ...));
float        	IMSL_PROTO(imsl_f_int_fcn_hyper_rect,
		      (float (*fcn) (int, float*), int, float*, float*, ...));
double       	IMSL_PROTO(imsl_d_int_fcn_hyper_rect,
		      (double (*fcn) (int, double*), int, double*,
		       double*, ...));
float        	IMSL_PROTO(imsl_f_int_fcn_2d,
		      (float (*fcn1) (float,float), float, float,
		       float (*fcn2) (float),float (*fcn3) (float), ...));
double       	IMSL_PROTO(imsl_d_int_fcn_2d,
		      (double (*fcn1) (double,double), double, double,
		       double (*fcn2) (double),double (*fcn3) (double), ...));
void         	IMSL_PROTO(imsl_f_gauss_quad_rule,
		      (int, float[],  float[], ...));
void         	IMSL_PROTO(imsl_d_gauss_quad_rule,
		      (int, double[], double[], ...));

	/* Chapter 5 --- Differential Equations */

void	     	IMSL_PROTO(imsl_f_ode_runge_kutta,
		      (int, float*, float, float*, char *,
		       void (*fcn)(int,float,float*,float*)));
void	     	IMSL_PROTO(imsl_d_ode_runge_kutta,
		      (int, double*, double, double*, char*,
		       void (*fcn)(int,double,double*,double*)));
void	     	IMSL_PROTO(imsl_f_ode_runge_kutta_mgr,
		      (int, char**, ...));
void	     	IMSL_PROTO(imsl_d_ode_runge_kutta_mgr,
		      (int, char**, ...));
void	     	IMSL_PROTO(imsl_f_ode_adams_gear,
		      (int,float*, float, float*, char*,
		       void (*fcn)(int,float,float*,float*)));
void	     	IMSL_PROTO(imsl_d_ode_adams_gear,
		      (int,double*,double,double*,char*,
		       void (*fcn)(int,double,double*,double*)));
void	     	IMSL_PROTO(imsl_f_ode_adams_gear_mgr,
		      (int, char**, ...));
void	     	IMSL_PROTO(imsl_d_ode_adams_gear_mgr,
		      (int, char**, ...));

	/* Chapter 6 --- Transforms */

float         * IMSL_PROTO(imsl_f_fft_real,
		      (int, float*, ...));
double        * IMSL_PROTO(imsl_d_fft_real,
		      (int, double*, ...));
float         * IMSL_PROTO(imsl_f_fft_real_init,
		      (int));
double        * IMSL_PROTO(imsl_d_fft_real_init,
		      (int));
float         * IMSL_PROTO(imsl_c_fft_complex_init,
		      (int));
double        * IMSL_PROTO(imsl_z_fft_complex_init,
		      (int));
f_complex     * IMSL_PROTO(imsl_c_fft_complex,
		      (int, f_complex*, ...));
d_complex     * IMSL_PROTO(imsl_z_fft_complex,
		      (int, d_complex*, ...));
f_complex     * IMSL_PROTO(imsl_c_fft_2d_complex,
		      (int, int, f_complex*, ...));
d_complex     * IMSL_PROTO(imsl_z_fft_2d_complex,
		      (int, int, d_complex*, ...));

	/* Chapter 7 --- Nonlinear Equations */

f_complex     * IMSL_PROTO(imsl_f_zeros_poly,
		      (int, float*, ...));
d_complex     * IMSL_PROTO(imsl_d_zeros_poly,
		      (int, double*, ...));
float         * IMSL_PROTO(imsl_f_zeros_fcn,
		      (float (*fcn)(float), ...));
double        * IMSL_PROTO(imsl_d_zeros_fcn,
		      (double (*fcn)(double), ...));
float         * IMSL_PROTO(imsl_f_zeros_sys_eqn,
		      (void (*fcn)(int, float[], float[]), int, ...));
double        * IMSL_PROTO(imsl_d_zeros_sys_eqn,
		      (void (*fcn)(int, double[], double[]), int, ...));

	/* Chapter 8 --- Optimization */

float        	IMSL_PROTO(imsl_f_min_uncon,
		      (float (*fcn)(float), float, float, ...));
double       	IMSL_PROTO(imsl_d_min_uncon,
		      (double (*fcn)(double), double, double, ...));
float        	IMSL_PROTO(imsl_f_min_uncon_deriv,
		      (float (*fcn) (float), float (*grad) (float), float,
		       float, ...));
double       	IMSL_PROTO(imsl_d_min_uncon_deriv,
		      (double (*fcn) (double), double (*grad) (double), double,
		       double, ...));
float         * IMSL_PROTO(imsl_f_min_uncon_multivar,
		      (float (*fcn)(int, float[]), int, ...));
double        * IMSL_PROTO(imsl_d_min_uncon_multivar,
		      (double (*fcn)(int, double[]), int, ...));
float         * IMSL_PROTO(imsl_f_nonlin_least_squares,
		      (void (*fcn)(int, int, float[], float[]),
		       int, int, ...));
double        * IMSL_PROTO(imsl_d_nonlin_least_squares,
		      (void (*fcn)(int, int, double[], double[]), 
		       int, int, ...));
float         * IMSL_PROTO(imsl_f_min_con_nonlin,
		      (void (*fcn)(int, int, int, float[], int[], float*, 
				   float[]), 
		       int, int, int, int, float[], float[], ...));
double        * IMSL_PROTO(imsl_d_min_con_nonlin,
		      (void (*fcn)(int, int, int, double[],
				   int[], double*, double[]), 
		       int, int, int, int, double[], double[], ...));
float         * IMSL_PROTO(imsl_f_lin_prog,
		      (int, int, float*, float[], float[], ...));
double        * IMSL_PROTO(imsl_d_lin_prog,
		      (int, int, double*, double[], double[], ...));
float         * IMSL_PROTO(imsl_f_quadratic_prog,
		      (int, int, int, float*, float*, float*, float*, ...));
double        * IMSL_PROTO(imsl_d_quadratic_prog,
		      (int, int, int, double*, double*, double*, double*, ...));

	/* Chapter 9 --- Special Functions */

int	    	IMSL_PROTO(imsl_i_min,
		      (int, int));
float	    	IMSL_PROTO(imsl_f_min,
		      (float, float));
double	    	IMSL_PROTO(imsl_d_min,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_vmin,
		      (int, ...));
double	    	IMSL_PROTO(imsl_d_vmin,
		      (int, ...));
int	    	IMSL_PROTO(imsl_i_max,
		      (int, int));
float	    	IMSL_PROTO(imsl_f_max,
		      (float, float));
double	    	IMSL_PROTO(imsl_d_max,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_vmax,
		      (int, ...));
double	    	IMSL_PROTO(imsl_d_vmax,
		      (int, ...));
int	    	IMSL_PROTO(imsl_ii_power,
		      (int, int));
float	    	IMSL_PROTO(imsl_fi_power,
		      (float, int));
double	    	IMSL_PROTO(imsl_di_power,
		      (double, int));
float	    	IMSL_PROTO(imsl_ff_power,
		      (float, float));
double	    	IMSL_PROTO(imsl_dd_power,
		      (double, double));
f_complex   	IMSL_PROTO(imsl_ci_power,
		      (f_complex, int));
d_complex   	IMSL_PROTO(imsl_zi_power,
		      (d_complex, int));
f_complex   	IMSL_PROTO(imsl_cf_power,
		      (f_complex, float));
d_complex   	IMSL_PROTO(imsl_zd_power,
		      (d_complex, double));
f_complex   	IMSL_PROTO(imsl_cc_power,
		      (f_complex, f_complex));
d_complex   	IMSL_PROTO(imsl_zz_power,
		      (d_complex, d_complex));
float       	IMSL_PROTO(imsl_f_erf,
		      (float));
double      	IMSL_PROTO(imsl_d_erf,
		      (double));
float       	IMSL_PROTO(imsl_f_erfc,
		      (float));
double      	IMSL_PROTO(imsl_d_erfc,
		      (double));
float       	IMSL_PROTO(imsl_f_erf_inverse,
		      (float));
double      	IMSL_PROTO(imsl_d_erf_inverse,
		      (double));
float       	IMSL_PROTO(imsl_f_erfc_inverse,
		      (float));
double      	IMSL_PROTO(imsl_d_erfc_inverse,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_J0,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_J0,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_J1,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_J1,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_I0,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_I0,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_I1,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_I1,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_Y0,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_Y0,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_Y1,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_Y1,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_K0,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_K0,
		      (double));
float       	IMSL_PROTO(imsl_f_bessel_K1,
		      (float));
double      	IMSL_PROTO(imsl_d_bessel_K1,
		      (double));
float       	IMSL_PROTO(imsl_f_gamma,
		      (float));
double      	IMSL_PROTO(imsl_d_gamma,
		      (double));
float	    	IMSL_PROTO(imsl_f_log_gamma,
		      (float));
double	    	IMSL_PROTO(imsl_d_log_gamma,
		      (double));
float	    	IMSL_PROTO(imsl_f_gamma_incomplete,
		      (float, float));
double	    	IMSL_PROTO(imsl_d_gamma_incomplete,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_t_inverse_cdf,
		      (float, float));
double      	IMSL_PROTO(imsl_d_t_inverse_cdf,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_normal_inverse_cdf,
		      (float));
double	    	IMSL_PROTO(imsl_d_normal_inverse_cdf,
		      (double));
float	    	IMSL_PROTO(imsl_f_binomial_cdf,
		      (int, int, float));
double      	IMSL_PROTO(imsl_d_binomial_cdf,
		      (int, int, double));
float       	IMSL_PROTO(imsl_f_normal_cdf,
		      (float));
double      	IMSL_PROTO(imsl_d_normal_cdf,
		      (double));
float	    	IMSL_PROTO(imsl_f_chi_squared_cdf,
		      (float, float));
double      	IMSL_PROTO(imsl_d_chi_squared_cdf,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_chi_squared_inverse_cdf,
		      (float, float));
double      	IMSL_PROTO(imsl_d_chi_squared_inverse_cdf,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_F_cdf,
		      (float, float, float));
double      	IMSL_PROTO(imsl_d_F_cdf,
		      (double, double, double));
float	    	IMSL_PROTO(imsl_f_gamma_cdf,
		      (float, float));
double      	IMSL_PROTO(imsl_d_gamma_cdf,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_t_cdf,
		      (float, float));
double      	IMSL_PROTO(imsl_d_t_cdf,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_F_inverse_cdf,
		      (float, float, float));
double      	IMSL_PROTO(imsl_d_F_inverse_cdf,
		      (double, double, double));
float	    	IMSL_PROTO(imsl_f_hypergeometric_cdf,
		      (int, int, int, int));
double      	IMSL_PROTO(imsl_d_hypergeometric_cdf,
		      (int, int, int, int));
float       	IMSL_PROTO(imsl_f_poisson_cdf,
		      (int, float));
double      	IMSL_PROTO(imsl_d_poisson_cdf,
		      (int, double));

	/* beta */
float	    	IMSL_PROTO(imsl_f_beta,
		      (float, float));
double      	IMSL_PROTO(imsl_d_beta,
		      (double, double));
float	    	IMSL_PROTO(imsl_f_beta_incomplete,
		      (float, float, float));
double      	IMSL_PROTO(imsl_d_beta_incomplete,
		      (double, double, double));
float	    	IMSL_PROTO(imsl_f_log_beta,
		      (float, float));
double      	IMSL_PROTO(imsl_d_log_beta,
		      (double, double));

	/* Chapter 10 --- Printing Routines */

void        	IMSL_PROTO(imsl_i_write_matrix,
		      (char*, int, int, int[], ...));
void        	IMSL_PROTO(imsl_f_write_matrix,
		      (char*, int, int, float[], ...));
void        	IMSL_PROTO(imsl_c_write_matrix,
		      (char*, int, int, f_complex[], ...));
void        	IMSL_PROTO(imsl_d_write_matrix,
		      (char*, int, int, double[], ...));
void        	IMSL_PROTO(imsl_z_write_matrix,
		      (char*, int, int, d_complex[], ...));
void        	IMSL_PROTO(imsl_page,
		      (int, int*));
void        	IMSL_PROTO(imsl_write_options,
		      (int, int*));

	/* Chapter 11 --- Statistics */

float	      * IMSL_PROTO(imsl_f_regression,
		      (int, int, float*, float[], ...));
double	      * IMSL_PROTO(imsl_d_regression,
		      (int, int, double*, double[], ...));
float	      * IMSL_PROTO(imsl_f_ranks,
		      (int, float[], ...));
double	      * IMSL_PROTO(imsl_d_ranks,
		      (int, double[], ...));
float	      * IMSL_PROTO(imsl_f_simple_statistics,
		      (int, int, float*, ...));
double        * IMSL_PROTO(imsl_d_simple_statistics,
		      (int, int, double*, ...));
float	     	IMSL_PROTO(imsl_f_chi_squared_test,
		      (float (*fcn)(float), int, int, float*, ...));
double       	IMSL_PROTO(imsl_d_chi_squared_test,
		      (double (*fcn)(double), int, int, double*, ...));
float	      * IMSL_PROTO(imsl_f_covariances,
		      (int, int, float*, ...));
double        * IMSL_PROTO(imsl_d_covariances,
		      (int, int, double*, ...));
void         	IMSL_PROTO(imsl_random_seed_set,
		      (int));
int          	IMSL_PROTO(imsl_random_seed_get,
		      (void));
void         	IMSL_PROTO(imsl_random_option,
		      (int));
float	      * IMSL_PROTO(imsl_f_random_uniform,
		      (int, ...));
double	      * IMSL_PROTO(imsl_d_random_uniform,
		      (int, ...));
int	      * IMSL_PROTO(imsl_random_poisson,
		      (int, float, ...));
float	      * IMSL_PROTO(imsl_f_random_normal,
		      (int, ...));
double        * IMSL_PROTO(imsl_d_random_normal,
		      (int, ...));
float	      * IMSL_PROTO(imsl_f_random_gamma,
		      (int, float, ...));
double        * IMSL_PROTO(imsl_d_random_gamma,
		      (int, double, ...));
float	      * IMSL_PROTO(imsl_f_random_beta,
		      (int, float, float, ...));
double        * IMSL_PROTO(imsl_d_random_beta,
		      (int, double, double, ...));

	/* Chapter 12 --- Utilities */

double      	IMSL_PROTO(imsl_ctime,
		      (void));
long        	IMSL_PROTO(imsl_error_code,
		      (void));
void        	IMSL_PROTO(imsl_error_options,
		      (int, ...));
int         	IMSL_PROTO(imsl_i_machine,
		      (int));
float       	IMSL_PROTO(imsl_f_machine,
		      (int));
double      	IMSL_PROTO(imsl_d_machine,
		      (int));
float	    	IMSL_PROTO(imsl_f_constant,
		      (char*, char*));
double	    	IMSL_PROTO(imsl_d_constant,
		      (char*, char*));
void        	IMSL_PROTO(imsl_days_to_date,
		      (int, int*, int*, int*));
int         	IMSL_PROTO(imsl_date_to_days,
		      (int, int, int));
char          * IMSL_PROTO(imsl_version,
		      (int));

	/* Complex function */

f_complex   	IMSL_PROTO(imsl_c_neg,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_neg,
		      (d_complex ));
f_complex   	IMSL_PROTO(imsl_c_add,
		      (f_complex, f_complex));
d_complex   	IMSL_PROTO(imsl_z_add,
		      (d_complex, d_complex));
f_complex   	IMSL_PROTO(imsl_c_sub,
		      (f_complex, f_complex));
d_complex   	IMSL_PROTO(imsl_z_sub,
		      (d_complex, d_complex));
f_complex   	IMSL_PROTO(imsl_c_mul,
		      (f_complex, f_complex));
d_complex   	IMSL_PROTO(imsl_z_mul,
		      (d_complex, d_complex));
f_complex   	IMSL_PROTO(imsl_c_div,
		      (f_complex, f_complex));
d_complex   	IMSL_PROTO(imsl_z_div,
		      (d_complex, d_complex));
int         	IMSL_PROTO(imsl_c_eq,
		      (f_complex, f_complex));
int         	IMSL_PROTO(imsl_z_eq,
		      (d_complex, d_complex));
f_complex   	IMSL_PROTO(imsl_cz_convert,
		      (d_complex));
d_complex   	IMSL_PROTO(imsl_zc_convert,
		      (f_complex));
float       	IMSL_PROTO(imsl_c_aimag,
		      (f_complex));
double      	IMSL_PROTO(imsl_z_aimag,
		      (d_complex));
float       	IMSL_PROTO(imsl_fc_convert,
		      (f_complex));
double      	IMSL_PROTO(imsl_dz_convert,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_cf_convert,
		      (float, float));
d_complex   	IMSL_PROTO(imsl_zd_convert,
		      (double, double));
f_complex   	IMSL_PROTO(imsl_c_conjg,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_conjg,
		      (d_complex));
float       	IMSL_PROTO(imsl_c_arg,
		      (f_complex));
double      	IMSL_PROTO(imsl_z_arg,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_c_sqrt,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_sqrt,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_c_log,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_log,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_c_exp,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_exp,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_c_sin,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_sin,
		      (d_complex));
f_complex   	IMSL_PROTO(imsl_c_cos,
		      (f_complex));
d_complex   	IMSL_PROTO(imsl_z_cos,
		      (d_complex));
float       	IMSL_PROTO(imsl_c_abs,
		      (f_complex));
double      	IMSL_PROTO(imsl_z_abs,
		      (d_complex));

#endif /* USE_IMSL */

	/* Keywords */

enum Imsl_keyword {
    IMSL_TRANSPOSE	            = 10001,
    IMSL_RESULT		            = 10002,
    IMSL_A_COL_DIM	            = 10003,
    IMSL_FACTOR		            = 10004,
    IMSL_FAC_COL_DIM	            = 10005,
    IMSL_FACTOR_ONLY	            = 10006,
    IMSL_SOLVE_ONLY	            = 10007,
    
    IMSL_BACKWARD	            = 10008,
    IMSL_PARAMS		            = 10009,

    IMSL_ERR_ABS                    = 10010,
    IMSL_ERR_REL                    = 10011,
    IMSL_ETA                        = 10012,
    IMSL_EPS                        = 10013,
    IMSL_GUESS                      = 10014,
    IMSL_ITMAX                      = 10016,
    IMSL_INFO                       = 10017,
    IMSL_NUM_ROOTS                  = 10018,
    
    IMSL_RULE		            = 10019,
    IMSL_ERR_EST	            = 10020,
    IMSL_MAX_SUBINTER	            = 10021,
    IMSL_N_SUBINTER	            = 10022,
    IMSL_N_EVALS	            = 10023,
    IMSL_ERR_LIST	            = 10024,
    IMSL_ERR_ORDER	            = 10025,
    
    IMSL_BREAKPOINTS                = 10026,
    IMSL_COEFS                      = 10027,
    IMSL_DERIV                      = 10028,
    IMSL_CONCAVE                    = 10029,
    IMSL_PERIODIC                   = 10030,
    IMSL_LEFT                       = 10031,
    IMSL_RIGHT                      = 10032,
    IMSL_WEIGHT                     = 10033,
    IMSL_SMPAR                      = 10034,
    IMSL_KNOTS                      = 10035,
    IMSL_ORDER                      = 10036,
    IMSL_OPT                        = 10037,
    IMSL_MIN_PROJECTION             = 10038,
    IMSL_KNOTS_USER                 = 10039,
    
    IMSL_PRINT_ALL                  = 10040,
    IMSL_ROW_LABELS                 = 10042,
    IMSL_COL_LABELS                 = 10043,
    IMSL_WRITE_FORMAT               = 10044,
    
    IMSL_X_COL_DIM                  = 10045,
    IMSL_COV_COL_DIM                = 10046,
    IMSL_COEF_COVARIANCES           = 10047,
    IMSL_COEF_COVARIANCES_USER      = 10048,
    IMSL_RANK                       = 10049,
    IMSL_NO_INTERCEPT               = 10050,
    IMSL_ANOVA_TABLE                = 10051,
    IMSL_ANOVA_TABLE_USER           = 10052,
    IMSL_TOLERANCE	            = 10053,
    
    IMSL_FREQUENCIES                = 10054,
    IMSL_FREQUENCIES_USER           = 10055,
    IMSL_BOUNDS                     = 10056,
    IMSL_N_PARAMETERS_ESTIMATED     = 10057,
    IMSL_CUTPOINTS                  = 10058,
    IMSL_CUTPOINTS_USER             = 10059,
    IMSL_CELL_COUNTS                = 10060,
    IMSL_CELL_COUNTS_USER           = 10061,
    IMSL_CELL_EXPECTED              = 10062,
    IMSL_CELL_EXPECTED_USER         = 10063,
    IMSL_CELL_CHI_SQUARED           = 10064,
    IMSL_CELL_CHI_SQUARED_USER      = 10065,
    IMSL_DEGREES_OF_FREEDOM         = 10066,
    
    IMSL_TIES_OPTION                = 10067,
    IMSL_FUZZ                       = 10068,
    IMSL_SCORE_OPTION               = 10069,
    
    IMSL_RESULT_USER                = 10070,
    
    IMSL_NORM                       = 10071,
    IMSL_TOL                        = 10072,
    IMSL_HINIT		            = 10073,
    IMSL_HMIN		            = 10074,
    IMSL_SCALE		            = 10075,
    IMSL_FLOOR		            = 10076,
    IMSL_MAX_NUMBER_STEPS           = 10077,
    IMSL_MAX_NUMBER_FCN_EVALS       = 10078,
    IMSL_INTERRUPT_1	            = 10079,
    IMSL_INTERRUPT_2	            = 10080,
    IMSL_NSTEP		            = 10081,
    IMSL_NFCN		            = 10082,
    IMSL_HTRIAL		            = 10083,
    IMSL_VNORM		            = 10084,
    IMSL_HMAX                       = 10085,
    
    IMSL_STAT_COL_DIM               = 10086,
    IMSL_CONFIDENCE_MEANS           = 10087,
    IMSL_CONFIDENCE_VARIANCES       = 10088,
    
    IMSL_MEANS                      = 10090,
    IMSL_MEANS_USER                 = 10091,
    IMSL_COVARIANCE_COL_DIM         = 10092,
    IMSL_COMPUTE_OPTION             = 10093,
    
    IMSL_VECTORS	            = 10094,
    IMSL_VECTORS_USER	            = 10095,
    IMSL_EVECU_COL_DIM	            = 10096,
    IMSL_RANGE		            = 10097,
    
    IMSL_INFO_USER                  = 10099,
    IMSL_XGUESS                     = 10100,
    IMSL_STEP                       = 10101,
    IMSL_BOUND                      = 10102,
    IMSL_MAX_FCN                    = 10103,
    
    IMSL_INIT_TRUST_REGION          = 10104,
    IMSL_GRAD                       = 10105,
    IMSL_XSCALE                     = 10106,
    IMSL_FSCALE                     = 10107,
    IMSL_GRAD_TOL                   = 10108,
    IMSL_STEP_TOL                   = 10109,
    IMSL_REL_FCN_TOL                = 10110,
    IMSL_MAX_STEP                   = 10111,
    IMSL_GOOD_DIGIT                 = 10112,
    IMSL_MAX_ITN                    = 10113,
    IMSL_MAX_GRAD                   = 10114,
    IMSL_INIT_HESSIAN               = 10115,
    IMSL_FVALUE                     = 10116,
    IMSL_GVALUE                     = 10117,
    IMSL_JACOBIAN                   = 10118,
    IMSL_FNORM                      = 10119,
    
    IMSL_RESULT_DUAL                = 10120,
    IMSL_UPPER_LIMIT                = 10121,
    IMSL_CONSTR_TYPE                = 10122,
    IMSL_LOWER_BOUND                = 10123,
    IMSL_UPPER_BOUND                = 10124,
    IMSL_OBJ                        = 10125,
    IMSL_DUAL_USER                  = 10126,
    IMSL_DUAL                       = 10127,
    IMSL_JAC_TOL                    = 10128,
    
    IMSL_H_COL_DIM                  = 10130,
    IMSL_ADD_TO_DIAG_H              = 10131,
    
    IMSL_FDATA_COL_DIM              = 10140,
    IMSL_WEIGHTS                    = 10141,
    IMSL_ERROR_SSQ                  = 10142,
    IMSL_OPTIMIZE                   = 10143,
    IMSL_INTERCEPT                  = 10144,
    IMSL_SSE                        = 10145,
    IMSL_SMOOTHING_PAR              = 10146,
    IMSL_SUR_COL_DIM                = 10147,
    IMSL_OPT_ITMAX                  = 10148,
    IMSL_CONCAVE_ITMAX              = 10149,
    
    IMSL_SOLUTION_USER              = 10150,
    IMSL_FACTOR_USER                = 10151,
    IMSL_INVERSE                    = 10152,
    IMSL_INVERSE_USER               = 10153,
    IMSL_INV_COL_DIM                = 10154,
    IMSL_INVERSE_ONLY               = 10155,
    
    IMSL_SSQ_POLY                   = 10156,
    IMSL_SSQ_LOF                    = 10157,
    IMSL_X_MEAN                     = 10158,
    IMSL_Y_MEAN                     = 10159,
    IMSL_X_VARIANCE                 = 10160,
    IMSL_Y_VARIANCE                 = 10161,
    
    IMSL_BASIS                      = 10170,
    IMSL_PIVOT                      = 10171,
    IMSL_RESIDUAL                   = 10172,
    IMSL_RESIDUAL_USER              = 10173,
    IMSL_Q                          = 10177,
    IMSL_Q_USER                     = 10178,
    IMSL_Q_COL_DIM                  = 10179,
    IMSL_P_COL_DIM                  = 10180,
    
    IMSL_A_MATRIX                   = 10181,
    IMSL_B_MATRIX                   = 10182,
    IMSL_X_VECTOR                   = 10183,
    IMSL_Y_VECTOR                   = 10184,
    IMSL_RETURN_COL_DIM             = 10185,
    IMSL_B_COL_DIM                  = 10186,
    
    IMSL_INVA_COL_DIM               = 10187,
    
    IMSL_SET_PRINT                  = 10188,
    IMSL_GET_PRINT                  = 10189,
    IMSL_SET_STOP                   = 10190,
    IMSL_GET_STOP                   = 10191,
    IMSL_GET_TRACEBACK              = 10192,
    IMSL_SET_TRACEBACK              = 10193,
    IMSL_SET_ERROR_FILE             = 10194,
    IMSL_ERROR_PRINT_PROC           = 10195,
    IMSL_ERROR_MSG_PATH             = 10196,
    IMSL_ERROR_MSG_NAME             = 10197,
    
    IMSL_S_USER                     = 10198,
    IMSL_U                          = 10199,
    IMSL_U_USER                     = 10200,
    IMSL_U_COL_DIM                  = 10201,
    IMSL_V                          = 10202,
    IMSL_V_USER                     = 10203,
    IMSL_V_COL_DIM                  = 10204,
    
    IMSL_SET_OUTPUT_FILE            = 10208,
    IMSL_GET_OUTPUT_FILE            = 10209,
    IMSL_GET_ERROR_FILE             = 10210,
    
    IMSL_PRINT_LOWER                = 10211,
    IMSL_PRINT_UPPER                = 10212,
    IMSL_PRINT_LOWER_NO_DIAG        = 10213,
    IMSL_PRINT_UPPER_NO_DIAG        = 10214,
    IMSL_ROW_NUMBER_ZERO            = 10215,
    IMSL_NO_ROW_LABELS              = 10216,
    IMSL_COL_NUMBER_ZERO            = 10217,
    IMSL_NO_COL_LABELS              = 10218,
    
    IMSL_VARIANCE_COVARIANCE_MATRIX = 10219,
    IMSL_CORRECTED_SSCP_MATRIX      = 10220,
    IMSL_CORRELATION_MATRIX         = 10221,
    IMSL_STDEV_CORRELATION_MATRIX   = 10222,
    
    IMSL_AVERAGE_TIE  		    = 10223,
    IMSL_HIGHEST                    = 10224,
    IMSL_LOWEST                     = 10225,
    IMSL_RANDOM_SPLIT               = 10226,
    
    IMSL_RANKS			    = 10227,
    IMSL_BLOM_SCORES		    = 10228,
    IMSL_TUKEY_SCORES		    = 10229,
    IMSL_VAN_DER_WAERDEN_SCORES     = 10230,
    IMSL_EXPECTED_NORMAL_SCORES     = 10231,
    IMSL_SAVAGE_SCORES              = 10232,
    
    IMSL_CHEBYSHEV_FIRST            = 10240,
    IMSL_CHEBYSHEV_SECOND           = 10241,
    IMSL_HERMITE                    = 10242,
    IMSL_COSH                       = 10243,
    IMSL_JACOBI                     = 10244,
    IMSL_GEN_LAGUERRE               = 10245,
    IMSL_FIXED_POINT                = 10246,
    IMSL_TWO_FIXED_POINTS           = 10247,
    IMSL_LEGENDRE                   = 10248,
    
    IMSL_GRADIENT                   = 10250,
    IMSL_PRINT                      = 10251,
    
    IMSL_FJAC_Q                     = 10253,
    IMSL_FJAC_Q_USER                = 10254,
    IMSL_FJAC_R                     = 10255,
    IMSL_FJAC_R_USER                = 10256,
    
    IMSL_RETURN_NUMBER	            = 10259,
    IMSL_RETURN_USER                = 10260,
    
    IMSL_X_MEAN_USER                = 10261,
    IMSL_DF_PURE_ERROR              = 10262,
    IMSL_SSQ_PURE_ERROR             = 10263,
    IMSL_SSQ_POLY_USER              = 10264,
    IMSL_SSQ_POLY_COL_DIM           = 10265,
    IMSL_SSQ_LOF_USER               = 10266,
    IMSL_SSQ_LOF_COL_DIM            = 10267,
    IMSL_ROW_NUMBER                 = 10268,
    IMSL_COL_NUMBER                 = 10269,
    
    IMSL_CONDITION	            = 10270,
    IMSL_MAX_MOMENTS                = 10271,
    
    IMSL_OS_VERSION                 = 10272,
    IMSL_COMPILER_VERSION           = 10273,
    IMSL_LIBRARY_VERSION            = 10274,
    
    IMSL_MAX_CYCLES                 = 10275,
    IMSL_N_CYCLES                   = 10276,
    IMSL_MAX_EVALS                  = 10277,
    
    IMSL_CUTPOINTS_EQUAL            = 10280,
    
    IMSL_ABS_FCN_TOL                = 10290,
    IMSL_MAX_JACOBIAN               = 10291,
    IMSL_INTERN_SCALE               = 10292,
    IMSL_FVEC                       = 10293,
    IMSL_FVEC_USER                  = 10294,
    IMSL_FJAC                       = 10295,
    IMSL_FJAC_USER                  = 10296,
    IMSL_FJAC_COL_DIM               = 10297,
    IMSL_JTJ_INVERSE                = 10298,
    IMSL_JTJ_INVERSE_USER           = 10299,
    IMSL_JTJ_INV_COL_DIM            = 10300,
    IMSL_FULL_TRACEBACK             = 10301,
    
    IMSL_MAX_ITER                   = 10305,
    IMSL_PRECOND                    = 10306,
    IMSL_REL_ERR                    = 10307,
    
    IMSL_LICENSE_NUMBER             = 10308,
    
    IMSL_METHOD		            = 10309,
    IMSL_MAXORD		            = 10310,
    IMSL_MITER		            = 10311,
    IMSL_NFCNJ		            = 10312,
    
    IMSL_MEDIAN                     = 10320,
    IMSL_MEDIAN_AND_SCALE           = 10321,
    IMSL_CHI_SQUARED                = 10322
};
#endif /* IMSL_H   */
