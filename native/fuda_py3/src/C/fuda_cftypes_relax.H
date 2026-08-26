/*Legal_notice:

Copyright 1996-2001 by Soren M. Kristensen, Department of Chemistry,
University of Copenhagen, DK-2100 Copenhagen O, Denmark.  All rights
reserved.

This software and its related software and documentation is provided 'as
is' without express or implied warrenty. The Departmet of Chemistry, The
University of Copenhagen and the author make no warrenties as to any matter
whatsoever with respect to the program and the related software and
documentation. In particular, any and all warranties of merchantability and
fitness for any particular purpose are expressly excluded.

Legal_notice*/

namespace RLX
{
  const double PI          =   3.14159265358979323846;
  const double GAMMA_1H    = 267.512e6;         /* Unit: 1/(T*s) */
  const double GAMMA_15N   = -27.107e6;         /* Unit: 1/(T*s) */
  const double GAMMA_13C   =  67.2640e6;        /* Unit: 1/(T*s) */
  const double H_PLANCK    =   6.62618e-34;     /* Unit: J/s    */
  const double MU_VACUUM   =  12.56637061e-7;   /* Unit: N/A^2 */
  const double K_BOLTZMAN  =   1.38066e-23;     /* Unit: J/K */ 
  const double N_AVOGADRO  =   6.02205e23;      /* Unit: 1/mol */
  const double R_GAS_CONST =   8.31441;         /* Unit: J/(K*mol) */

  double calc_cos2 (double theta,
		    double phi,
		    double rx,
		    double ry,
		    double rz);

  double sdens_iso_mf (double omega,
		       double tau_c,
		       double S2_f,
		       double tau_f,
		       double S2_s,
		       double tau_s);
  
  void A_tau_ax (double tau_c,
		 double anisotropy,
		 double cos2,
		 double *A0,
		 double *A1,
		 double *A2,
		 double *tau0,
		 double *tau1,
		 double *tau2);
  
  double sdens_ax_mf (double omega,
		      double tau_c,
		      double S2_f,
		      double tau_f,
		      double S2_s,
		      double tau_s,
		      double anisotropy,
		      double cos2);
  
  double c_dipole (double gamma_i,
		   double gamma_s,
		   double r_is);
  
  double c_csa (double omega,
		double delta_csa);
  
  double r1_s (double j0,
	       double ji,
	       double js,
	       double ji_minus_s,
	       double ji_plus_s,
	       double c_dipole,
	       double c_csa);
  
  double r2_s (double j0,
	       double ji,
	       double js,
	       double ji_minus_s,
	       double ji_plus_s,
	       double c_dipole,
	       double c_csa);
  
  double rho_is (double j0,
		 double ji,
		 double js,
		 double ji_minus_s,
		 double ji_plus_s,
		 double c_dipole,
		 double c_csa);
  
  double noe_s (double j0,
		double ji,
		double js,
		double ji_minus_s,
		double ji_plus_s,
		double c_dipole,
		double c_csa,
		double gamma_i,
		double gamma_s);
  
  double r2r1_s_iso_mf0 (double omega_i,
			 double omega_s,
			 double c_dipole,
			 double c_csa,
			 double tau_m);
  
  double r2r1_s_ax_mf0 (double omega_i,
			double omega_s,
			double c_dipole,
			double c_csa,
			double tau_m,
			double anisotropy,
			double cos2);
  
  double r2r1_s_ax_mf0_b (double tau_c,
			  double anisotropy,
			  double theta,
			  double phi,
			  double rx,
			  double ry,
			  double rz,
			  double B0,
			  double delta_csa,
			  double r_is,
			  double gamma_i,
			  double gamma_s);
  
  double r1_s_ax_mf (double tau_c,
		     double S2_f,
		     double tau_f,
		     double S2_s,
		     double tau_s,
		     double aniso,
		     double cos2,
		     double B0,
		     double delta_csa,
		     double r_is,
		     double gamma_i,
		     double gamma_s);
  
  double r2_s_ax_mf (double tau_c,
		     double S2_f,
		     double tau_f,
		     double S2_s,
		     double tau_s,
		     double aniso,
		     double cos2,
		     double B0,
		     double delta_csa,
		     double r_is,
		     double gamma_i,
		     double gamma_s);
  
  double noe_s_ax_mf (double tau_c,
		      double S2_f,
		      double tau_f,
		      double S2_s,
		      double tau_s,
		      double aniso,
		      double cos2,
		      double B0,
		      double delta_csa,
		      double r_is,
		      double gamma_i,
		      double gamma_s);
  
  double r2_cpmg_fast_1(double tau_cp,
                        double r2_0,
                        double disp,
                        double k_ex);
  
  double r2_cpmg_2site(double tau_cp,
                       double R2_a,
                       double R2_b,
                       double p_a,
                       double k_ex,
                       double delta_cs,
		       double B0,
		       double gamma);

  double r2_cpmg_2site_slow(double tau_cp,
			    double R2_a,
			    double k_a,
			    double delta_omega);

  double exsy_two_site(int model,
		       double tm,
                       double M0a,
                       double M0b,
                       double R1a,
                       double R1b,
                       double kba);
  

}
 
namespace FUDA
{
  int sdens_ax_mf(void *fs, double p[], int dp_flg[],
		   double dp[], double *value);
  int r1_s_ax_mf(void *fs, double p[], int dp_flg[],
		  double dp[], double *value);
  int r2_s_ax_mf(void *fs, double p[], int dp_flg[],
		  double dp[], double *value);
  int noe_s_ax_mf(void *fs, double p[], int dp_flg[],
		   double dp[], double *value);
  int r2_cpmg_fast_1(void *fs, double p[], int dp_flg[],
		   double dp[], double *value);
  int r2_cpmg_2site(void *fs, double p[], int dp_flg[],
		   double dp[], double *value);
  int r2_cpmg_2site_slow(void *fs, double p[], int dp_flg[],
			 double dp[], double *value);
  int r2_cpmg_2site_fast(void *fs, double p[], int dp_flg[],
			 double dp[], double *value);
  int r_ex_dimer_Kd(void *fs, double p[], int dp_flg[],
			 double dp[], double *value);
  int exsy_two_site(void *fs, double p[], int dp_flg[],
		    double dp[], double *value);

  void declare_ftypes_relax(Fuda *fuda, std::string tag);
} 
