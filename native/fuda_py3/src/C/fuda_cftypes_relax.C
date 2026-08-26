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
#include <cmath>
#include <cstdio>
#include <string>
#include <list>
#include <vector>
#include <complex>
#include "fuda_classes.H"
#include "fuda_utils.H"
#include "fuda_cftypes_relax.H"

/* Relaxation utility functions */
namespace RLX
{
  /* General Model-free spectral density function of Lipari & Szabo
     with all three terms. See e.g. Evenaes, Forsen, Malmendal & Akke,
     J.Mol.Biol 289, 603-17 (1999).

     The simpler model-free spectral density functions can be made
     from this as follows:

     Reduced spectral density function (omega,tau_c,S2_f): tau_f=0.0,
     S2_s=1.0, tau_s=0.0.

     Simple model-free spectral density function
     (omega,tau_c,S2_f,tau_f): S2_s=1.0, tau_s=0.0.

     Extended model-free spectral density function
     (omega,tau_c,S2_f,S2_s,tau_s): tau_f=0.0.

  */

  double sdens_iso_mf (double omega,
		       double tau_c,
		       double S2_f,
		       double tau_f,
		       double S2_s,
		       double tau_s)
  {
    double tau_f1 = fabs(tau_f*tau_c)/(fabs(tau_f)+fabs(tau_c));
    double tau_s1 = fabs(tau_s*tau_c)/(fabs(tau_s)+fabs(tau_c));
      
    double sdens = (2.0/5.0)
      *( S2_f*S2_s*tau_c/(1.0+pow(omega*tau_c,2))
	 +(1.0-S2_f)*tau_f1/(1.0+pow(omega*tau_f1,2))
	 +S2_f*(1.0-S2_s)*tau_s1/(1.0+pow(omega*tau_s1,2)) );
    
    return sdens;
  }

    
  /* Calculate cos^2 to angle between a direction vector
     (e.g. principal axis of rotational diffusion tensor) given as
     spherical polar angles, theta and phi, and a vector (e.g. a bond
     vector). */
  double calc_cos2(double theta,
		   double phi,
		   double rx,
		   double ry,
		   double rz)
  {
    // Return square of the cosine to the angle between the bond
    // vector and the principal axis of the diffusion tensor as the
    // square of the dot product of the two normalized vectors.
    return pow(rx*sin(theta)*cos(phi)
	       +ry*sin(theta)*sin(phi)
	       +rz*cos(theta),2)/(pow(rx,2)+pow(ry,2)+pow(rz,2));    
  }
  

  /* Calculate (A) scale factors and correlation times (tau) for axial
     symmetic top spectral density functions from the isotropic
     correlation time, the anisotropy and cos^2 of the angle, theta,
     between the bond vector and the principal axis. The notation used
     is identical to that of: M.Andrec, G.T.Montelione and R.M.Levy;
     J. Magn. Reson. 139, 408-21 (1999).  */
  void A_tau_ax (double tau_c,
		 double anisotropy,
		 double cos2,
		 double *A0,
		 double *A1,
		 double *A2,
		 double *tau0,
		 double *tau1,
		 double *tau2)
  {
    double sin2, D_iso, D_vertical, D_parallel;
    
    /* First calculate theta-dependent scale factors A0, A1 and A2. */
    sin2 = 1.0-cos2;
    *A0 = 0.25*pow(3.0*cos2-1,2);
    *A1 = 3.0*sin2*cos2;
    *A2 = 0.75*pow(sin2,2);
    
    /* Then the corresponding correlation times. */
    D_iso = 1.0/(6.0*tau_c);
    D_vertical = 3.0*D_iso/(2+anisotropy);
    D_parallel = anisotropy*D_vertical;
    *tau0 = 1.0/(6.0*D_vertical);
    *tau1 = 1.0/(5.0*D_vertical+D_parallel);
    *tau2 = 1.0/(2.0*D_vertical+4.0*D_parallel);
  }

  /* General Model-free spectral density function of Lipari & Szabo
     with all three terms for axial symmetric rotational
     diffusion. The theoretical validity of this general density
     function has not been verified for all physical regeimes of the
     input parameters.

     The simpler model-free spectral density functions can be made
     from this as follows:

     Reduced spectral density function (omega,tau_c,S2_f): tau_f=0.0,
     S2_s=1.0, tau_s=0.0.

     Simple model-free spectral density function
     (omega,tau_c,S2_f,tau_f): S2_s=1.0, tau_s=0.0.

     Extended model-free spectral density function
     (omega,tau_c,S2_f,S2_s,tau_s): tau_f=0.0.

     See M.Andrec, G.T.Montelione and R.M.Levy; J. Magn. Reson. 139,
     408-21 (1999) for a discussion of the validity for the model-free
     spectral density function regeime: tau_s=0.0, S2_s=1.0.

     See Kristensen et al. Mol. Biol. 299, 771-788 (2000) for a
     discussion of the validity for the extended model-free spectral
     density function regeime: tau_f=0. */

  double sdens_ax_mf (double omega,
		      double tau_c,
		      double S2_f,
		      double tau_f,
		      double S2_s,
		      double tau_s,
		      double anisotropy,
		      double cos2)
  {
    double sdens, A0, A1, A2, tau0, tau1, tau2;
    
    /* Convert to scale factors and correlation times */
    A_tau_ax (tau_c, anisotropy, cos2, &A0, &A1, &A2,
	      &tau0, &tau1, &tau2);
        
    /* Calculate spectral density. */
    sdens = ( A0*sdens_iso_mf(omega, tau0, S2_f, tau_f, S2_s, tau_s)
	      +A1*sdens_iso_mf(omega, tau1, S2_f, tau_f, S2_s, tau_s)
	      +A2*sdens_iso_mf(omega, tau2, S2_f, tau_f, S2_s, tau_s)
	      );
    
    return sdens;
  }
    

  /* Calculate the c_dipole scale factor (called d_00 in Cavanagh et al.
     (1996), table 5.5) */

  double c_dipole (double gamma_i,
		   double gamma_s,
		   double r_is)
  {
    double c = (MU_VACUUM/(4.0*PI))
      *gamma_i*gamma_s
      *(H_PLANCK/(2.0*PI))
      *(1.0/pow(r_is,3));
    c = pow(c,2);
    
    return c;
  }
  
  
  double c_csa (double omega,
		double delta_csa)
    
    /* Calculate the c_csa scale factor (called d_00 in Cavanagh et al.
       (1996), table 5.7) */    
  {
    double c =  (1.0/3.0)*pow(omega*delta_csa,2);
    
    return c;
  }


  /* Calculate the R1 relaxation rate of the S spin including contributions
     from dipole-dipole interaction and csa relaxation (Cavanagh et al.
     (1996), table 5.5+5.7) */

  double r1_s (double j0,
	       double ji,
	       double js,
	       double ji_minus_s,
	       double ji_plus_s,
	       double c_dipole,
	       double c_csa)
  {
    double r1_s;
    
    r1_s = (c_dipole/4.0)*(ji_minus_s+3.0*js+6.0*ji_plus_s) + c_csa*js;
    
    return r1_s;
  }
  
  
  /* Calculate the R2 relaxation rate of the S spin including contributions
     from dipole-dipole interaction and csa relaxation (Cavanagh et al.
     (1996), table 5.5+5.7) */

  double r2_s (double j0,
	       double ji,
	       double js,
	       double ji_minus_s,
	       double ji_plus_s,
	       double c_dipole,
	       double c_csa)
  {
    double r2_s;
    
    r2_s = (c_dipole/8.0)*( 4.0*j0+ji_minus_s+3.0*js+6.0*ji+6.0*ji_plus_s )
      + (c_csa/6.0)*( 4.0*j0+3.0*js );
    
    return r2_s;
  }


  /* Calculate the cross-relaxation rate, rho_is, (Cavanagh et al.  (1996),
     table 5.5) */

  double rho_is (double j0,
		 double ji,
		 double js,
		 double ji_minus_s,
		 double ji_plus_s,
		 double c_dipole,
		 double c_csa)
  {
    double rho_is;
    
    rho_is = (c_dipole/4.0)*( -ji_minus_s+6.0*ji_plus_s );
    
    return rho_is;
  }
  

  /* Calculate the NOE of the S spin defined as <Sz(steady-state)>/<S_0>,
     where <Sz(steady-state)> is the steady-state S_z magnetization in the
     presence of saturation of the I magnetization and <S_0> is the
     equlibrium S_z magnetization in the absence of I saturation (Cavanagh et
     al.  (1996), p.287-290) */
  double noe_s (double j0,
		double ji,
		double js,
		double ji_minus_s,
		double ji_plus_s,
		double c_dipole,
		double c_csa,
		double gamma_i,
		double gamma_s)
  {
    double noe, r1, rho;
    
    rho = rho_is (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    r1 = r1_s (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    noe = 1.0 + (rho/r1)*(gamma_i/gamma_s);
    
    return noe;
  }
  

  /* Calculate R2/R1 ratio for the S spin with the isotropic reduced
     model-free spectral density function */

  double r2r1_s_iso_mf0 (double omega_i,
			 double omega_s,
			 double c_dipole,
			 double c_csa,
			 double tau_c)
  {
    double r1, r2, ratio, j0, ji, js, ji_minus_s, ji_plus_s;
    
    /* Calculate spectral density functions. */
    j0 = sdens_iso_mf(0.0,tau_c,1.0,0.0,1.0,0.0);
    ji = sdens_iso_mf(omega_i,tau_c,1.0,0.0,1.0,0.0);
    js = sdens_iso_mf(omega_s,tau_c,1.0,0.0,1.0,0.0);
    ji_minus_s = sdens_iso_mf(omega_i-omega_s,tau_c,1.0,0.0,1.0,0.0);
    ji_plus_s = sdens_iso_mf(omega_i+omega_s,tau_c,1.0,0.0,1.0,0.0);
    
    /* Canclulate r1 and r2 and ratio */
    r1 = r1_s (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    r2 = r2_s (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    ratio = r2/r1;
    
    return ratio;
  }
  
  
  /* Calculate R2/R1 ratio for the S spin with the reduced model-free spectral
     density function. */
  double r2r1_s_ax_mf0 (double omega_i,
			double omega_s,
			double c_dipole,
			double c_csa,
			double tau_c,
			double anisotropy,
			double cos2)
  {
    double r1, r2, ratio, j0, ji, js, ji_minus_s, ji_plus_s;
    
    /* Calculate spectral density functions. */
    j0 = sdens_ax_mf(0.0,tau_c,1.0,0.0,1.0,0.0,anisotropy,cos2);
    ji = sdens_ax_mf(omega_i,tau_c,1.0,0.0,1.0,0.0,anisotropy,cos2);
    js = sdens_ax_mf(omega_s,tau_c,1.0,0.0,1.0,0.0,anisotropy,cos2);
    ji_minus_s = sdens_ax_mf(omega_i-omega_s,tau_c,1.0,0.0,1.0,0.0,
			      anisotropy,cos2);
    ji_plus_s = sdens_ax_mf(omega_i+omega_s,tau_c,1.0,0.0,1.0,0.0,
			     anisotropy,cos2);
    
    /* Canclulate r1 and r2 and ratio */
    r1 = r1_s (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    r2 = r2_s (j0,ji,js,ji_minus_s,ji_plus_s,c_dipole,c_csa);
    ratio = r2/r1;
    
    return ratio;
  }
  

  /* Calculate R2/R1 ratio for the S spin with the reduced model-free
     spectral density function as a function of the direction of the
     principal axis of the diffusion tensor in the molecular frame
     given by the angles theta and phi and the direction vector of the
     I-S bond vector. */
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
			  double gamma_s)
  {
    // Calculate angular frequencies.
    double omega_i = -gamma_i*B0;
    double omega_s = -gamma_s*B0;

    // Calculate c_dipole and c_csa.
    double C_dipole = c_dipole(gamma_i, gamma_s, r_is);
    double C_csa = c_csa(omega_s, delta_csa);

    // Calculate square of the cosine to the angle between the bond
    // vector and the principal axis of the diffusion tensor as the
    // square of the dot product of the two normalized vectors.
    double cos2 = pow(rx*sin(theta)*cos(phi)
		      +ry*sin(theta)*sin(phi)
		      +rz*cos(theta),2)
      /(pow(rx,2)+pow(ry,2)+pow(rz,2));
    
    // Call r2r1_s_ax_mf0 with the calculated value of cos2.
    double ratio = r2r1_s_ax_mf0(omega_i, omega_s,
				 C_dipole, C_csa,
				 tau_c, anisotropy, cos2);
    
    return ratio;
  }
  

  
  /* Calculate r1 from model-free spectral density-function
     parameters */

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
		     double gamma_s)
  {
    //Calculate help constants.
    double om_i = -gamma_i*B0;
    double om_s = -gamma_s*B0;
    double C_dipole = c_dipole(gamma_i, gamma_s, r_is);
    double C_csa = c_csa(om_s, delta_csa);
  
    double j0 =
      sdens_ax_mf(0.0,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji =
      sdens_ax_mf(om_i,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double js = 
      sdens_ax_mf(om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_m_s = 
      sdens_ax_mf(om_i-om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_p_s = 
      sdens_ax_mf(om_i+om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);

    double r1 = r1_s (j0,ji,js,ji_m_s,ji_p_s,C_dipole,C_csa);

    return (r1);
  }


  /* Calculate r2 from model-free spectral density-function
     parameters */

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
		     double gamma_s)
  {
    //Calculate help constants.
    double om_i = -gamma_i*B0;
    double om_s = -gamma_s*B0;
    double C_dipole = c_dipole(gamma_i, gamma_s, r_is);
    double C_csa = c_csa(om_s, delta_csa);
  
    double j0 =
      sdens_ax_mf(0.0,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji =
      sdens_ax_mf(om_i,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double js = 
      sdens_ax_mf(om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_m_s = 
      sdens_ax_mf(om_i-om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_p_s = 
      sdens_ax_mf(om_i+om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);

    double r2 = r2_s (j0,ji,js,ji_m_s,ji_p_s,C_dipole,C_csa);

    return (r2);
  }


  /* Calculate noe from model-free spectral density-function
     parameters */

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
		      double gamma_s)
  {
    //Calculate help constants.
    double om_i = -gamma_i*B0;
    double om_s = -gamma_s*B0;
    double C_dipole = c_dipole(gamma_i, gamma_s, r_is);
    double C_csa = c_csa(om_s, delta_csa);
  
    double j0 =
      sdens_ax_mf(0.0,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji =
      sdens_ax_mf(om_i,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double js = 
      sdens_ax_mf(om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_m_s = 
      sdens_ax_mf(om_i-om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    double ji_p_s = 
      sdens_ax_mf(om_i+om_s,tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2);
    
    double noe = noe_s (j0,ji,js,ji_m_s,ji_p_s,
			C_dipole,C_csa,gamma_i,gamma_s);

    return (noe);
  }


  /* Calculate r2 for cpmg sequence in the fast-exchange limit according to
     formula 10 in: Mandel, Akke and Palmer, Biochemistry 35, 16009-23
     (1996) */

  double r2_cpmg_fast_1(double tau_cp,
                        double r2_0,
                        double disp,
                        double k_ex)
    
  {
    // Calculate exponent and branch out according to magnitude.
    double z = fabs(k_ex*tau_cp/2.0);
    double r2;
    if (z < 1.0e-6)
      {
        r2 = r2_0;
      }
    else if (z < 100.0)
      {
        r2 = r2_0 + (disp/k_ex)*(1.0-tanh(z)/z);
      }
    else
      {
        r2 = r2_0 + (disp/k_ex)*(1.0-1.0/z);
      }
    return (r2);
  }


  /* Calculate r2 for cpmg sequence according exact expression for all
     timescales as expressed in Eq. (3) in Millet, Loria, Kroenke, Pons and
     Palmer, JACS 122, 2867 (2000).

     See also Davis et al. JMR B104 266-275 (1994) referenced in
     Millet et al. 

     NOTE: According to Tollinger, Skrynnikov, Mulder, Forman-Kay &
     Kay, JACS 123, 11341 (2001), this expression does NOT give
     correct oscillatory behaviour in the slow-exchange region and
     they present a better alternative which is here coded as
     r2_cpmg_2site_slow.  

     NOTE: tau_cp is the distance between the centers of two
     subsequent 180 degree refocussing pulses in the CPMG sequence. */
  
  double r2_cpmg_2site(double tau_cp,
                       double R2_a,
                       double R2_b,
                       double p_a,
                       double k_ex,
                       double delta_cs,
		       double B0,
		       double gamma)
  {
    // Calculate some help variables defined in Millet et al., eqs. (4)-(7).
    const double sqrt2 = sqrt(2.0);
    double p_b = 1.0-p_a;
    double delta_omega = fabs(1.0e-6*delta_cs*B0*gamma);
    double psi = pow(R2_a-R2_b-p_a*k_ex+p_b*k_ex,2)
      - pow(delta_omega,2) + 4.0*p_a*p_b*pow(k_ex,2);
    double zeta = 2.0*delta_omega*(R2_a-R2_b-p_a*k_ex+p_b*k_ex);
    double eta_p = (tau_cp/sqrt2)*sqrt(psi+sqrt(pow(psi,2)+pow(zeta,2)));
    double eta_m = (tau_cp/sqrt2)*sqrt(-psi+sqrt(pow(psi,2)+pow(zeta,2)));
    double D_p = 0.5*(1.0+(psi+2.0*pow(delta_omega,2))
		      /sqrt(pow(psi,2)+pow(zeta,2)));
    double D_m = 0.5*(-1.0+(psi+2.0*pow(delta_omega,2))
		      /sqrt(pow(psi,2)+pow(zeta,2)));
    
    // Branch out according to eta_p magnitude.
    double z, acosh_z;
    if (fabs(eta_p)<100.0)
      {
	// Calculate exponent of acosh in eq. (3).
	z = D_p*cosh(eta_p)-D_m*cos(eta_m);
	
	// Calculate acosh to exponent.
	acosh_z = log(z+sqrt(pow(z,2)-1.0));
      }
    else
      {
	// We approximate to avoid overflow.
	z = 0.0;
	acosh_z = log(D_p) + fabs(eta_p);
      }
    
    // Calculate and return R2 according to eq. (3).
    double R2 = 0.5*(R2_a+R2_b+k_ex-acosh_z/tau_cp);

    /*

    printf("tau_cp: %.16g psi: %.16g zeta: %.16g eta_p: %.16g eta_m: %.16g\n",
	   tau_cp, psi, zeta, eta_p, eta_m);
    printf("D_p: %.16g D_m: %.16g z: %.16g acosh_z: %.16g\n",
	   D_p, D_m, z, acosh_z);
    printf("approx: %.16g\n", approx);

    */
    
    return (R2);
  }  


  /* Calculate r2 for cpmg sequence for slow-exchange limit according
     to Tollinger, Skrynnikov, Mulder, Forman-Kay & Kay, JACS 123,
     11341 (2001). 

     NOTE: tau_cp is the distance between the centers of two
     subsequent 180 degree refocussing pulses in the CPMG sequence. */
  
  double r2_cpmg_2site_slow(double tau_cp,
			    double R2_a,
			    double k_a,
			    double delta_omega)
  {
    double x = delta_omega*tau_cp/2.0;
    return (R2_a + k_a*(1.0 - sin(x)/x));
  }  


  double exsy_two_site(int model,
		       double tm,
                       double M0a,
                       double M0b,
                       double R1a,
                       double R1b,
                       double kba)
  {
    /* EXSY intensities for a slow two-site exchange system. model==1
       returns Iaa, model==2 returns Ibb and model==3 returns
       Iab=Iba. 

       Palmers expressions (Protein NMR spectroscopy p.294).

       The function has been tested in fuda by comparing numerical
       derivatives of the three return model values to
       the derivatives calculated from the McConnel equations.

    */

    // Define some helper values.
    double kab = (kba*M0b)/M0a; // Imposed by chemical equlibrium.
    double sigma = R1a+R1b+kab+kba;
    double delta = R1a-R1b+kab-kba;    
    double gamma_p = 0.5*(sigma + sqrt(pow(delta,2) + 4.0*kab*kba));
    double gamma_m = 0.5*(sigma - sqrt(pow(delta,2) + 4.0*kab*kba));
    
    
    // Branch out according to model.
    double intens = 0.0;
    if (model==1)
      {
	double a11 = 0.5 * ((1.0-delta/(gamma_p-gamma_m))
			    *exp(-gamma_m*tm)
			    +(1.0+delta/(gamma_p-gamma_m))
			    *exp(-gamma_p*tm));
	intens = a11*M0a;
      }
    else if (model==2)
      {
	double a22 = 0.5 * ((1.0+delta/(gamma_p-gamma_m))
			    *exp(-gamma_m*tm)
			    +(1.0-delta/(gamma_p-gamma_m))
			    *exp(-gamma_p*tm));
	intens = a22*M0b;
      }
    else if (model==3)
      {
	double a12 = (kba/(gamma_p-gamma_m))
	  *(exp(-gamma_m*tm)-exp(-gamma_p*tm));
	intens = a12*M0b;
      }
    else exit(1);
    
    return intens;
  }
  

} // End of RLX namespace



/* Relaxation functions in fuda */
namespace FUDA
{
  // Lipari-Szabo general model-free spectral density function with
  // isotropic symmetric rotational diffusion tensor.
  int sdens_iso_mf(void *fs, double p[], int dp_flg[],
		   double dp[], double *value) 
  {      
    double& I = p[0];
    double& omega = p[1];
    double& tau_c = p[2];
    double& S2_f = p[3];
    double& tau_f = p[4];
    double& S2_s = p[5];
    double& tau_s = p[6];
    double& y = *value;
    
    y = I*RLX::sdens_iso_mf(omega, tau_c, S2_f, tau_f, S2_s, tau_s);

    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    
    return (0);
  }
  

  // Lipari-Szabo general model-free spectral density function with
  // axially symmetric rotational diffusion tensor.
  int sdens_ax_mf(void *fs, double p[], int dp_flg[],
		  double dp[], double *value) 
  {      
    double& I = p[0];
    double& omega = p[1];
    double& tau_c = p[2];
    double& S2_f = p[3];
    double& tau_f = p[4];
    double& S2_s = p[5];
    double& tau_s = p[6];
    double& aniso = p[7];
    double& cos2 = p[8];    
    double& y = *value;
    
    y = I*RLX::sdens_ax_mf(omega, tau_c, S2_f, tau_f,
			   S2_s, tau_s, aniso, cos2);

    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    
    return (0);
  }
  

  // R2/R1 ratio for Lipari-Szabo simple model-free spectral density
  // function with isotropic tumbling.
  int r2r1_s_iso_mf0(void *fs, double p[], int dp_flg[],
		     double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& B0 = p[1];
    double& delta_csa = p[2];
    double& r_is = p[3];
    double& gamma_i = p[4];
    double& gamma_s = p[5];
    double& y = *value;
    
    double omega_i = -gamma_i*B0;
    double omega_s = -gamma_s*B0;
    double c_dipole = RLX::c_dipole(gamma_i, gamma_s, r_is);
    double c_csa = RLX::c_csa(omega_s, delta_csa);
    

    y = RLX::r2r1_s_iso_mf0(omega_i, omega_s, c_dipole, c_csa, tau_c);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    
    return (0);
  }
  

  // R2/R1 ratio for Lipari-Szabo simple model-free spectral density
  // function with axially symmetric rotational diffusion tensor.
  int r2r1_s_ax_mf0_b(void *fs, double p[], int dp_flg[],
		      double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& anisotropy = p[1];
    double& theta = p[2];
    double& phi = p[3];
    double& rx = p[4];
    double& ry = p[5];
    double& rz = p[6];    
    double& B0 = p[7];
    double& delta_csa = p[8];
    double& r_is = p[9];
    double& gamma_i = p[10];
    double& gamma_s = p[11];
    double& y = *value;
    
    y = RLX::r2r1_s_ax_mf0_b(tau_c,anisotropy,theta,phi,rx,ry,rz,
			     B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    
    return (0);
  }
  

  // R1 for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int r1_s_ax_mf(void *fs, double p[], int dp_flg[],
		 double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& cos2 = p[6];    
    double& B0 = p[7];
    double& delta_csa = p[8];
    double& r_is = p[9];
    double& gamma_i = p[10];
    double& gamma_s = p[11];
    double& y = *value;
    
    y = RLX::r1_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    
    return (0);
  }
  

  // R2 for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int r2_s_ax_mf(void *fs, double p[], int dp_flg[],
		 double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& cos2 = p[6];    
    double& B0 = p[7];
    double& delta_csa = p[8];
    double& r_is = p[9];
    double& gamma_i = p[10];
    double& gamma_s = p[11];
    double& y = *value;
    
    y = RLX::r2_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    
    return (0);
  }
  

  // NOE for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int noe_s_ax_mf(void *fs, double p[], int dp_flg[],
		  double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& cos2 = p[6];    
    double& B0 = p[7];
    double& delta_csa = p[8];
    double& r_is = p[9];
    double& gamma_i = p[10];
    double& gamma_s = p[11];
    double& y = *value;
    
    y = RLX::noe_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			 B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    
    return (0);
  }
  

  // R1 for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int r1_s_ax_mf_b(void *fs, double p[], int dp_flg[],
		   double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& theta = p[6];
    double& phi = p[7];
    double& rx = p[8];
    double& ry = p[9];
    double& rz = p[10];    
    double& B0 = p[11];
    double& delta_csa = p[12];
    double& r_is = p[13];
    double& gamma_i = p[14];
    double& gamma_s = p[15];
    double& y = *value;

    double cos2 = RLX::calc_cos2(theta,phi,rx,ry,rz);

    y = RLX::r1_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    if (dp_flg[12]) return(256*12+1);
    if (dp_flg[13]) return(256*13+1);
    if (dp_flg[14]) return(256*14+1);
    if (dp_flg[15]) return(256*15+1);
    
    return (0);
  }
  

  // R2 for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int r2_s_ax_mf_b(void *fs, double p[], int dp_flg[],
		   double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& theta = p[6];
    double& phi = p[7];
    double& rx = p[8];
    double& ry = p[9];
    double& rz = p[10];    
    double& B0 = p[11];
    double& delta_csa = p[12];
    double& r_is = p[13];
    double& gamma_i = p[14];
    double& gamma_s = p[15];
    double& y = *value;
    
    double cos2 = RLX::calc_cos2(theta,phi,rx,ry,rz);

    y = RLX::r2_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    if (dp_flg[12]) return(256*12+1);
    if (dp_flg[13]) return(256*13+1);
    if (dp_flg[14]) return(256*14+1);
    if (dp_flg[15]) return(256*15+1);
    
    return (0);
  }
  

  // NOE for Lipari-Szabo general model-free spectral density function
  // with axially symmetric rotational diffusion tensor.
  int noe_s_ax_mf_b(void *fs, double p[], int dp_flg[],
		    double dp[], double *value) 
  {      
    double& tau_c = p[0];
    double& S2_f = p[1];
    double& tau_f = p[2];
    double& S2_s = p[3];
    double& tau_s = p[4];
    double& aniso = p[5];
    double& theta = p[6];
    double& phi = p[7];
    double& rx = p[8];
    double& ry = p[9];
    double& rz = p[10];    
    double& B0 = p[11];
    double& delta_csa = p[12];
    double& r_is = p[13];
    double& gamma_i = p[14];
    double& gamma_s = p[15];
    double& y = *value;
    
    double cos2 = RLX::calc_cos2(theta,phi,rx,ry,rz);

    y = RLX::noe_s_ax_mf(tau_c,S2_f,tau_f,S2_s,tau_s,aniso,cos2,
			 B0,delta_csa,r_is,gamma_i,gamma_s);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    if (dp_flg[9]) return(256*9+1);
    if (dp_flg[10]) return(256*10+1);
    if (dp_flg[11]) return(256*11+1);
    if (dp_flg[12]) return(256*12+1);
    if (dp_flg[13]) return(256*13+1);
    if (dp_flg[14]) return(256*14+1);
    if (dp_flg[15]) return(256*15+1);
    
    return (0);
  }
  

  // R2 with two-site exchange in CPMG, all timescales.
  int r2_cpmg_2site(void *fs, double p[], int dp_flg[],
		    double dp[], double *value) 
  {      
    double& tau_cp = p[0];
    double& R2_a = p[1];
    double& R2_b = p[2];
    double& delta_cs = p[3];
    double& p_a = p[4];
    double& k_ex = p[5];    
    double& B0 = p[6];
    double& gamma = p[7];
    double& y = *value;
    
    y = RLX::r2_cpmg_2site(tau_cp, R2_a, R2_b,
			   delta_cs, p_a, k_ex, B0, gamma);
    
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    if (dp_flg[7]) return(256*7+1);
    if (dp_flg[8]) return(256*8+1);
    
    return (0);
  }
  

  // R2 with two-site slow exchange in CPMG.
  int r2_cpmg_2site_slow(void *fs, double p[], int dp_flg[],
			 double dp[], double *value) 
  {      
    double& tau_cp = p[0];
    double& R2_a = p[1];
    double& k_a = p[2];
    double& delta_omega = p[3];
    double& y = *value;
    
    y =  RLX::r2_cpmg_2site_slow(tau_cp,
				 R2_a,
				 k_a,
				 delta_omega);
  
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    
    return (0);
  }
  

  // R2 with two-site slow exchange in CPMG.
  int r2_cpmg_2site_fast(void *fs, double p[], int dp_flg[],
			 double dp[], double *value) 
  {      
    double& tau_cp = p[0];
    double& R2 = p[1];
    double& disp = p[2];
    double& k_ex = p[3];
    double& y = *value;
    
    y =  RLX::r2_cpmg_fast_1(tau_cp,
			     R2,
			     disp,
			     k_ex);
  
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    
    return (0);
  }
  

  // R2 with two-site slow exchange between monomer and dimer state as
  // a function of total (protein) concentration c and dissociation
  // constant Kd. The Exchange term is modelled as R_ex = s * p_m *
  // p_d, where p_m and p_d are the populations of species in
  // monomeric and dimeric states and s is a scaling factor. This
  // scaling factor will depend on the chemical shift difference and
  // the exchange rate constant. p_m and p_d are calculated from the
  // equilibrium expression for dissociation as a function of the
  // total concentration c: dimer <--> 2*monomer, Kd =
  // [monomer]^2/[dimer] with c = 2*[dimer] + [monomer] and p_m =
  // [monomer]/c, p_d = 2[dimer]/c.

  int r_ex_dimer_Kd(void *fs, double p[], int dp_flg[],
			 double dp[], double *value) 
  {      
    double& scale = p[0];
    double& Kd = p[1];
    double& c = p[2];
    double& y = *value;
    
    double p_m = (-Kd+sqrt(Kd*Kd+8.0*c*Kd))/(4.0*c);
    y = scale*p_m*(1.0-p_m);
  
    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    
    return (0);
  }
  

  // EXSY two-site slow exchange.
  int exsy_two_site(void *fs, double p[], int dp_flg[],
		    double dp[], double *value) 
  {      
    double& tm = p[0];
    double& M0a = p[1];
    double& M0b = p[2];
    double& R1a = p[3];
    double& R1b = p[4];
    double& kba = p[5];    
    double& model = p[6];
    double& y = *value;
    
    // Get value.
    y = RLX::exsy_two_site((int) model, tm, M0a, M0b, R1a, R1b, kba);

    // No derrivatives are available.
    if (dp_flg[0]) return(256*0+1);
    if (dp_flg[1]) return(256*1+1);
    if (dp_flg[2]) return(256*2+1);
    if (dp_flg[3]) return(256*3+1);
    if (dp_flg[4]) return(256*4+1);
    if (dp_flg[5]) return(256*5+1);
    if (dp_flg[6]) return(256*6+1);
    
    return (0);
  }
  

  
  // Here we declare function types for the functions above.
  void declare_cftypes_relax(Fuda *fuda_ptr, std::string tag)
  {
    std::string name, descr;
    std::vector<std::string> p_name, p_descr;
    std::vector<unsigned int> p_var;
    Ftype *ft;

    // Check tag.
    if (tag!="relax" && tag!="all") return;
    
    // Get ref. to fuda object.
    Fuda& fuda = *fuda_ptr;
    
    // Declare function types one by one.
    
    // sdens_iso_mf.
    ft = fuda.ftype_add_cfunc("sdens_iso_mf", 7, sdens_iso_mf, 0, 0);
    ft->set_descr("Spectral density func. for isotropic rot. dif. tensor");
    ft->set_param(0,1,0,"scale","Scale factor");
    ft->set_param(1,1,0,"omega","angular velocity");
    ft->set_param(2,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(3,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(4,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(5,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(6,1,0,"tau_s","Eff. correlation time for slow motions");
    
    // sdens_ax_mf.
    ft = fuda.ftype_add_cfunc("sdens_ax_mf", 9, sdens_ax_mf, 0, 0);
    ft->set_descr("Spectral density func. for axially symmetric rot. dif. tensor");
    ft->set_param(0,1,0,"scale","Scale factor");
    ft->set_param(1,1,0,"omega","angular velocity");
    ft->set_param(2,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(3,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(4,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(5,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(6,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(7,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(8,1,0,"cos2","square of the cosine to angle "
		  "between r_is and principal axis");
    
    // r1_s_ax_mf.
    ft = fuda.ftype_add_cfunc("r1_s_ax_mf", 12, r1_s_ax_mf, 0, 0);
    ft->set_descr("R1 with axially symmetric rot. dif. tensor");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"cos2","square of the cosine to angle "
		  "between r_is and principal axis");
    ft->set_param(7,1,0,"B0","Magnetic field strength");
    ft->set_param(8,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(9,1,0,"r_is","IS inter-spin distance");
    ft->set_param(10,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(11,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r2_s_ax_mf.
    ft = fuda.ftype_add_cfunc("r2_s_ax_mf", 12, r2_s_ax_mf, 0, 0);
    ft->set_descr("R2 with axially symmetric rot. dif. tensor");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"cos2","square of the cosine to angle "
		  "between r_is and principal axis");
    ft->set_param(7,1,0,"B0","Magnetic field strength");
    ft->set_param(8,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(9,1,0,"r_is","IS inter-spin distance");
    ft->set_param(10,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(11,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // noe_s_ax_mf.
    ft = fuda.ftype_add_cfunc("noe_s_ax_mf", 12, noe_s_ax_mf, 0, 0);
    ft->set_descr("Heteronuclear NOE with axially symmetric "
		  "rot. dif. tensor");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"cos2","square of the cosine to angle "
		  "between r_is and principal axis");
    ft->set_param(7,1,0,"B0","Magnetic field strength");
    ft->set_param(8,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(9,1,0,"r_is","IS inter-spin distance");
    ft->set_param(10,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(11,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r1_s_ax_mf_b.
    ft = fuda.ftype_add_cfunc("r1_s_ax_mf_b", 16, r1_s_ax_mf_b, 0, 0);
    ft->set_descr("R1 with ax. sym. rot. dif. tensor and bond vector");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"theta","angle between z-axis and "
		  "rot. dif. tensor principal axis ");
    ft->set_param(7,1,0,"phi","angle between x-axis and "
		  "rot. dif. tensor principal axis projection on xy-plane");
    ft->set_param(8,1,0,"rx","bond vector x composant");
    ft->set_param(9,1,0,"ry","bond vector y composant");
    ft->set_param(10,1,0,"rz","bond vector z composant");
    ft->set_param(11,1,0,"B0","Magnetic field strength");
    ft->set_param(12,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(13,1,0,"r_is","IS inter-spin distance");
    ft->set_param(14,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(15,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r2_s_ax_mf_b.
    ft = fuda.ftype_add_cfunc("r2_s_ax_mf_b", 16, r2_s_ax_mf_b, 0, 0);
    ft->set_descr("R2 with ax. sym. rot. dif. tensor and bond vector");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"theta","angle between z-axis and "
		  "rot. dif. tensor principal axis ");
    ft->set_param(7,1,0,"phi","angle between x-axis and "
		  "rot. dif. tensor principal axis projection on xy-plane");
    ft->set_param(8,1,0,"rx","bond vector x composant");
    ft->set_param(9,1,0,"ry","bond vector y composant");
    ft->set_param(10,1,0,"rz","bond vector z composant");
    ft->set_param(11,1,0,"B0","Magnetic field strength");
    ft->set_param(12,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(13,1,0,"r_is","IS inter-spin distance");
    ft->set_param(14,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(15,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // noe_s_ax_mf_b.
    ft = fuda.ftype_add_cfunc("noe_s_ax_mf_b", 16, noe_s_ax_mf_b, 0, 0);
    ft->set_descr("Heteronuclear NOE with ax. sym. rot. dif. tensor and bond vector");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"S2_f","Squared order parameter for fast motions");
    ft->set_param(2,1,0,"tau_f","Eff. correlation time for fast motions");
    ft->set_param(3,1,0,"S2_s","Squared order parameter for slow motions");
    ft->set_param(4,1,0,"tau_s","Eff. correlation time for slow motions");
    ft->set_param(5,1,0,"aniso","Rot. dif. anisotropy D_par/D_per");
    ft->set_param(6,1,0,"theta","angle between z-axis and "
		  "rot. dif. tensor principal axis ");
    ft->set_param(7,1,0,"phi","angle between x-axis and "
		  "rot. dif. tensor principal axis projection on xy-plane");
    ft->set_param(8,1,0,"rx","bond vector x composant");
    ft->set_param(9,1,0,"ry","bond vector y composant");
    ft->set_param(10,1,0,"rz","bond vector z composant");
    ft->set_param(11,1,0,"B0","Magnetic field strength");
    ft->set_param(12,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(13,1,0,"r_is","IS inter-spin distance");
    ft->set_param(14,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(15,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r2r1_s_iso_mf0.
    ft = fuda.ftype_add_cfunc("r2r1_s_iso_mf0", 6, r2r1_s_iso_mf0, 0, 0);
    ft->set_descr("R2/R1 ratio with isotropic tumbling");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"B0","Magnetic field strength");
    ft->set_param(2,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(3,1,0,"r_is","IS inter-spin distance");
    ft->set_param(4,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(5,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r2r1_s_ax_mf0_b.
    ft = fuda.ftype_add_cfunc("r2r1_s_ax_mf0_b", 12, r2r1_s_ax_mf0_b, 0, 0);
    ft->set_descr("R2/R1 ratio with axially symmetric rot. dif. tensor");
    ft->set_param(0,1,0,"tau_c","Overall isotropic correlation time");
    ft->set_param(1,1,0,"aniso","Rotationa diff. tensor anisotropy");
    ft->set_param(2,1,0,"theta","Spherical polar theta coordiante "
		  "of tensor");
    ft->set_param(3,1,0,"phi","Spherical polar phi coordiante of tensor");
    ft->set_param(4,1,0,"rx","x-coordinate of I-S inter-spin vector");
    ft->set_param(5,1,0,"ry","y-coordinate of I-S inter-spin vector");
    ft->set_param(6,1,0,"rz","z-coordinate of I-S inter-spin vector");
    ft->set_param(7,1,0,"B0","Magnetic field strength");
    ft->set_param(8,1,0,"delta_csa","difference between axial and "
		  "perpendicular components of chemical shift tensor");
    ft->set_param(9,1,0,"r_is","IS inter-spin distance");
    ft->set_param(10,0,0,"gamma_i","magnetogyric ratio for I spin");
    ft->set_param(11,0,0,"gamma_s","magnetogyric ratio for S spin");
    
    // r2_cpmg_2site.
    ft = fuda.ftype_add_cfunc("r2_cpmg_2site", 8, r2_cpmg_2site, 0, 0);
    ft->set_descr("R2 with 2-site exchange in CPMG, all timescales");
    ft->set_param(0,1,0,"tau_cp","CPMG delay (delay between 180 pulses)");
    ft->set_param(1,1,0,"R2_a","R2 of site A");
    ft->set_param(2,1,0,"R2_b","R2 of site B");
    ft->set_param(3,1,0,"delta_cs","difference in chemical shift");
    ft->set_param(4,1,0,"p_a","population of site A");
    ft->set_param(5,1,0,"k_ex","exchage rate constant");
    ft->set_param(6,1,0,"B0","Magnetic field strength");
    ft->set_param(7,0,0,"gamma","Magnetogyric ratio");

    // r2_cpmg_2site_fast.
    ft = fuda.ftype_add_cfunc("r2_cpmg_2site_fast",\
			      4, r2_cpmg_2site_fast, 0, 0);
    ft->set_descr("R2 with 2-site slow exchange in CPMG");
    ft->set_param(0,1,0,"tau_cp",\
		  "CPMG delay (delay between 180 pulses)");
    ft->set_param(1,1,0,"R2","R2 of site A");
    ft->set_param(2,1,0,"disp","Dispersion factor");
    ft->set_param(3,1,0,"k_ex","exchange-rate constant");

    // r2_cpmg_2site_slow.
    ft = fuda.ftype_add_cfunc("r2_cpmg_2site_slow",\
			      4, r2_cpmg_2site_slow, 0, 0);
    ft->set_descr("R2 with 2-site slow exchange in CPMG");
    ft->set_param(0,1,0,"tau_cp",\
		  "CPMG delay (delay between 180 pulses)");
    ft->set_param(1,1,0,"R2_a","R2 of site A");
    ft->set_param(2,1,0,"k_a","forward exchange rate constant, k_a=p_b*k_ex");
    ft->set_param(3,1,0,"delta_omega","difference in lamour frequency");

    // r_ex_dimer_Kd.
    ft = fuda.ftype_add_cfunc("r_ex_dimer_Kd",\
			      3, r_ex_dimer_Kd, 0, 0);
    ft->set_descr("Monomer-dimer exchange term");
    ft->set_param(0,1,0,"scale",\
		  "Scale factor");
    ft->set_param(1,1,0,"Kd","Dissociation constant");
    ft->set_param(2,1,0,"c","Concentration");

    // exsy_two_site.
    ft = fuda.ftype_add_cfunc("exsy_two_site", 7, exsy_two_site, 0, 0);
    ft->set_descr("EXSY two-site slow exchange intensity model");
    ft->set_param(0,1,0,"tm","mixing time");
    ft->set_param(1,1,0,"M0a","equlibrium magnitization of site a");
    ft->set_param(2,1,0,"M0b","equlibrium magnetization of site b");
    ft->set_param(3,1,0,"R1a","R1 for site a");
    ft->set_param(4,1,0,"R1b","R1 for site b");
    ft->set_param(5,1,0,"kba","exchange rate constant b->a");
    ft->set_param(6,0,0,"model","Intensity model: 1:Iaa, 2:Ibb, 3:Iab");

  }
}






