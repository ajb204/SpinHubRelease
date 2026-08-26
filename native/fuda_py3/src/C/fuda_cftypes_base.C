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
#include <string>
#include <list>
#include <vector>
#include <complex>
#include "fuda_classes.H"
#include "fuda_utils.H"

typedef std::complex<double> dcomplex;

/* Predefined functions in fuda */
namespace FUDA
{
  // Normalized polynomial term pow(x,power).
  int norm_poly(void *fs, double p[], int dp_flg[], double dp[], double *value) 
    {
      double& x = p[0];
      double& power = p[1];
      double& y = *value;

      y = pow(x,power);
      
      if (dp_flg[0]) dp[0] = power*pow(x,power-1.0);
      if (dp_flg[1]) dp[1] = y*log(x);
      
      return (0);
    }
  
  // Normalized exponential term.
  int norm_exp(void *fs, double p[], int dp_flg[], double dp[], double *value) 
    {
      double& a = p[0];
      double& x = p[1];
      double& y = *value;

      y = exp(a*x);
      
      if (dp_flg[0]) dp[0] = x*y;
      if (dp_flg[1]) dp[1] = a*y;
      
      return (0);
    }
  
  // Normalized exponential term with base number a.
  int norm_pow(void *fs, double p[], int dp_flg[], double dp[], double *value) 
    {
      double& a = p[0];
      double& b = p[1];
      double& c = p[2];
      double& x = p[3];
      double& y = *value;

      y = pow(a,b*c*x);
      
      if (dp_flg[0]) dp[0] = b*c*x*pow(a,b*c*x-1.0);
      if (dp_flg[1]) dp[1] = y*c*x*log(a);
      if (dp_flg[2]) dp[2] = y*b*x*log(a);
      if (dp_flg[3]) dp[3] = y*b*c*log(a);
      
      return (0);
    }
  
  // Normalized decaying exponential term exp(-a*x).
  int norm_exp_decay(void *fs, double p[], int dp_flg[],
		     double dp[], double *value) 
    {
      double& a = p[0];
      double& x = p[1];
      double& y = *value;

      y = exp(-a*x);
      
      if (dp_flg[0]) dp[0] = -x*y;
      if (dp_flg[1]) dp[1] = -a*y;
      
      return (0);
    }
  
  // cos(omega*t+phase). Note: in radians and angular velocity.
  int norm_cos(void *fs, double p[], int dp_flg[], double dp[], double *value) 
    {
      double& omega = p[0];
      double& t = p[1];
      double& phase = p[2];
      double& y = *value;

      double angle = omega*t+phase;
      y = cos(angle);
      
      if (dp_flg[0]) dp[0] = -t*sin(angle);
      if (dp_flg[1]) dp[1] = -omega*sin(angle);
      if (dp_flg[2]) dp[2] = -sin(angle);
      
      return (0);
    }
  
  // cos(omega*t+phase)*exp(-R*t). Note: in radians and angular velocity.
  int norm_cos_exp_decay(void *fs, double p[], int dp_flg[],
			 double dp[], double *value) 
    {
      double& omega = p[0];
      double& t = p[1];
      double& phase = p[2];
      double& R = p[3];
      double& y = *value;

      y = cos(omega*t+phase)*exp(-R*t);
      
      if (dp_flg[0]) dp[0] = -t*sin(omega*t+phase)*exp(-R*t);
      if (dp_flg[1]) dp[1] = -omega*sin(omega*t+phase)*exp(-R*t)-R*y;
      if (dp_flg[2]) dp[2] = -sin(omega*t+phase)*exp(-R*t);
      if (dp_flg[3]) dp[3] = -t*y;
      
      return (0);
    }
  
  // Nomalized Lorentzian.
  int norm_lore(void *fs, double p[], int dp_flg[], 
		double dp[], double *value) 
    {
      double& f = p[0];    // frequency variable 
      double& f0 = p[1];   // frequency
      double& w = p[2];    // width
      double& ph = p[3];     // phase
      double& y = *value;

      /* First we calculate some temporary values for calculating
	 s(f,r,ph) and the derrivatives of s. */
      double d_omega = 2.0*pi*(f0-f);
      double r2 = w*pi;
      double cos_ph = cos(ph*pi/180.0);
      double sin_ph = sin(ph*pi/180.0);

      double s_a = 1.0/(1.0+pow((d_omega/r2),2.0));
      double s_b = s_a*d_omega/r2;

      // Then we calculate the point value y.
      y = (cos_ph/r2)*s_a - (sin_ph/r2)*s_b;

      // Derivatives.

      // Frequency variable.
      if (dp_flg[0])
	dp[0] = -2.0*pi* ( -(cos_ph/r2)*(2.0*d_omega/(pow(r2,2)))
			   /(pow((1.0+pow((d_omega/r2),2)),2))
			   -(sin_ph/r2)*(1.0/r2-(pow((d_omega/r2),2))/r2)
			   /(pow((1.0+pow((d_omega/r2),2)),2)) );
      // Frequency.
      if (dp_flg[1])
	dp[1] = 2.0*pi* ( -(cos_ph/r2)*(2.0*d_omega/(pow(r2,2)))
			  /(pow((1.0+pow((d_omega/r2),2)),2))
			  -(sin_ph/r2)*(1.0/r2-(pow((d_omega/r2),2))/r2)
			  /(pow((1.0+pow((d_omega/r2),2)),2)) );

      // Line width, w = r2/pi.
      if (dp_flg[2])
	dp[2] = pi*( -cos_ph*(1.0-pow((d_omega/r2),2.0))
		     /pow((r2+pow(d_omega,2)/r2),2)
		     +sin_ph*(2.0*d_omega*r2)
		     /pow((pow(r2,2)+pow(d_omega,2)),2) );

      // Phase.
      if (dp_flg[3])
	dp[3] = (pi/180.0) * ( -(sin_ph/r2)*s_a - (cos_ph/r2)*s_b );
      
      return (0);
    }

  // Nomalized Gaussian-Lorentzian.
  int norm_gausslore(void *fs, double p[], int dp_flg[], 
		     double dp[], double *value) 
    {
      double& f = p[0];    // frequency variable 
      double& f0 = p[1];   // frequency
      double& w = p[2];    // width
      double& cg = p[3];   // Gaussian fraction
      double& y = *value;

      // The gaussian is only considered within +-gaussrange*w and 
      // is taken as zero outside this interval (se evaluations below).
      const double gaussrange = 4.0;

      /* First we calculate some temporary values for calculating
	 s(f,f0,w,cg) and the derrivatives of s. */
      double g1, g2, g;
      if (fabs(f0-f)<gaussrange*fabs(w))
	{  
	  g1 = (2.0/w)*sqrt(log2/pi);
	  g2 = exp(-log2*pow(2.0*(f0-f)/w,2));
	  g = g1*g2;	  
	}
      else
	g = g1 = g2 = 0.0;
	  
      double l1 = 2.0/(pi*w);
      double l2 = 1.0/(1+pow(2.0*(f0-f)/w,2));
      double l = l1*l2;
      
      //Then we calculate the point value s(f).
      y = cg*g + (1.0-cg)*l;

      // Derivatives.

      // Frequency variable. dp[0] = -dp[1].
      if (dp_flg[0])
	{
	  double dg_df0;
	  if (fabs(f0-f)<gaussrange*fabs(w))
	    dg_df0 = -log2*pow(2.0/w,2)*2.0*(f0-f)*g;
	  else
	    dg_df0 = 0.0;
      
	  double dl_df0 = -l1*pow(l2,2)*pow(2.0/w,2)*2*(f0-f);	  
	  dp[0] = -(cg*dg_df0 + (1.0-cg)*dl_df0);
	}
      
      // Frequency.
      if (dp_flg[1])
	{
	  double dg_df0;
	  if (fabs(f0-f)<gaussrange*fabs(w))
	    dg_df0 = -log2*pow(2.0/w,2)*2.0*(f0-f)*g1*g2;
	  else
	    dg_df0 = 0.0;
      
	  double dl_df0 = -l1*pow(l2,2)*pow(2.0/w,2)*2*(f0-f);	  
	  dp[1] = cg*dg_df0 + (1.0-cg)*dl_df0;
	}
      
      // Line width.
      if (dp_flg[2])
	{
	  double dg_dw;
	  if (fabs(f0-f)<gaussrange*fabs(w))
	    {
	      double dg1_dw = -g1/w;
	      double dg2_dw = log2*pow(f0-f,2)*g2*8.0/pow(w,3);
	      dg_dw = g1*dg2_dw + dg1_dw*g2;
	    }
	  else
	    dg_dw = 0.0;

	  double dl1_dw = -l1/w;
	  double dl2_dw = 8.0*pow(l2,2)*pow((f0-f)/w,2)/w;
	  double dl_dw = l1*dl2_dw + dl1_dw*l2;
	  dp[2] = cg*dg_dw + (1.0-cg)*dl_dw;
	}
      

      // Gaussian fraction.
      if (dp_flg[3])
	dp[3] = g-l;
      
      return (0);
    }

  /* Nomalized two-site exchange lineshape covering from slow exchange
     through intermediate exchange to fast exchange. The expression is
     taken from Lu-Yun Lian and Gordon C. K. Roberts: Effects of
     chemical exchange on NMR spectra, in NMR of Macromolecules, a
     practical approach, edited by G. C. K. Roberts, IRL Press,
     pp. 159 (1995). Note that there is a typo in the formula in the
     book. */

  int norm_exch_lshape(void *fs, double p[], int dp_flg[], 
		       double dp[], double *value) 
    {
      double& f = p[0];    // frequency variable 
      double& tau = p[1];  // exchange time constant
      double& pa = p[2];   // population in site a
      double& fa = p[3];   // frequency in site a
      double& R2a = p[4];  // R2 rate in site a
      double& fb = p[5];   // frequency in site b
      double& R2b = p[6];  // R2 rate in site b
      double& y = *value;

      double pb = 1.0-pa;  // population in site b
      dcomplex alpha_a(R2a+pb/tau, 2.0*pi*(fa-f));
      dcomplex alpha_b(R2b+pa/tau, 2.0*pi*(fb-f));
      dcomplex G = dcomplex(0.0, -1.0)
	*(2.0*pa*pb*tau + pow(tau,2)*(pa*alpha_b+pb*alpha_a))
	/(pa*pb - pow(tau,2)*alpha_a*alpha_b);
      
      // Return value
      y = G.imag();
      
      // No derivatives are supported.
      for(int i=0;i<7;i++)
	{
	  if (dp_flg[i]) return(256*i+1);
	}
      
      return (0);
    }

  // Nomalized discrete Lorentzian.
  // sim==0 not implemented.
  // offset not implemented.
  int norm_dlore(void *fs, double p[], int dp_flg[], 
		double dp[], double *value) 
    {
      enum {F, F0, W, PH, OFFSET, TACQ, SW, SIM, IM_SIGN, IMAG, BLVL, NPARAM};
      double& f = p[F];    // frequency variable 
      double& f0 = p[F0];   // frequency
      double& w = p[W];    // width
      double& ph = p[PH];     // phase
      // double& offset = p[OFFSET]; // frequency offset
      double& tacq = p[TACQ]; // acquisition time
      double& sw = p[SW]; // spectral width
      double& sim = p[SIM]; // flag for simultaneous quadrature.
      double& im_sign = p[IM_SIGN]; // Sign of imaginary part.
      double& imag = p[IMAG]; // if imag!=0.0, return imaginary part.
      double& blvl = p[BLVL]; // base plane level.
      double& y = *value;

      // Derrived values.
      double imag_sign;
      if (im_sign<0.0) imag_sign = -1.0;
      else imag_sign = 1.0;

      // Branch out according to sim.
      if (sim!=0.0)
	{
	  // Simultaneous quadrature.

	  // The expression of the discrete Fourier transform of a complex
	  // exponentially damped sinusoid.

	  // First we calculate some temporary values for calculating
	  // s(f0,w,p) and the derrivatives of s.

	  // Function of ph.
	  dcomplex q(exp(dcomplex(0.0,pi*ph/180.0)));

	  // Functions of w.
	  dcomplex a(exp(dcomplex(-w*pi*tacq,0.0)));
	  dcomplex b(exp(dcomplex(-w*pi/sw,0.0)));

	  // Functions of f0.
	  dcomplex g(exp(dcomplex(0.0, imag_sign*2.0*pi*(f0-f)*tacq)));
	  dcomplex h(exp(dcomplex(0.0, imag_sign*2.0*pi*(f0-f)/sw)));

	  // Calculate the point value.
	  dcomplex cy = q*((dcomplex(1.0)-a*g)
                           /(dcomplex(1.0)-b*h)
                           -dcomplex(blvl))
            /dcomplex(sw,0.0);

	  // return either real or imaginary part.
	  if (int(imag)) y = cy.imag();
	  else y = cy.real();

	  //	  std::cout << q << a << b << g << h << cy << "\n";      

	  // Calculate partial derrivatives of help functions.
	  dcomplex dq(dcomplex(0.0, pi/180.0)*q);
	  dcomplex da(dcomplex(-pi*tacq)*a);
	  dcomplex db(dcomplex(-pi/sw)*b);
	  dcomplex dg(dcomplex(0.0, imag_sign*2.0*pi*tacq)*g);
	  dcomplex dh(dcomplex(0.0, imag_sign*2.0*pi/sw)*h);

	  // Derrivatives.
	  if (dp_flg[F])
	    {
	      // Gives the negative of dp[F0] (see below).
	      dcomplex ds_df0((q/sw) *
			     ((-a*dg)*(1.0-b*h)-(-b*dh)*(1.0-a*g))/
			     pow(1.0-b*h, 2.0));
	      if (int(imag)) dp[F] = -ds_df0.imag();
	      else dp[F] = -ds_df0.real();
	    }
	  if (dp_flg[F0])
	    {
	      dcomplex ds_df0((q/sw) *
			     ((-a*dg)*(1.0-b*h)-(-b*dh)*(1.0-a*g))/
			     pow(1.0-b*h, 2.0));
	      if (int(imag)) dp[F0] = ds_df0.imag();
	      else dp[F0] = ds_df0.real();
	    }
	  if (dp_flg[W])
	    {
	      dcomplex ds_dw((q/sw) *
			     ((-da*g)*(1.0-b*h)-(-db*h)*(1.0-a*g))/
			     pow(1.0-b*h, 2.0));
	      if (int(imag)) dp[W] = ds_dw.imag();
	      else dp[W] = ds_dw.real();
	    }
	  if (dp_flg[PH])
	    {
	      dcomplex ds_dph((dq/sw) * ((1.0-a*g)/(1.0-b*h)-blvl));
	      if (int(imag)) dp[PH] = ds_dph.imag();
	      else dp[PH] = ds_dph.real();
	    }
          if (dp_flg[OFFSET]) return(256*OFFSET+1);
          if (dp_flg[TACQ]) return(256*TACQ+1);
          if (dp_flg[SW]) return(256*SW+1);
          if (dp_flg[SIM]) return(256*SIM+1);
          if (dp_flg[IM_SIGN]) return(256*IM_SIGN+1);
          if (dp_flg[IMAG]) return(256*IMAG+1);
          if (dp_flg[BLVL]) return(256*BLVL+1);
	}
      else
	{
	  // Invalid/unimplemented SIM value.
	  return(256*SIM+2);
	}
      
      return (0);
    }
  

  // Titration curve for single ionizable group.
  int titrate1(void *fs, double p[], int dp_flg[],
	       double dp[], double *value) 
    {
      enum {F_A, F_B, PH, PKA, NPARAM};
      double& f_a = p[F_A];
      double& f_b = p[F_B];
      double& pH = p[PH];
      double& pKa = p[PKA];
      double& y = *value;
      
      // H+ and Ka dependence.
      double z = pow(10.0,(pH-pKa));

      // Calculate the expression for the acid fraction x_a.
      double x_a = 1.0/(z+1);

      // Finally the value.
      y = f_a*x_a + f_b*(1.0-x_a);
      
      // todo: calculate derivatives.
      if (dp_flg[F_A]) dp[F_A] = x_a;
      if (dp_flg[F_B]) dp[F_B] = 1.0-x_a;
      if (dp_flg[PH]) return(256*PH+1);
      if (dp_flg[PKA]) return(256*PKA+1);
      
      return (0);
    }



  // f(p1) = c
  int constant_value(void *fs, double p[], int dp_flg[],
		     double dp[], double *value) 
    {
      double& c = p[0];
      double& y = *value;
      
      y = c;
      
      if (dp_flg[0]) dp[0] = 1.0;
      
      return (0);
    }



  // f(p1) = a*x+b
  int linear(void *fs, double p[], int dp_flg[], double dp[], double *value) 
    {
      double& x = p[0];
      double& a = p[1];
      double& b = p[2];
      double& y = *value;
      
      y = a*x+b;
      
      if (dp_flg[0]) dp[0] = a;
      if (dp_flg[1]) dp[1] = x;
      if (dp_flg[2]) dp[2] = 1.0;
      
      return (0);
    }



  // Here we declare function types for the functions above.
  void declare_cftypes_base(Fuda *fuda_ptr, std::string tag)
  {
    std::string name, descr;
    std::vector<std::string> p_name, p_descr;
    std::vector<unsigned int> p_var;
    Ftype *ft;
    Fuda& fuda = *fuda_ptr;

    // Check tag.
    if (tag!="base" && tag!="all") return;
    
    // Declare function types one by one.
      
    // constant
    ft = fuda.ftype_add_cfunc("constant", 1, constant_value, 0, 0);
    ft->set_descr("Constant function: c");
    ft->set_param(0,1,1,"c","constant value");
    
    // linear
    ft = fuda.ftype_add_cfunc("linear", 3, linear, 0, 0);
    ft->set_descr("Linear function: a*x+b");
    ft->set_param(0,1,1,"a","slope");
    ft->set_param(1,1,1,"x","x-axis variable");
    ft->set_param(2,1,1,"b","y-interception");
    
    // norm_poly
    ft = fuda.ftype_add_cfunc("norm_poly", 2, norm_poly, 0, 0);
    ft->set_descr("Polynomial term: pow(x,p)");
    ft->set_param(0,1,1,"x","x-axis variable");
    ft->set_param(1,1,1,"p","power of polynomial term");
    
    // Natural exponential function.
    ft = fuda.ftype_add_cfunc("norm_exp", 2, norm_exp, 0, 0);
    ft->set_descr("Exponential term: exp(a*x)");
    ft->set_param(0,1,1,"a","exponential factor");
    ft->set_param(1,1,1,"x","exponential variable");
    
    // General exponential function
    ft = fuda.ftype_add_cfunc("norm_pow", 4, norm_pow, 0, 0);
    ft->set_descr("Exponential term: pow(a,b*c*x)");
    ft->set_param(0,1,1,"a","exponential base number");
    ft->set_param(1,1,1,"b","exponential factor b");
    ft->set_param(2,1,1,"c","exponential factor c");
    ft->set_param(3,1,1,"x","exponential variable x");
    
    // Exponential decay.
    ft = fuda.ftype_add_cfunc("norm_exp_decay", 2, norm_exp_decay, 0, 0);
    ft->set_descr("Exponential term: exp(-k*t)");
    ft->set_param(0,1,1,"k","exponential decay constant");
    ft->set_param(1,1,1,"t","time variable");
    
    // Cosine function.
    ft = fuda.ftype_add_cfunc("norm_cos", 3, norm_cos, 0, 0);
    ft->set_descr("Cosine: cos(omega*t+phase)");
    ft->set_param(0,1,1,"omega","angular frequency");
    ft->set_param(1,1,1,"t","time variable");
    ft->set_param(2,1,1,"phase","phase");
    
    // Lorentz.
    ft = fuda.ftype_add_cfunc("norm_lore", 4, norm_lore, 0, 0);
    ft->set_descr("Lorentzian line shape");
    ft->set_param(0,1,1,"f","frequency variable");
    ft->set_param(1,1,1,"f0","center frequency");
    ft->set_param(2,1,1,"w","linewidth");
    ft->set_param(3,1,1,"p","phase");
    
    // Gaussian-Lorentzian line shape.
    ft = fuda.ftype_add_cfunc("norm_gausslore", 4, norm_gausslore, 0, 0);
    ft->set_descr("Gaussian-Lorentzian line shape");
    ft->set_param(0,1,1,"f","frequency variable");
    ft->set_param(1,1,1,"f0","center frequency");
    ft->set_param(2,1,1,"w","linewidth");
    ft->set_param(3,1,1,"cg","gaussian fraction");
    
    // two-site exchange line shape.
    ft = fuda.ftype_add_cfunc("norm_exch_lshape", 7, norm_exch_lshape, 0, 0);
    ft->set_descr("two-site exchange line shape");
    ft->set_param(0,1,1,"f","frequency variable");
    ft->set_param(1,1,1,"tau","time constant");
    ft->set_param(2,1,1,"pa","population for site a");
    ft->set_param(3,1,1,"fa","center frequency for site a");
    ft->set_param(4,1,1,"R2a","R2 for site a");
    ft->set_param(5,1,1,"fb","center frequency for site b");
    ft->set_param(6,1,1,"R2b","R2 for site b");
    
    // Discrete Lorentzian line shape.
    ft = fuda.ftype_add_cfunc("norm_dlore", 11, norm_dlore, 0, 0);
    ft->set_descr("Discrete Lorentzian line shape");
    ft->set_param(0,1,1,"f","frequency variable");
    ft->set_param(1,1,1,"f0","center frequency");
    ft->set_param(2,1,1,"w","linewidth");
    ft->set_param(3,1,1,"p","phase");
    ft->set_param(4,0,0,"offset","frequency offset");
    ft->set_param(5,0,0,"t_acq","acquisition time");
    ft->set_param(6,0,0,"sw","spectral width");
    ft->set_param(7,0,0,"sim","simultaneous quadrature flag");
    ft->set_param(8,0,0,"im_sign","sign of imaginary part");
    ft->set_param(9,0,0,"imag","imaginary flag");
    ft->set_param(10,0,0,"b_lvl","baseline level");

    // linear
    ft = fuda.ftype_add_cfunc("titrate1", 4, titrate1, 0, 0);
    ft->set_descr("titration curve: f_a*x_a(pH,pKa) + f_b*x_b(pH,pKa)");
    ft->set_param(0,1,0,"f_a","acid limiting value");
    ft->set_param(1,1,0,"f_b","basic limiting value");
    ft->set_param(2,1,0,"pH","pH value");
    ft->set_param(3,1,0,"pKa","-log to acid constant");
    
  }
}






