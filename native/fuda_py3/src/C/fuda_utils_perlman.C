
#include <cstdlib>
#include <ctime>
#include <math.h>
#include "fuda_utils_perlman.H"

///////////////////////////////////////////////////////////////
// FUDA::perlman namespace utility functions.
///////////////////////////////////////////////////////////////

/* All functions in the FUDA::perlman namespace is taken from the
   NETLIB a/perlman statistics package. They have been slightly
   altered to conform to c++ and placed within a namespace.

   Functions outside in the FUDA namespace are wrappers for some of
   the Perlman functions + other related utility functions. */



namespace FUDA
{

  // Return normal alpha critical value.
  double z_distrib_crit(double alpha)
  {
    return (perlman::critz(alpha));
  }
  

  // Return accumulated propability from -oo to z.
  double z_distrib_p(double z)
  {
    return (perlman::poz(z));
  }
  

  // Return chi2 alpha critical value for n degrees of freedom.
  double chi2_distrib_crit(double alpha, int n)
  {
    return (perlman::critchi(alpha,n));
  }
  

  // Return accumulated propability from 0 to x with n degrees of freedom.
  double chi2_distrib_p(double x, int n)
  {
    return (perlman::pochisq(x,n));
  }
  

  // Return F alpha critical value for f1, f2 degrees of freedom.
  double f_distrib_crit(double alpha, int f1, int f2)
  {
    return (perlman::critf(alpha,f1,f2));
  }
  

  // Return accumulated probability from 0 to f with f1 and f2 degrees
  // of freedom.
  double f_distrib_p(double f, int f1, int f2)
  {
    return (perlman::pof(f,f1,f2));
  }
  



  namespace perlman
  {
    

/*HEADER
	Module:       z.c
	Purpose:      compute approximations to normal z distribution probabilities
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DZTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

/*LINTLIBRARY*/

#define	Z_EPSILON      0.000001       /* accuracy of critz approximation */
#define	Z_MAX          6.0            /* maximum meaningful z value */


#ifdef	ZTEST
int main ()
	{
	double	z;
	printf ("%4s  %10s  %10s  %10s\n",
		"z", "poz(z)", "poz(-z)", "z'");
	for (z = 0.0; z <= Z_MAX; z += .01)
		{
		printf ("%4.2f  %10.6f  %10.6f  %10.6f\n",
			z, poz (z), poz (-z), critz (poz (z)));
		}
	}
#endif	//ZTEST

/*FUNCTION poz: probability of normal z value */
/*ALGORITHM
	Adapted from a polynomial approximation in:
		Ibbetson D, Algorithm 209
		Collected Algorithms of the CACM 1963 p. 616
	Note:
		This routine has six digit accuracy, so it is only useful for absolute
		z values < 6.  For z values >= to 6.0, poz() returns 0.0.
*/

/* Returns cumulative probability from -oo to z */
double poz (double z)
	{
	double	y, x, w;
	
	if (z == 0.0)
		x = 0.0;
	else
		{
		y = 0.5 * fabs (z);
		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0)
			{
			w = y*y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w +0.005198775019) * w
				-0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
				-0.531923007300) * w +0.797884560593) * y * 2.0;
			}
		else
			{
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
				+0.000535310849) * y +0.999936657524;
			}
		}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
	}

/*FUNCTION critz: compute critical z value to produce given probability */
/*ALGORITHM
	Begin with upper and lower limits for z values (maxz and minz)
	set to extremes.  Choose a z value (zval) between the extremes.
	Compute the probability of the z value.  Set minz or maxz, based
	on whether the probability is less than or greater than the
	desired p.  Continue adjusting the extremes until they are
	within Z_EPSILON of each other.
*/

/*VAR returns z such that fabs (poz(p) - z) <= .000001 */
double critz (double p) /* p: critical probability level */
	{
	double	minz = -Z_MAX;    /* minimum of range of z */
	double	maxz = Z_MAX;     /* maximum of range of z */
	double	zval = 0.0;       /* computed/returned z value */
	double	pval;     /* prob (z) function, pval := poz (zval) */
	
	if (p <= 0.0 || p >= 1.0)
		return (0.0);
	
	while (maxz - minz > Z_EPSILON)
		{
		pval = poz (zval);
		if (pval > p)
			maxz = zval;
		else
			minz = zval;
		zval = (maxz + minz) * 0.5;
		}
	return (zval);
	}




/*
	Module:       chisq.c
	Purpose:      compute approximations to chisquare distribution probabilities
	Contents:     pochisq(), critchi()
	Uses:         poz() in z.c (Algorithm 209)
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DCHISQTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

/*LINTLIBRARY*/

#define	CHI_EPSILON     0.000001    /* accuracy of critchi approximation */
#define	CHI_MAX     99999.0         /* maximum chi square value */

#define	LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define	BIGX           20.0         /* max value to represent exp (x) */
#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))


#ifdef	CHISQTEST
double	Prob[] = { .10, .05, .01, .005, .001, -1.0 };
main ()
	{
	int 	df;
	int 	p;
	printf ("%-4s ", "df");
	for (p = 0; Prob[p] > 0.0; p++)
		printf ("%8.3f ", Prob[p]);
	putchar ('\n');
	for (df = 1; df < 30; df++)
		{
		printf ("%4d ", df);
		for (p = 0; Prob[p] > 0.0; p++)
			printf ("%8.3f ", critchi (Prob[p], df));
		putchar ('\n');
		}
	}
#endif	/* CHISQTEST */

/*FUNCTION pochisq: probability of chi sqaure value */
/*ALGORITHM Compute probability of chi square value.
	Adapted from:
		Hill, I. D. and Pike, M. C.  Algorithm 299
		Collected Algorithms for the CACM 1967 p. 243
	Updated for rounding errors based on remark in
		ACM TOMS June 1985, page 185
*/
double pochisq (double x, int df)
/* x: obtained chi-square value */
/* df: degrees of freedom */
	{
	double y = 0.0;
	double	a, s;
	double	e, c, z;
	int 	even;     /* true if df is an even number */
	
	if (x <= 0.0 || df < 1)
		return (1.0);
	
	a = 0.5 * x;
	even = (2*(df/2)) == df;
	if (df > 1)
		y = ex (-a);
	s = (even ? y : (2.0 * poz (-sqrt (x))));
	if (df > 2)
		{
		x = 0.5 * (df - 1.0);
		z = (even ? 1.0 : 0.5);
		if (a > BIGX)
			{
			e = (even ? 0.0 : LOG_SQRT_PI);
			c = log (a);
			while (z <= x)
				{
				e = log (z) + e;
				s += ex (c*z-a-e);
				z += 1.0;
				}
			return (s);
			}
		else
			{
			e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
			c = 0.0;
			while (z <= x)
				{
				e = e * (a / z);
				c = c + e;
				z += 1.0;
				}
			return (c * y + s);
			}
		}
	else
		return (s);
	}

/*FUNCTION critchi: compute critical chi square value to produce given p */
double critchi (double p, int df)
	{
	double	minchisq = 0.0;
	double	maxchisq = CHI_MAX;
	double	chisqval;
	
	if (p <= 0.0)
		return (maxchisq);
	else if (p >= 1.0)
		return (0.0);
	
	chisqval = df / sqrt (p);    /* fair first value */
	while (maxchisq - minchisq > CHI_EPSILON)
		{
		if (pochisq (chisqval, df) < p)
			maxchisq = chisqval;
		else
			minchisq = chisqval;
		chisqval = (maxchisq + minchisq) * 0.5;
		}
	return (chisqval);
	}



/*
	Module:       f.c
	Purpose:      compute approximations to F distribution probabilities
	Contents:     pof(), critf()
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DFTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

/*LINTLIBRARY*/

#ifndef	I_PI        /* 1 / pi */
#define	I_PI        0.3183098861837906715377675
#endif
#define	F_EPSILON     0.000001       /* accuracy of critf approximation */
#define	F_MAX      9999.0            /* maximum F ratio */

#ifdef	FTEST

int 	DFs[] = { 1, 2, 5, 10, 20, 40, 60, 120, -1 };
double	Prob[] = { .10, .05, .01, .005, .001, -1.0 };

main ()
	{
	int 	dfnum;
	int 	dfdenom;
	int 	p;
	
	for (p = 0; Prob[p] > 0.0; p++)
		{
		printf ("alpha = %g                      df1\n", Prob[p]);
		printf ("%-4s\\", "df2");
		for (dfnum = 0; DFs[dfnum] > 0; dfnum++)
			printf ("%7d ", DFs[dfnum]);
		putchar ('\n');
		for (dfdenom = 0; DFs[dfdenom] > 0; dfdenom++)
			{
			printf ("%4d ", DFs[dfdenom]);
			for (dfnum = 0; DFs[dfnum] > 0; dfnum++)
				printf ("%7.2f ", critf (Prob[p], DFs[dfnum], DFs[dfdenom]));
			putchar ('\n');
			}
		putchar ('\n');
		}
	}
#endif	// FTEST

/*FUNCTION pof: probability of F */
/*ALGORITHM Compute probability of F ratio.
	Adapted from Collected Algorithms of the CACM
	Algorithm 322
	Egon Dorrer
*/
double pof (double F, int df1, int df2)
	{
	int	i, j;
	int	a, b;
	double	w, y, z, d, p;
	
	if (F < F_EPSILON || df1 < 1 || df2 < 1)
		return (1.0);
	a = df1%2 ? 1 : 2;
	b = df2%2 ? 1 : 2;
	w = (F * df1) / df2;
	z = 1.0 / (1.0 + w);
	if (a == 1)
		if (b == 1)
			{
			p = sqrt (w);
			y = I_PI; /* 1 / 3.14159 */
			d = y * z / p;
			p = 2.0 * y * atan (p);
			}
		else
			{
			p = sqrt (w * z);
			d = 0.5 * p * z / w;
			}
	else if (b == 1)
		{
		p = sqrt (z);
		d = 0.5 * z * p;
		p = 1.0 - p;
		}
	else
		{
		d = z * z;
		p = w * z;
		}
	y = 2.0 * w / z;
#ifdef	REMARK /* speedup modification suggested by Tolman (wrong answer!) */
	if (a == 1)
		for (j = b + 2; j <= df2; j += 2)
			{
			d *= (1.0 + a / (j - 2.0)) * z;
			p += d * y / (j - 1.0);
			}
	else
		{
		double	zk = 1.0;
		for (j = (df2 - 1) / 2; j; j--)
			zk *= z;
		d *= zk * df2/b;
		p *= zk + w * z * (zk - 1.0)/(z-1.0);
		}
#else /* original version */
	for (j = b + 2; j <= df2; j += 2)
		{
		d *= (1.0 + a / (j - 2.0)) * z;
		p = (a == 1 ? p + d * y / (j - 1.0) : (p + w) * z);
		}
#endif	// REMARK
	y = w * z;
	z = 2.0 / z;
	b = df2 - 2;
	for (i = a + 2; i <= df1; i += 2)
		{
		j = i + b;
		d *= y * j / (i - 2.0);
		p -= z * d / j;
		}
	/* correction for approximation errors suggested in certification */
	if (p < 0.0)
		p = 0.0;
	else if (p > 1.0)
		p = 1.0;
	return (1.0-p);
	}

/*FUNCTION critf: compute critical F value t produce given probability */
/*ALGORITHM
	Begin with upper and lower limits for F values (maxf and minf)
	set to extremes.  Choose an f value (fval) between the extremes.
	Compute the probability of the f value.  Set minf or maxf, based
	on whether the probability is less than or greater than the
	desired p.  Continue adjusting the extremes until they are
	within F_EPSILON of each other.
*/
double critf (double p, int df1, int df2)
	{
	double	fval;
	double	maxf = F_MAX;     /* maximum F ratio */
	double	minf = 0.0;       /* minimum F ratio */
	
	if (p <= 0.0 || p >= 1.0)
		return (0.0);
	
	fval = 1.0 / p;             /* the smaller the p, the larger the F */
	
	while (fabs (maxf - minf) > F_EPSILON)
		{
		if (pof (fval, df1, df2) < p)     /* F too large */
			maxf = fval;
		else                              /* F too small */
			minf = fval;
		fval = (maxf + minf) * 0.5;
		}
	
	return (fval);
	}




 
  } // perlman
} // FUDA
