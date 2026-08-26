
namespace FUDA
{
  double z_distrib_crit(double alpha);
  double chi2_distrib_crit(double alpha, int nf);
  double f_distrib_crit(double alpha, int f1, int f2);
  double z_distrib_p(double z);
  double chi2_distrib_p(double x, int nf);
  double f_distrib_p(double f, int f1, int f2);

  namespace perlman
  {    
    double poz(double z);
    double critz(double z);
    double pochisq(double x, int df);
    double critchi(double p, int df);
    double pof(double f, int f1, int f2);
    double critf(double p, int f1, int f2);
  }
}

