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

/* User-defined functions in fuda */
namespace FUDA
{
  // second order polynomial a*x*x+b*x+c.
  int user_example(void *fs, double p[], int dp_flg[],
		   double dp[], double *value) 
    {
      double& x = p[0];
      double& a = p[1];
      double& b = p[2];
      double& c = p[3];
      double& y = *value;

      y = a*x*x+b*x+c;
      
      if (dp_flg[0]) dp[0] = 2.0*a*x+b;
      if (dp_flg[1]) dp[1] = x*x;
      if (dp_flg[2]) dp[2] = x;
      if (dp_flg[3]) dp[3] = 1.0;
      
      return (0);
    }
  

  // Here we declare function types for the functions above.
  void declare_cftypes_user(Fuda *fuda_ptr, std::string tag)
  {
    std::string name, descr;
    std::vector<std::string> p_name, p_descr;
    std::vector<unsigned int> p_var;
    Ftype *ft;
    Fuda& fuda = *fuda_ptr;

    // Check tag.
    if (tag!="user" && tag!="all") return;
    
    // Register function types one by one.
      
    // Second order polynomial.
    ft = fuda.ftype_add_cfunc("user_example", 4, user_example, 0, 0);
    ft->set_descr("Second-order polynomial: y = a*x*x+b*x+c");
    ft->set_param(0,1,1,"x","explanatory variable");
    ft->set_param(1,1,1,"a","2. order coefficient");
    ft->set_param(2,1,1,"b","1. order coefficient");
    ft->set_param(3,1,1,"c","0. order coefficient");    
  }
}






