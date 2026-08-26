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

#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdio>
#include "fuda_classes.H"
#include "fuda_utils.H"

/* Ftype class member functions */

/* Constructors */
Ftype::Ftype(Fuda *a_fuda, std::string a_name)
  : fuda(a_fuda),
    name(a_name),
    nparam(0),
    p_name(0),
    p_descr(0),
    p_var(0),
    var_index(0) {}


FtypeCFUNC::FtypeCFUNC(Fuda *a_fuda, std::string a_name,
		       unsigned int a_nparam,
		       FUDA::CFUNC a_func, void *a_fs, unsigned int a_fs_size)
  : Ftype(a_fuda, a_name),
    func(a_func),
    fs(a_fs),
    fs_size(a_fs_size)
{
  // Dimension the ftype.
  set_nparam(a_nparam);
}


FtypePROD::FtypePROD(Fuda *a_fuda, std::string a_name,
		     bool a_scale_flg, std::vector<Ftype*>& a_ftypes)
  : Ftype(a_fuda, a_name),
    scale_flg(a_scale_flg),
    ftypes(a_ftypes)
{
  // Dimension normalization vector.
  norm.resize(ftypes.size());

  // Get number of parameters from ftypes.
  int np = 0;

  // This one is for the scaling factor.
  if (scale_flg) np++;

  // Now, add parameter count from ftypes.
  for (unsigned int i=0; i<get_nftype(); i++)
    np += get_ftype(i)->get_nparam();
  
  // Set the number of parameters (dimension the base class).
  set_nparam(np);

  // Setup p_var vector.
  unsigned int p_offset = 0;

  // The scale is the first variable.
  if (scale_flg)
    {
      // Set the p_* entries.
      set_param(p_offset, 1, 1, "scale", "Scale factor");
      p_offset++;
    }

  // Loop over parameters of all ftypes in product to setup p_var.
  for (unsigned int ift=0; ift<ftypes.size(); ift++)
    {
      // Loop over parameters of this ftype.
      for (unsigned int ip=0; ip<ftypes[ift]->get_nparam(); ip++)
	{
	  std::string p_name;
	  std::string p_descr;
	  char add_str[64];

	  // If more than one ftype in product, we add '__%d'.
	  if (ftypes.size()>1) snprintf(add_str, sizeof(add_str),"__%d",ift);
	  else strcpy(add_str,"");
	  
	  // Get the p_name from the contributing ftype and modify.
	  ftypes[ift]->get_p_name(ip, p_name);
	  p_name += add_str;

	  // Get the descr string from the contributing ftype and modify.
	  ftypes[ift]->get_p_descr(ip, p_descr);
	  p_descr += add_str;

	  // Set the p_* entries.
	  set_param(p_offset+ip, ftypes[ift]->is_p_var(ip),
		    ftypes[ift]->is_p_deriv(ip), p_name, p_descr);
	}
      
      // Add up p_offset.
      p_offset += ftypes[ift]->get_nparam();
    }

  // Set the ftype descr.
  std::string descr;
  if (a_scale_flg)
    descr = "Scaled product: scale";
  else
    descr = "Product: ";

  std::string ftype_name;
  for (unsigned int ift=0; ift<ftypes.size(); ift++)
    {
      ftypes[ift]->get_name(ftype_name);
      if (a_scale_flg || ift>0) descr += " * ";
      descr += ftype_name;
    }
  set_descr(descr);
}


FtypeCOMP::FtypeCOMP(Fuda *a_fuda,
		     std::string a_name,
		     bool a_scale_flg, 
		     unsigned int a_p_index,
		     Ftype *a_f,
		     Ftype *a_g)
  : Ftype(a_fuda, a_name),
    scale_flg(a_scale_flg),
    f(a_f),
    g(a_g),
    p_index(a_p_index)
{
  // Get number of parameters from ftypes.
  unsigned int np_f = f->get_nparam();
  unsigned int np_g = g->get_nparam();

  // Check validity of p_index.
  if (p_index<0 || p_index>=np_f)
    throw Uferr::ParamIndexInvalid(p_index);

  // Check that the p_intex parameter in f() is a variable parameter.
  if (!f->is_p_var(p_index))
    throw Uferr::ParamNotVar(p_index);

  // Set number of parameters for ftype.
  int np = 0;

  // Add one for the scaling factor.
  if (scale_flg) np++;

  // Add the paraemter of f() and g() (subtract 1 as we replace one in f().
  np += np_f - 1 + np_g;

  // Set the number of parameters (dimension the base class).
  set_nparam(np);

  // Setup parameter entries for new ftype.

  // Some help variables.
  std::string p_name, p_descr;
  int p_offset = 0;

  // The scale is the first variable.
  if (scale_flg)
    {
      set_param(p_offset, 1, 1, "scale", "scale factor");
      p_offset++;
    }

  // Loop over parameters in f() before p_index.
  for (unsigned int ip=0; ip<p_index; ip++)
    {      
      f->get_p_name(ip,p_name);
      p_name += "__f";
      f->get_p_descr(ip,p_descr);
      p_descr += "__f";
      set_param(ip+p_offset,
		f->is_p_var(ip),
		f->is_p_deriv(ip),
		p_name, p_descr);
    }

  // Add up p_offset.
  if (p_index>=0) p_offset += p_index;

  // Loop over parameters in g().
  for (unsigned int ip=0; ip<np_g; ip++)
    {
      g->get_p_name(ip,p_name);
      p_name += "__g";
      g->get_p_descr(ip,p_descr);
      p_descr += "__g";
      set_param(ip+p_offset,
		g->is_p_var(ip),
		g->is_p_deriv(ip),
		p_name, p_descr);
    }
  
  // Add up p_offset.
  p_offset += np_g;

  // Loop over parameters in f() after p_index.
  for (unsigned int ip=p_index+1; ip<np_f; ip++)
    {
      f->get_p_name(ip,p_name);
      p_name += "__f";
      f->get_p_descr(ip,p_descr);
      p_descr += "__f";
      set_param(ip+p_offset-(p_index+1),
		f->is_p_var(ip),
		f->is_p_deriv(ip),
		p_name, p_descr);
    }

  // Setup and allocate f and g function call workspaces.
  fcall = new FuncCall();
  fcall->set_size(np_f);
  gcall = new FuncCall();
  gcall->set_size(np_g);

  // Set the ftype descr.
  std::string descr;
  if (a_scale_flg)
    descr = "Scaled composite function: scale * ";
  else
    descr = "Composite function: ";

  std::string ftype_name;
  f->get_name(ftype_name);
  descr += ftype_name;
  descr += "( ";
  g->get_name(ftype_name);
  for (unsigned int ip=0; ip<p_index; ip++) descr += ", ";  
  descr += ftype_name;
  for (unsigned int ip=p_index; ip<f->get_nparam()-1; ip++) descr += ", ";  
  descr += ")";
  set_descr(descr);
}


FtypeSUM::FtypeSUM(Fuda *a_fuda, std::string a_name,
		   bool a_scale_flg, std::vector<Ftype*>& a_ftypes)
  : Ftype(a_fuda, a_name),
    scale_flg(a_scale_flg),
    ftypes(a_ftypes)
{
  // Get number of parameters from ftypes.
  int np = 0;

  // This one is for the scaling factor.
  if (scale_flg) np++;

  // Now, add parameter count from ftypes.
  for (unsigned int i=0; i<get_nftype(); i++)
    np += get_ftype(i)->get_nparam();
  
  // Set the number of parameters (dimension the base class).
  set_nparam(np);

  // Setup p_var vector.
  unsigned int p_offset = 0;

  // The scale is the first variable.
  if (scale_flg)
    {
      // Set the p_* entries.
      set_param(p_offset, 1, 1, "scale", "Scale factor");
      p_offset++;
    }

  // Loop over parameters of all ftypes in product to setup p_var.
  for (unsigned int ift=0; ift<ftypes.size(); ift++)
    {
      // Loop over parameters of this ftype.
      for (unsigned int ip=0; ip<ftypes[ift]->get_nparam(); ip++)
	{
	  std::string p_name;
	  std::string p_descr;
	  char add_str[64];

	  // If more than one ftype in product, we add '__%d'.
	  if (ftypes.size()>1) snprintf(add_str, sizeof(add_str),"__%d",ift);
	  else strcpy(add_str,"");
	  
	  // Get the p_name from the contributing ftype and modify.
	  ftypes[ift]->get_p_name(ip, p_name);
	  p_name += add_str;

	  // Get the descr string from the contributing ftype and modify.
	  ftypes[ift]->get_p_descr(ip, p_descr);
	  p_descr += add_str;

	  // Set the p_* entries.
	  set_param(p_offset+ip, ftypes[ift]->is_p_var(ip),
		    ftypes[ift]->is_p_deriv(ip), p_name, p_descr);
	}
      
      // Add up p_offset.
      p_offset += ftypes[ift]->get_nparam();
    }

  // Set the ftype descr.
  std::string descr;
  if (a_scale_flg)
    descr = "Scaled sum: scale * ( ";
  else
    descr = "Sum: ";

  std::string ftype_name;
  for (unsigned int ift=0; ift<ftypes.size(); ift++)
    {
      ftypes[ift]->get_name(ftype_name);
      if (ift>0) descr += " + ";
      descr += ftype_name;
    }
  if (a_scale_flg) descr += " )";

  set_descr(descr);
}


Ftype::~Ftype()
{
  std::cout << "Ftype::~Ftype() not implemented\n";
  //this->print();
  exit(0);  
}


FtypeCFUNC::~FtypeCFUNC()
{
  std::cout << "FtypeCFUNC::~Ftype() not implemented\n";
  exit(0);  
}


FtypeSUM::~FtypeSUM()
{
  std::cout << "FtypeSUM::~Ftype() not implemented\n";
  exit(0);  
}


FtypePROD::~FtypePROD()
{
  std::cout << "FtypePROD::~Ftype() not implemented\n";
  exit(0);  
}


FtypeCOMP::~FtypeCOMP()
{
  std::cout << "FtypeCOMP::~Ftype() not implemented\n";
  exit(0);  
}


Fuda *Ftype::get_fuda() 
{
  return (fuda);
}


void Ftype::set_nparam(unsigned int a_nparam) 
{
  // Dimension the Ftype for nparam parameters.
  nparam = a_nparam;
  p_name.resize(nparam);
  p_descr.resize(nparam);
  p_deriv.resize(nparam);
  p_var.resize(nparam);

  // The number of variables is set to the number of parameters.
  nvar = nparam;
  var_index.resize(nvar);  

  // Initialize the vectors.
  for (unsigned int i=0; i<nparam; i++)
    {
      // We fill p_name and p_descr with "p%d" and "p_descr%d", where
      // "%d" is the number of the parameter.
      char str[64];
      snprintf(str, sizeof(str), "p%d", i);
      p_name[i] = str;
      snprintf(str, sizeof(str), "p_descr%d", i);
      p_descr[i] = str;
      
      // Set all parameters to be variables by default.
      p_var[i] = 1;
      var_index[i] = i;
    }
}


void Ftype::get_name(std::string& a_name)
{
  a_name = name;
}


unsigned int Ftype::get_nparam()
{
  return (nparam);
}

unsigned int Ftype::get_nvar() 
{
  return (nvar); 
}

void Ftype::get_p_name(unsigned int i, std::string& p_nm)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_nm = p_name[i];
}

void Ftype::set_p_name(unsigned int i, std::string p_nm)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_name[i] = p_nm;
}

void Ftype::get_p_descr(unsigned int i, std::string& p_ds)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_ds = p_descr[i];
}

void Ftype::set_p_descr(unsigned int i, std::string p_ds)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_descr[i] = p_ds;
}

void Ftype::set_p_deriv(unsigned int i, bool deriv_flg)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_deriv[i] = deriv_flg;
}

void Ftype::set_p_var(unsigned int i, bool var_flg)
{
  if (i<0 || i>=nparam) throw Uferr::ParamIndexInvalid(i);
  p_var[i] = var_flg;
  set_var_index();
}

void Ftype::set_param(unsigned int i, bool var_flg, bool deriv_flg,
		      std::string a_name, std::string a_descr)
{
  // Set p_var, p_deriv, p_name and p_descr.
  set_p_var (i, var_flg);
  set_p_deriv (i, deriv_flg);
  set_p_name (i, a_name);
  set_p_descr (i, a_descr);
}

bool Ftype::is_p_deriv(unsigned int i)
{
  if (i>=nparam) throw Uferr::ParamIndexInvalid(i);
  return (p_deriv[i]);
}

bool Ftype::is_p_var(unsigned int i)
{
  if (i>=nparam) throw Uferr::ParamIndexInvalid(i);
  return(p_var[i]);
}

unsigned int Ftype::get_var_index(unsigned int i)
{
  if (i>=nvar) throw Uferr::VarIndexInvalid(i);
  return(var_index[i]);
}

void Ftype::set_var_index()
{
  /* This is a privat function which sets up var_index to reflect the
     current p_var vector */

  // Get the new nvar.
  nvar = 0;
  for (unsigned int i=0; i<nparam; i++)
    if (p_var[i]) nvar++;

  // Resize var_index.
  var_index.resize(nvar);

  //if (new_nvar<1 || new_nvar>nparam)
  //  throw Uferr::VarIndexSizeInvalid(new_nvar);

  // Setup new var_index.
  unsigned int var_count = 0;
  for (unsigned int i=0; i<nparam; i++)
    if (p_var[i])
      {
        var_index[var_count] = i;
        var_count++;
      }
}


void Ftype::get_descr(std::string& a_descr)
{
  a_descr=descr;
}


void Ftype::set_descr(std::string a_descr)
{
  descr=a_descr;
}


FUDA::CFUNC FtypeCFUNC::get_func() 
{
  return (func);
}


void *FtypeCFUNC::get_fs()
{
  return (fs);
}


unsigned int FtypeCFUNC::get_fs_size()
{
  return (fs_size);
}


FUDA::Ftype_type FtypeCFUNC::get_type()
{
  return (FUDA::FTYPE_CFUNC);
}


FUDA::Ftype_type FtypePROD::get_type()
{
  return (FUDA::FTYPE_PROD);
}


FUDA::Ftype_type FtypeCOMP::get_type()
{
  return (FUDA::FTYPE_COMP);
}


FUDA::Ftype_type FtypeSUM::get_type()
{
  return (FUDA::FTYPE_SUM);
}


// Make direct call to underlying Ftype function.
void FtypeCFUNC::call(double p[], int dp_flg[], double dp[], double *value)
{
  // Supply the func structure and call the function.
  int rtnval = func(fs, p, dp_flg, dp, value);
  if (rtnval != 0)
    {
      std::string msg = "FtypeCFUNC call returned non-zero value for ftype ";
      std::string name;
      get_name(name);
      msg+=name+"\n";
      throw Uferr::FtypeCallError(msg);      
    }
}


// Calculate ftype function as a function of all parameters.
double FtypeCFUNC::call(std::vector<double>& p_vec)
{
  // Vector size must match number of parameters in function.
  if (p_vec.size()!=get_nparam())
    throw Uferr::NumParamInvalid(p_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(get_fuda()->fcall_get_ref());

  // Setup func_call structure.

  // Set  values.
  for (unsigned int i=0; i<get_nparam(); i++)
    func_call.p[i] = p_vec[i];

  // Setup function call.
  double value;
  
  // Call the function.
  int rtnval = func(fs, func_call.p,
		    func_call.dp_flg_false,
		    func_call.dp, &value);

  // Check return value.
  if (rtnval != 0)
    {
      std::string msg = "FtypeCFUNC function returned non-zero value\n";
      throw Uferr::FtypeCallError(msg);      
    }

  return (value);
}

/* Print object contents */
void Ftype::print()
{
  std::cout << "name  : " << name << "\n";
  std::cout << "descr : " << descr << "\n";
  if (get_type()==FUDA::FTYPE_CFUNC)
    std::cout << "type  : CFUNC\n";
  else if (get_type()==FUDA::FTYPE_PROD)
    std::cout << "type  : PROD\n";
  else if (get_type()==FUDA::FTYPE_COMP)
    std::cout << "type  : COMP\n";
  else if (get_type()==FUDA::FTYPE_SUM)
    std::cout << "type  : SUM\n";
  else
    std::cout << "type  : unknown\n";
  std::cout << "nparam: " << nparam << "\n";
  for(unsigned int i=0; i<p_name.size(); i++)
    std::cout << "p_name[" << i << "] : " << p_name[i] << "\n";
  for(unsigned int i=0; i<p_descr.size(); i++)
    std::cout << "p_descr[" << i << "] : " << p_descr[i] << "\n";
  for(unsigned int i=0; i<p_var.size(); i++)
    std::cout << "p_var[" << i << "] : " << p_var[i] << "\n";
  for(unsigned int i=0; i<p_deriv.size(); i++)
    std::cout << "p_deriv[" << i << "] : " << p_deriv[i] << "\n";
  for(unsigned int i=0; i<var_index.size(); i++)
    std::cout << "var_index[" << i << "] : " << var_index[i] << "\n";
}

/* Print object contents */
void FtypeCFUNC::print()
{
  Ftype::print();
  std::cout << "fs_size: " << fs_size << "\n";
  std::cout << "\n";
}

/* Print object contents */
void FtypePROD::print()
{
  Ftype::print();
  for (unsigned int i=0; i<get_nftype(); i++)
    {
      std::string name;
      get_ftype(i)->get_name(name);
      printf("Ftype[%d] : %s\n", i, name.c_str());
    }  
  std::cout << "\n";
}


/* Print object contents */
void FtypeCOMP::print()
{
  Ftype::print();
  std::string fname, gname;
  get_f()->get_name(fname);
  get_g()->get_name(gname);
  printf("f(g()) ftypes : %s(%s())\n", 
	 fname.c_str(), gname.c_str());
  printf("g() replaces parameter %d in f()\n", p_index);

  printf("\n");
}


/* Print object contents */
void FtypeSUM::print()
{
  Ftype::print();
  for (unsigned int i=0; i<get_nftype(); i++)
    {
      std::string name;
      get_ftype(i)->get_name(name);
      printf("Ftype[%d] : %s\n", i, name.c_str());
    }  
  std::cout << "\n";
}


/* Get number of ftypes in product */
unsigned int FtypePROD::get_nftype()
{
  return(ftypes.size());
}


/* Get number of ftypes in sum */
unsigned int FtypeSUM::get_nftype()
{
  return(ftypes.size());
}


/* Get reference to ftype in product */
Ftype *FtypePROD::get_ftype(unsigned int i)
{
  if (i<0 || i>=ftypes.size())
    std::exit(0);
  return(ftypes[i]);
}


/* Get reference to ftype in sum */
Ftype *FtypeSUM::get_ftype(unsigned int i)
{
  if (i<0 || i>=ftypes.size())
    std::exit(0);
  return(ftypes[i]);
}


void FtypePROD::call(double p[], int dp_flg[], double dp[], double *value)
{
  // Get scale.
  double scale;
  if (scale_flg) scale = p[0];
  else scale = 1.0;

  /* Initialize normalized return value. Later we multiply with scale
     to get final return value */
  double f_norm = 1.0;

  // Initialize normalization factor for all ftypes in product.
  for (unsigned int ift=0; ift<get_nftype(); ift++)
    norm[ift] = 1.0;
  
  /* Initialize parameter offset to 0. */
  unsigned int p_offset = 0;
  if (scale_flg) p_offset++;

  // Loop over ftypes.
  for (unsigned int ift=0; ift<get_nftype(); ift++)
    {
      // Calculate function values for each dimension.
      double dvalue;
      ftypes[ift]->call(&p[p_offset],
			&dp_flg[p_offset],
			&dp[p_offset],
			&dvalue);
      
      // Update f_norm.
      f_norm *= dvalue;
      
      // Update normalization factors.
      for (unsigned int jft=0; jft<get_nftype(); jft++)
	if (jft!=ift) norm[jft] *= dvalue; 
      
      // Count up offset.
      p_offset += ftypes[ift]->get_nparam();
    }
  
  // Calculate the function value.
  *value = scale*f_norm;
  
  // We now turn to calculating the derrivatives.
  p_offset = 0;
  
  // First the scale. 
  if (scale_flg)
    {
      // dy/dscale is simply f_norm.
      if (dp_flg[0]) dp[0] = f_norm;
      p_offset++;
    }
  
  /* The rest have already been calculated by the normalized
     function, but needs scaling with the normalization constants
     for the respective dimensions. */
  for (unsigned int ift=0; ift<get_nftype(); ift++)
    {
      for (unsigned int iparam=0; iparam<ftypes[ift]->get_nparam(); iparam++)
	if (dp_flg[p_offset+iparam])
	  dp[p_offset+iparam] *= scale*norm[ift];
      
      // Count up parameter offset.
      p_offset += ftypes[ift]->get_nparam();
    }

  // Here goes the debug spam.
  if (0)
    {
      for(unsigned int i=0; i<get_nparam(); i++)
	{
	  std::cout << "prod: " << i << " " 
		    << p[i] << " ";
	  
	  if (dp_flg[i])
	    std::cout << " deriv " << dp[i];
	  
	  std::cout << "\n";
	}
    }
  
}


double FtypePROD::call(std::vector<double>& p_vec)
{
  // Vector size must match number of parameters in function.
  if (p_vec.size()!=get_nparam())
    throw Uferr::NumParamInvalid(p_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(get_fuda()->fcall_get_ref());

  // Setup func_call structure.

  // Set  values.
  for (unsigned int i=0; i<get_nparam(); i++)
    func_call.p[i] = p_vec[i];

  // Setup function call.
  double value;

  // Call the function.
  call(func_call.p,
       func_call.dp_flg_false,
       func_call.dp, &value);

  return (value);
}


void FtypeSUM::call(double p[], int dp_flg[], double dp[], double *value)
{
  // Get scale.
  double scale;
  if (scale_flg) scale = p[0];
  else scale = 1.0;

  /* Initialize normalized return value. Later we multiply with scale
     to get final return value */
  double f_norm = 0.0;

  /* Initialize parameter offset to 0. */
  unsigned int p_offset = 0;
  if (scale_flg) p_offset++;

  // Loop over ftypes.
  for (unsigned int ift=0; ift<get_nftype(); ift++)
    {
      // Calculate function values for each dimension.
      double dvalue;
      ftypes[ift]->call(&p[p_offset],
			&dp_flg[p_offset],
			&dp[p_offset],
			&dvalue);
      
      // Update f_norm.
      f_norm += dvalue;
      
      // Count up offset.
      p_offset += ftypes[ift]->get_nparam();
    }
  
  // Calculate the function value.
  *value = scale*f_norm;
  
  // We now turn to calculating the derrivatives.
  p_offset = 0;
  
  // First the scale. 
  if (scale_flg)
    {
      // dy/dscale is simply f_norm.
      if (dp_flg[0]) dp[0] = f_norm;
      p_offset++;
    }
  
  /* The rest have already been calculated by the normalized
     function, but eventually needs scaling. */
  if (scale_flg)
    {
      for (unsigned int ift=0; ift<get_nftype(); ift++)
	{
	  for (unsigned int iparam=0; iparam<ftypes[ift]->get_nparam();
	       iparam++)
	    if (dp_flg[p_offset+iparam])
	      dp[p_offset+iparam] *= scale;
	  
	  // Count up parameter offset.
	  p_offset += ftypes[ift]->get_nparam();
	}
    }
}


double FtypeSUM::call(std::vector<double>& p_vec)
{
  // Vector size must match number of parameters in function.
  if (p_vec.size()!=get_nparam())
    throw Uferr::NumParamInvalid(p_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(get_fuda()->fcall_get_ref());

  // Setup func_call structure.

  // Set  values.
  for (unsigned int i=0; i<get_nparam(); i++)
    func_call.p[i] = p_vec[i];

  // Setup function call.
  double value;

  // Call the function.
  call(func_call.p,
       func_call.dp_flg_false,
       func_call.dp, &value);

  return (value);
}


/* Get reference to f ftype */
Ftype *FtypeCOMP::get_f()
{
  return(f);
}


/* Get reference to g ftype */
Ftype *FtypeCOMP::get_g()
{
  return(g);
}


/* Get scale_flg */
bool FtypePROD::is_scaled()
{
  return(scale_flg);
}


/* Get scale_flg */
bool FtypeSUM::is_scaled()
{
  return(scale_flg);
}


/* Get scale_flg */
bool FtypeCOMP::is_scaled()
{
  return(scale_flg);
}


/* Get p_index */
unsigned int FtypeCOMP::get_p_index()
{
  return(p_index);
}


void FtypeCOMP::call(double p[], int dp_flg[], double dp[], double *value)
{
  // Get/set scale and scale parameter offset (used below).
  double scale = 1.0;
  unsigned int sp_offset = 0;
  if (scale_flg) 
    {
      scale = p[0];
      sp_offset = 1;
    }

  // Get number of parameters in f and g.
  unsigned int np_f = f->get_nparam();
  unsigned int np_g = g->get_nparam();

  // Calculate index in p[] where first f parameters begin.
  unsigned int if1 = sp_offset;
  
  // Calculate index in p[] where g parameters begin.
  unsigned int ig1 = if1 + p_index;

  // Copy g function parameters to gcall workspace and set g_dp_flg to
  // true if one or more g parmeters require derivatives.
  bool g_dp_flg = 0;
  for (unsigned int ip=0; ip<np_g; ip++)
    {
      gcall->p[ip] = p[ig1+ip];
      gcall->dp_flg[ip] = dp_flg[ig1+ip];
      g_dp_flg = g_dp_flg or gcall->dp_flg[ip];
    }

  // Evaluate g function.
  double gvalue;
  g->call(gcall->p, gcall->dp_flg, gcall->dp, &gvalue);

  // Copy f function parameters before g to fcall workspace.
  for (unsigned int ip=0; ip<p_index; ip++)
    {
      fcall->p[ip] = p[if1+ip];
      fcall->dp_flg[ip] = dp_flg[if1+ip];
    }

  // Copy gvalue and g_dp_flg to parameter p_index in f() workspace.
  fcall->p[p_index] = gvalue;
  fcall->dp_flg[p_index] = g_dp_flg;

  // Copy f function parameters after g to fcall workspace.
  unsigned int p_offset = sp_offset + np_g -1;
  for (unsigned int ip=p_index+1; ip<np_f; ip++)
    {
      fcall->p[ip] = p[ip+p_offset];
      fcall->dp_flg[ip] = dp_flg[ip+p_offset];
    }

  // Evaluate f function.
  double fvalue;
  f->call(fcall->p, fcall->dp_flg, fcall->dp, &fvalue);

  // Get df/dg value.
  double df_dg = fcall->dp[p_index];

  // Set return function value.
  *value = scale * fvalue;

  // Setup dp[] derivative array.
  p_offset = sp_offset;

  // Set scale derivative.
  if (scale_flg && dp_flg[0]) dp[0] = fvalue;

  // Move derivatives back for first f() parameters to dp[].
  for (unsigned int ip=0; ip<p_index; ip++)
    if (fcall->dp_flg[ip])
      dp[ip+p_offset] = scale * fcall->dp[ip];

  // Move derivatives back for g() parameters to dp[]. We scale with
  // df/dg according to the chain rule: dy/dx = dy/du * du/dx.
  p_offset += p_index;
  for (unsigned int ip=0; ip<np_g; ip++)
    if (gcall->dp_flg[ip])
      dp[ip+p_offset] = scale * df_dg * gcall->dp[ip];

  // Move derivatives back for last f() parameters to dp[].
  p_offset += np_g - (p_index + 1);
  for (unsigned int ip=p_index+1; ip<np_f; ip++)
    if (fcall->dp_flg[ip])
      dp[ip+p_offset] = scale * fcall->dp[ip];

  // Historical derivative debug dump removed: it was permanently disabled by `&& 0`.

  
  if (0)
    {
      
      for(unsigned int i=0; i<np_f; i++)
	{
	  std::cout << "f: " << i << " " 
		    << fcall->p[i] << " ";
	  
	  if (fcall->dp_flg[i])
	    std::cout << " deriv " << fcall->dp[i];

	  std::cout << "\n";
	}
      for(unsigned int i=0; i<np_g; i++)
	{
	  std::cout << "g: " << i << " " 
		    << gcall->p[i] << " ";
	  
	  if (gcall->dp_flg[i])
	    std::cout << " deriv " << gcall->dp[i];

	  std::cout << "\n";
	}
      for(unsigned int i=0; i<get_nparam(); i++)
	{
	  std::cout << "fog: " << i << " " 
		    << p[i] << " ";
	  
	  if (dp_flg[i])
	    std::cout << " deriv " << dp[i];

	  std::cout << "\n";
	}
    }
  

}

  
double FtypeCOMP::call(std::vector<double>& p_vec)
{
  // Vector size must match number of parameters in function.
  if (p_vec.size()!=get_nparam())
    throw Uferr::NumParamInvalid(p_vec.size());

  // Get reference to fcall record.
  FuncCall& func_call = *(get_fuda()->fcall_get_ref());

  // Setup func_call structure.

  // Set  values.
  for (unsigned int i=0; i<get_nparam(); i++)
    func_call.p[i] = p_vec[i];

  // Setup function call.
  double value;

  // Call the function.
  call(func_call.p,
       func_call.dp_flg_false,
       func_call.dp, &value);

  return (value);
}


