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



/* Param class member functions */

/* Construct with name and ptype */
Param::Param(Fuda *a_fuda,
	     std::string& a_name,
	     Ptype *a_ptype)
{
  // Initialize fundamental data entries.
  fuda = a_fuda;
  name = a_name;
  ptype = a_ptype;
  kind = ptype->get_kind();

  // Zero eval struct.
  eval.active = 0;
  eval.index = 0;

  // Branch out according to kind.
  if (kind==FUDA::PKIND_EXPL)
    ;
  else if (kind==FUDA::PKIND_PARAM)
    {
      // Initialize dependent parameter specific entries.
      init = ptype->get_init();
      value = init;
      free_flg = ptype->is_free();
      bounds = ptype->get_bounds();
      lower = ptype->get_lower();
      upper = ptype->get_upper();
      delta = ptype->get_delta();
      esd = ptype->get_esd();
    }
  else if (kind==FUDA::PKIND_CONST)
    // Initialize constant value.
    value = ptype->get_init();
  else throw Uferr::ParamKindInvalid();
}


// Construct with name and kind.
Param::Param(Fuda *a_fuda,
	     std::string& a_name,
	     FUDA::Param_kind a_kind)
{
  // Initialize fundamental data entries.
  fuda = a_fuda;
  name = a_name;
  kind =a_kind;
  ptype = NULL;

  // Zero eval struct.
  eval.active = 0;
  eval.index = 0;

  // Zero entries.
  init = 0.0;
  value = 0.0;
  free_flg = 0;
  bounds = 0;
  lower = 0.0;
  upper = 0.0;
  delta = 0.0;
  esd = 0.0;

  // Branch out according to kind.
  if (kind==FUDA::PKIND_EXPL)
    ;
  else if (kind==FUDA::PKIND_PARAM)
    {
      // Set the parameter free by default.
      free_flg = 1;
    }
  else if (kind==FUDA::PKIND_CONST)
    ;
  else throw Uferr::ParamKindInvalid();
}


// Destructor.
Param::~Param()
{
  // Empty.
}


Fuda *Param::get_fuda() 
{
  return (fuda);
}

FUDA::Param_kind Param::get_kind()
{
  return (kind);
}

void Param::get_name(std::string& a_name)
{
  a_name = name;
}

void Param::get_ptype_name(std::string& a_name)
{
  if (ptype) ptype->get_name(a_name);
  else a_name = "";
}

Ptype *Param::get_ptype()
{
  return (ptype);
}

bool Param::is_free()
{
  return (free_flg);
}

void Param::set_free(bool a_free)
{
  free_flg=a_free;
  fuda->param_clear_sync();
}

void Param::free()
{
  free_flg=1;
  fuda->param_clear_sync();
}

void Param::fix() 
{
  free_flg=0;
  fuda->param_clear_sync();
}

bool Param::is_referenced()
{
  // Return true if the parameter is referenced by at least one Func
  // object.

  bool referenced = 0;
  
  // Loop over all functions.
  for(Func_iterator fi = fuda->func_begin(); fi != fuda->func_end(); fi++)
    {
      // Get number of parameters for function.
      int nparam = (*fi)->get_ftype()->get_nparam();

      // Loop over all parameters of function.
      for(int ip = 0; ip<nparam; ip++)
	{
	  // Compare reference to self refererence.
	  if (this == (*fi)->get_param(ip))
	    {
	      referenced = 1;
	      break;
	    }
	  if (referenced) break;
	}
    }
  return (referenced);
}

double Param::get_init()
{
  return (init);
}

void Param::set_init(double a_init)
{
  init=a_init;
  fuda->param_clear_sync();
}

double Param::get_value()
{
  return (value);
}

void Param::set_value(double a_value)
{
  value=init=a_value; 
  fuda->param_clear_sync();
}

void Param::init_value()
{ 
  value = init;
  fuda->param_clear_sync();
}

unsigned int Param::get_bounds()
{
  return (bounds);
}

void Param::set_bounds(unsigned int a_bounds) 
{
  bounds=a_bounds; 
  fuda->param_clear_sync();
}

double Param::get_lower()
{
  return (lower);
}

void Param::set_lower(double a_lower)
{
  lower=a_lower;
  fuda->param_clear_sync();
}

double Param::get_upper()
{
  return (upper);
}

void Param::set_upper(double a_upper)
{
  upper=a_upper; 
  fuda->param_clear_sync();
}

double Param::get_delta()
{
  return (delta);
}

void Param::set_delta(double a_delta)
{
  delta=a_delta; 
  fuda->param_clear_sync();
}

double Param::get_esd() 
{
  return (esd);
}

void Param::set_esd(double a_esd)
{
  esd=a_esd;
  fuda->param_clear_sync();
}

bool Param::is_eval_active()
{
  return (eval.active);
}

unsigned int Param::get_eval_index()
{
  return (eval.index);
}



/* Print object contents */
void Param::print()
{
  std::cout << "name  : " << name << "\n";
  if (kind==FUDA::PKIND_CONST)
    std::cout << "kind  : CONST\n";
  else if (kind==FUDA::PKIND_PARAM)
    std::cout << "kind  : PARAM\n";
  else if (kind==FUDA::PKIND_EXPL)
    std::cout << "kind  : EXPL\n";
  else
    std::cout << "kind  unknown : " << kind << "\n";
  std::string str;
  if (ptype) ptype->get_name(str);
  else str = "No ptype";
  std::cout << "ptype : " << str << "\n";

  if (kind==FUDA::PKIND_EXPL)
    ;
  else if (kind==FUDA::PKIND_PARAM)
    {
      std::cout << "value   : " << value << "\n";
      std::cout << "init    : " << init << "\n";
      std::cout << "free_flg: " << free_flg << "\n";
      std::cout << "bounds  : " << bounds << "\n";
      std::cout << "lower   : " << lower << "\n";
      std::cout << "upper   : " << upper << "\n";
      std::cout << "delta   : " << delta << "\n";
      std::cout << "esd     : " << esd << "\n";
    }
  else if (kind==FUDA::PKIND_CONST)
    std::cout << "value : " << value << "\n";
  else throw Uferr::ParamKindInvalid();
  
  std::cout << "eval.active : " << eval.active << "\n";
  std::cout << "eval.index  : " << eval.index << "\n";
  std::cout << "\n";
}
