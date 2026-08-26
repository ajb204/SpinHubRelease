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



/* Ptype class member functions */

/* Constructor */
Ptype::Ptype(Fuda *a_fuda, std::string& a_name, FUDA::Param_kind a_kind)
{
  fuda = a_fuda;
  name = a_name;
  kind = a_kind;
  init = 0.0;

  /* If this is a dependent parameter, we initialize the relevant
     entries. */
  if (kind==FUDA::PKIND_PARAM)
    {
      free_flg = 1;
      bounds = 0;
      lower = -1.5;
      upper = 1.5;
      delta = 1.0;
      esd = 0.1;
    }
}


// Destructor.
Ptype::~Ptype()
{
  // Empty.
}


/* Object data get/set functions */
Fuda *Ptype::get_fuda() 
{
  return (fuda);
}

void Ptype::get_name(std::string& a_name) { a_name = name; }

FUDA::Param_kind Ptype::get_kind() { return (kind); }

bool Ptype::is_free() { return (free_flg); }

void Ptype::set_free(bool a_free) { free_flg=a_free; }

void Ptype::free() { free_flg=1; }

void Ptype::fix() { free_flg=0; }

double Ptype::get_init() { return (init); }

void Ptype::set_init(double a_init) { init=a_init; }

double Ptype::get_value() { return (init); }

void Ptype::set_value(double a_init) { init=a_init; }

unsigned int Ptype::get_bounds() { return (bounds); }

void Ptype::set_bounds(unsigned int a_bounds) { bounds=a_bounds; }

double Ptype::get_lower() { return (lower); }

void Ptype::set_lower(double a_lower) { lower=a_lower; }

double Ptype::get_upper() { return (upper); }

void Ptype::set_upper(double a_upper) { upper=a_upper; }

double Ptype::get_delta() { return (delta); }

void Ptype::set_delta(double a_delta) { delta=a_delta; }

double Ptype::get_esd() { return (esd); }

void Ptype::set_esd(double a_esd) { esd=a_esd; }


bool Ptype::is_referenced()
{
  // Return true if the ptype is referenced by at least one Param
  // object.

  bool referenced = 0;
  
  // Loop over all params.
  for(Param_iterator pi = fuda->param_begin(); pi != fuda->param_end(); pi++)
    {
      // Compare reference to self refererence.
      if (this == (*pi)->get_ptype())
	{
	  referenced = 1;
	  break;
	}
    }
  return (referenced);
}


/* Print object contents */
void Ptype::print()
{
  std::cout << "name    : " << name << "\n";
  std::cout << "kind    : " << kind << "\n";
  std::cout << "init    : " << init << "\n";
  std::cout << "free_flg: " << free_flg << "\n";
  std::cout << "bounds  : " << bounds << "\n";
  std::cout << "lower   : " << lower << "\n";
  std::cout << "upper   : " << upper << "\n";
  std::cout << "delta   : " << delta << "\n";
  std::cout << "\n";
}
